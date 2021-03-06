#include <vector>
#include <string>
using namespace std;
//#include <regex>
//vector<string> split(const string str, const string regex_str) { return {sregex_token_iterator(str.begin(), str.end(), regex(regex_str), -1), sregex_token_iterator()}; }
vector<string> split (string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}
//bool is_number(const std::string &s) { return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit); }
#include <fstream>
string read(string path) {
    ifstream input_file(path);
    return string((istreambuf_iterator<char>(input_file)), istreambuf_iterator<char>());
}

#include <cassert>
#include <unistd.h>

#include <sys/stat.h>
bool exists(std::string name) { struct stat buffer; return (stat (name.c_str(), &buffer) == 0); }

const double R = 1.380649e-23 * 6.02214076e23;

#include <occa.hpp>
#include "nekrk.h"

namespace {
    occa::kernel production_rates_kernel;
    occa::kernel transportCoeffs_kernel;
    occa::kernel number_of_species_kernel;
    occa::kernel mean_specific_heat_at_CP_R_kernel;
    occa::kernel molar_mass_kernel;

    occa::device device;

    double reference_pressure;
    double reference_temperature;
    double reference_mass_rate;
    double reference_energy_rate;

    std::vector<std::string> species_names;
    std::optional<uint> number_of_active_species;
    std::vector<double> species_molar_mass;

    MPI_Comm comm;

}

void setup(const char* mechanism_code_file_path_ptr, occa::device _device, occa::properties kernel_properties,
       const int group_size, MPI_Comm _comm, bool transport, bool verbose)
{
    assert(mechanism_code_file_path_ptr);
    comm   = _comm;
    device = _device;

    string mechanism_code_file_path = mechanism_code_file_path_ptr;
    if (!exists(mechanism_code_file_path)) mechanism_code_file_path = string(getenv("NEKRK_PATH") ?: ".") + "/mechanisms/" + mechanism_code_file_path + ".c";
    else {
#if __cpp_lib_starts_ends_with
        if(mechanism_code_file_path.ends_with(".yaml")) cerr << "Warning: Mechanism code generation is not integrated in nekRK. First generate the code from a YAML chemical kinetics mechanism using main.py. And then pass the resulting code file path here";
        if(!mechanism_code_file_path.ends_with(".c")) cerr << "Warning: File path given for mechanism code does not end in .c. This path should be used to specify the mechanism code file generated by main.py from a YAML chemical kinetics mechanism";
#endif
    }
    assert(exists(mechanism_code_file_path)); // FIXME: OCCA seems to create an empty file otherwise
    // Remark: Using absolute paths can avoid some reparse as OCCA will reparse files if given with another (relative) path

    vector<string> lines = split(read(mechanism_code_file_path),"\n");
    assert(lines.size() > 2);
    //assert(is_number(lines[0].substr(2)));
    number_of_active_species = stoi(lines[0].substr(2));
    species_names = split(lines[1].substr(2)," ");
    //for(auto&& value: species_names) cout <<"\""<< value << "\"\n";
    auto values = split(lines[2].substr(3), " ");
    for(auto&& value: values) species_molar_mass.push_back(std::stof(value));
    //for(auto&& value: species_molar_mass) cout <<"\""<< value << "\"\n";

    kernel_properties["includes"].asArray();
    kernel_properties["includes"] += mechanism_code_file_path;
    kernel_properties["defines/float"] = "double";
    kernel_properties["defines/fg_exp2"] = "exp2"; // "__expf"
    kernel_properties["defines/p_blockSize"] = to_string(group_size);
    kernel_properties["flags"].asObject(); // ?
    //kernel_properties["compiler_flags"] += " --prec-div=false --prec-sqrt=false";
    //kernel_properties["compiler_flags"] += " --use_fast_math";
    //kernel_properties["okl/enabled"] = false;

    string okl_path = string(getenv("NEKRK_PATH") ?: ".")+"/okl/lib.okl";

    // ?
    int rank;
    MPI_Comm_rank(comm, &rank);
    for (int r = 0; r < 2; r++) {
      if ((r == 0 && rank == 0) || (r == 1 && rank > 0)) {
            if (transport) {
                kernel_properties["defines/CFG_FEATURE_TRANSPORT"] = "1";
                if (verbose) { cerr<<"Transport\n"; }
                transportCoeffs_kernel           = device.buildKernel(okl_path.c_str(), "transport", kernel_properties);
            }
            kernel_properties["defines/CFG_FEATURE_TRANSPORT"] = "0";
            if (verbose) { cerr<<"Rates\n"; }
            production_rates_kernel           = device.buildKernel(okl_path.c_str(), "production_rates", kernel_properties);
            mean_specific_heat_at_CP_R_kernel = device.buildKernel(okl_path.c_str(), "mean_specific_heat_at_CP_R", kernel_properties); // FIXME: Should always be host CPU
      }
      MPI_Barrier(comm);
    }

    if(rank==0 && verbose) { cerr<<"mechanism file: "<<mechanism_code_file_path<<"\n"<<nekRK::number_of_species()<<" species\n"; }
}

#include "nekrk.h"

/* API */

void nekRK::init(const char* model_path, occa::device device,
      occa::properties kernel_properties, int group_size, MPI_Comm comm, bool transport, bool verbose)
{
  setup(model_path, device, kernel_properties, group_size, comm, transport, verbose);
}

double nekRK::mean_specific_heat_at_CP_R(double T, vector<double> mole_fractions)
{
    double mcp[1];
    auto o_mcp = device.malloc<double>(1);
  auto o_mole_fractions = device.malloc<double>(number_of_species(), mole_fractions.data());
    // This is not a kernel, just an interface to fg_molar_heat_capacity_at_constant_pressure_R (for the single reference state) (FIXME: should be on CPU)
    mean_specific_heat_at_CP_R_kernel(T, o_mole_fractions, o_mcp);
    o_mcp.copyTo(mcp);
    return mcp[0];
}

void nekRK::set_reference_parameters(
    double reference_pressure_in,
    double reference_temperature_in,
    double reference_velocity_in,
    double reference_length_in,
    double* reference_mass_fractions_in)
{
    double sum_rcp_molar_mass = 0.;
    for(int k=0;k<number_of_species();k++)
      sum_rcp_molar_mass += reference_mass_fractions_in[k] / species_molar_mass()[k];
    double reference_molar_mass = 1./sum_rcp_molar_mass;
    double reference_concentration = reference_pressure_in / R / reference_temperature_in;
    double reference_density = reference_concentration * reference_molar_mass;

    vector<double> reference_mole_fractions(number_of_species());
    for(int k=0;k<number_of_species();k++)
        reference_mole_fractions[k] = 1./species_molar_mass()[k] * reference_molar_mass * reference_mass_fractions_in[k];

    double reference_molar_heat_capacity_R =
      nekRK::mean_specific_heat_at_CP_R(reference_temperature_in, reference_mole_fractions);

    reference_pressure = reference_pressure_in;
    reference_temperature = reference_temperature_in;
    const double reference_time = reference_length_in / reference_velocity_in;
    reference_mass_rate = reference_density / reference_time;
    reference_energy_rate =
      -(reference_molar_heat_capacity_R * reference_pressure_in) / reference_time;
}

void nekRK::production_rates(const int n_states, double pressure,
                            occa::memory o_temperature, occa::memory o_mass_fractions,
                    occa::memory o_mass_rates, occa::memory o_energy_rate)
{
  const double pressure_R = pressure * reference_pressure / R;
  production_rates_kernel(
       n_states,
       pressure_R,
       o_temperature,
       o_mass_fractions,
       o_mass_rates,
       o_energy_rate,
       reference_temperature,
       1./reference_mass_rate,
       R/reference_energy_rate
    );
}

int nekRK::number_of_species()
{
    assert(species_names);
    return ::species_names.size();
}

uint nekRK::number_of_active_species()
{
    assert(::number_of_active_species != ~0);
    return ::number_of_active_species;
}

const std::vector<double> nekRK::species_molar_mass()
{
  return ::species_molar_mass;
}

const std::vector<std::string> nekRK::species_names()
{
    return ::species_names;
}

void nekRK::transportCoeffs(int nStates, double pressure_Pa, occa::memory T, occa::memory Yi, occa::memory lambda, occa::memory mue, occa::memory rho_Di, double reference_temperature)
{
    transportCoeffs_kernel(
        nStates,
        pressure_Pa,
        T,
        Yi,
        lambda,
        mue,
        rho_Di,
        reference_temperature
    );
}
