#include <vector>
#include <string>
using namespace std;
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

#include <fstream>
string read(string path) {
    ifstream input_file(path);
    return string((istreambuf_iterator<char>(input_file)), istreambuf_iterator<char>());
}

#undef NDEBUG
#include <cassert>
#include <unistd.h>

#include <occa.hpp>
#include "nekrk.h"

static const double R = 1.380649e-23 * 6.02214076e23;

namespace {
    occa::kernel production_rates_kernel, production_rates_fp32_kernel;
        occa::kernel transportCoeffs_kernel;
    occa::kernel number_of_species_kernel;
    occa::kernel mean_specific_heat_at_CP_R_kernel;
    occa::kernel molar_mass_kernel;

    occa::device device;

    double reference_pressure = 1;
    double reference_temperature = 1;
    double reference_mass_rate = 1;
    double reference_energy_rate = 1;

    int n_species = -1;
    double *m_molar;

    MPI_Comm comm;
    int initialized = 0;

        std::vector<std::string> species_names;
};

void set_number_of_species()
{
  auto o_tmp = device.malloc<int>(1);
  number_of_species_kernel(o_tmp);
  o_tmp.copyTo(&n_species);
}

void set_molar_mass()
{
  if(!m_molar) m_molar = new double[n_species];
  auto o_tmp = device.malloc<double>(n_species);
  molar_mass_kernel(o_tmp);
  o_tmp.copyTo(m_molar);
}

void setup(const char* mech, occa::device _device, occa::properties kernel_properties,
       int _group_size, MPI_Comm _comm/*, bool transport*/)
{
    comm   = _comm;
    device = _device;
    int group_size = std::max(_group_size, /*32*/1); // Why ?

    const std::string mechFile = std::string(mech) + ".c";
    if(!getenv("NEKRK_PATH")) {
      std::string path = std::string(getenv("HOME")) + "/.local/nekRK";
      setenv("NEKRK_PATH", path.c_str(),0);
    }
    if (!getenv("OCCA_DIR")) occa::env::OCCA_DIR = std::string(getenv("NEKRK_PATH")) + "/";
    occa::env::OCCA_INSTALL_DIR = occa::env::OCCA_DIR;

    kernel_properties["includes"].asArray();
    kernel_properties["flags"].asObject();

    { // workaround to bypass occa parser
      kernel_properties["compiler_flags"] += " -D__FG_ENABLE__";
      kernel_properties["okl/strict_headers"] = false;
      kernel_properties["includes"] += mechFile;
      const std::string incStatement = " -I " + std::string(getenv("NEKRK_PATH") ?: ".") + "/share/mechanisms";
      kernel_properties["compiler_flags"] += incStatement;

      if(device.mode() == "CUDA") {
        setenv("OCCA_CXXFLAGS", incStatement.c_str(), 1); // required for launcher
                kernel_properties["compiler_flags"] += " --fmad=true";
        kernel_properties["compiler_flags"] += " -D__FG_DEVICE__=__device__";
        kernel_properties["compiler_flags"] += " -D__FG_CONST__=__constant__";
        kernel_properties["compiler_flags"] += " -D__FG_EXP_APPROX__=expf";
        kernel_properties["compiler_flags"] += " -D__FG_POW_APPROX__=powf";
        kernel_properties["compiler_flags"] += " -D__FG_LOG10_APPROX__=log10f";
        kernel_properties["compiler_flags"] += " --use_fast_math";
      } else if(device.mode() == "HIP") {
        setenv("OCCA_CXXFLAGS", incStatement.c_str(), 1); // required for launcher
        kernel_properties["compiler_flags"] += " -O3 ";
        kernel_properties["compiler_flags"] += " -ffp-contract=fast ";
        kernel_properties["compiler_flags"] += " -D__FG_DEVICE__=__device__";
        kernel_properties["compiler_flags"] += " -D__FG_CONST__=__constant__";
        kernel_properties["compiler_flags"] += " -D__FG_EXP_APPROX__=expf";
        kernel_properties["compiler_flags"] += " -D__FG_POW_APPROX__=powf";
        kernel_properties["compiler_flags"] += " -D__FG_LOG10_APPROX__=log10f";
        kernel_properties["compiler_flags"] += " -funsafe-math-optimizations";
        kernel_properties["compiler_flags"] += " -ffast-math";
      } else if(device.mode() == "OPENCL") {
        setenv("OCCA_CXXFLAGS", incStatement.c_str(), 1); // required for launcher
        kernel_properties["compiler_flags"] += " -D__FG_DEVICE__=__device__";
        kernel_properties["compiler_flags"] += " -D__FG_CONST__=__constant__";
        kernel_properties["compiler_flags"] += " -cl-std=CL1.2";
        kernel_properties["compiler_flags"] += " -cl-mad-enable";
        kernel_properties["compiler_flags"] += " -D__FG_EXP_APPROX__=native_exp";
        kernel_properties["compiler_flags"] += " -D__FG_POW_APPROX__=native_pow";
        kernel_properties["compiler_flags"] += " -D__FG_LOG10_APPROX__=native_log10";
        kernel_properties["compiler_flags"] += " -cl-unsafe-math-optimizations";
        kernel_properties["compiler_flags"] += " -cl-fast-relaxed-math";
      } else {
        kernel_properties["compiler_flags"] += " -D__FG_DEVICE__=";
        kernel_properties["compiler_flags"] += " -D__FG_CONST__=const";
        kernel_properties["compiler_flags"] += " -D__FG_EXP_APPROX__=expf";
        kernel_properties["compiler_flags"] += " -D__FG_POW_APPROX__=powf";
        kernel_properties["compiler_flags"] += " -D__FG_LOG10_APPROX__=log10f";

                group_size = 1;
      }
    }

    kernel_properties["defines/p_BLOCKSIZE"] = std::to_string(group_size);

    occa::properties kernel_properties_fp32 = kernel_properties;

        {
            std::string dfloatType = "double";
            kernel_properties["defines/dfloat"] = dfloatType;
            kernel_properties["compiler_flags"] += " -Ddfloat=" + dfloatType;
            kernel_properties["compiler_flags"] += " -D__FG_EXP__=exp";
        }
        {
            std::string dfloatType = "float";
            kernel_properties_fp32["defines/dfloat"] = dfloatType;
            kernel_properties_fp32["compiler_flags"] += " -Ddfloat=" + dfloatType;
            if(device.mode() == "CUDA" || device.mode() == "HIP") {
                kernel_properties_fp32["compiler_flags"] += " -D__FG_EXP__=expf";
            } else if(device.mode() == "HIP") {
                kernel_properties_fp32["compiler_flags"] += " -D__FG_EXP__=expf";
            } else if(device.mode() == "OPENCL") {
                kernel_properties_fp32["compiler_flags"] += " -D__FG_EXP__=native_exp";
            } else {
                kernel_properties_fp32["compiler_flags"] += " -D__FG_EXP__=expf";
            }
        }

    int rank;
    MPI_Comm_rank(comm, &rank);
    for (int r = 0; r < 2; r++) {
      if ((r == 0 && rank == 0) || (r == 1 && rank > 0)) {
        const std::string okl_path = std::string(getenv("NEKRK_PATH") ?: ".")+"/okl/wrapper.okl";
                /*if (verbose)*/ { std::cerr<<"Rates\n"; }
        production_rates_kernel           = device.buildKernel(okl_path.c_str(), "production_rates", kernel_properties);
        production_rates_fp32_kernel      = device.buildKernel(okl_path.c_str(), "production_rates", kernel_properties_fp32);
        number_of_species_kernel          = device.buildKernel(okl_path.c_str(), "number_of_species", kernel_properties);
        mean_specific_heat_at_CP_R_kernel = device.buildKernel(okl_path.c_str(), "mean_specific_heat_at_CP_R", kernel_properties);
        molar_mass_kernel                 = device.buildKernel(okl_path.c_str(), "molar_mass", kernel_properties);
                /*if (transport)*/ {
                    kernel_properties["defines/CFG_FEATURE_TRANSPORT"] = "1";
                    kernel_properties["compiler_flags"] += " -DCFG_FEATURE_TRANSPORT=1";
                    /*if (verbose)*/ { std::cerr<<"Transport\n"; }
                    transportCoeffs_kernel           = device.buildKernel(okl_path.c_str(), "transport", kernel_properties);
                }
                kernel_properties["defines/CFG_FEATURE_TRANSPORT"] = "0";
                kernel_properties["compiler_flags"] += " -DCFG_FEATURE_TRANSPORT=0";
      }
      MPI_Barrier(comm);
    }
    set_number_of_species();
    set_molar_mass();

    if(rank==0) {
      std::cerr << "nekRK initialized successfully\n"
                << "mechanism file: " << mechFile << '\n'
                << "nSpecies: "<< n_species << '\n'
                << "active occa mode: " << device.mode() << '\n'
                << "blockSize: " << group_size << '\n';
    }

    initialized = 1;

        auto mechanism_code_file_path = std::string(getenv("NEKRK_PATH") ?: ".") + "/share/mechanisms/" + mechFile;
        vector<string> lines = split(read(mechanism_code_file_path),"\n");
        if (lines.size() < 2) { cerr<<"Invalid mechanism code file:"<<mechanism_code_file_path; abort(); }
        //assert(is_number(lines[0].substr(2)));
        //number_of_active_species = stoi(lines[0].substr(2));
        species_names = split(lines[1].substr(2)," ");
}

#include "nekrk.h"

/* API */

void nekRK::init(const char* model_path, occa::device device,
      occa::properties kernel_properties, int group_size, MPI_Comm comm/*, bool transport*/)
{
    setup(model_path, device, kernel_properties, group_size, comm/*, bool transport*/);
}

double nekRK::mean_specific_heat_at_CP_R(double T, double* mole_fractions)
{
  assert(initialized);
  auto tmp = new double[1];
  auto o_tmp = device.malloc<double>(1);
  auto o_mole_fractions = device.malloc<double>(n_species, mole_fractions);
  mean_specific_heat_at_CP_R_kernel(T, o_mole_fractions, o_tmp);
  o_tmp.copyTo(tmp);
  return tmp[0];
}

void nekRK::set_reference_parameters(
    double reference_pressure_in,
    double reference_temperature_in,
    double reference_velocity_in,
    double reference_length_in,
    double* reference_mass_fractions_in)
{
  assert(initialized);
  double sum_rcp_molar_mass = 0.;
  for(int k=0;k<number_of_species();k++)
    sum_rcp_molar_mass += reference_mass_fractions_in[k] / m_molar[k];
  double reference_molar_mass = 1./sum_rcp_molar_mass;
  double reference_concentration = reference_pressure_in / R / reference_temperature_in;
  double reference_density = reference_concentration * reference_molar_mass;

  auto reference_mole_fractions = new double[number_of_species()];
  for(int k=0;k<number_of_species();k++)
    reference_mole_fractions[k] = 1./m_molar[k] * reference_molar_mass * reference_mass_fractions_in[k];

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
                    occa::memory o_mass_rates, occa::memory o_energy_rate,
                int fp32)
{
  assert(initialized);
  const double pressure_R = pressure * reference_pressure / R;
  occa::kernel kernel = production_rates_kernel;
  if(fp32) kernel =  production_rates_fp32_kernel;
  kernel(
    n_states,
    pressure_R,
    o_temperature,
    o_mass_fractions,
    o_mass_rates,
    o_energy_rate,
    reference_temperature,
    1./reference_mass_rate,
    R/reference_energy_rate);
}

int nekRK::number_of_species()
{
  assert(initialized);
  return n_species;
}

const double* nekRK::molar_mass()
{
  assert(initialized);
  return (const double*) m_molar;
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

const std::vector<std::string> nekRK::species_names()
{
    return ::species_names;
}
