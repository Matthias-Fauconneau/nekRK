#include <string>
using namespace std;
#include <cassert>
#include <unistd.h>

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

    int n_species = -1;
    double *m_molar = 0;

    MPI_Comm comm;

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
       const int group_size, MPI_Comm _comm, bool transport, bool verbose)
{
    comm   = _comm;
    device = _device;

    // FIXME: Assert exists. OCCA seems to create an empty file otherwise
    std::string mechFile = string(getenv("NEKRK_PATH") ?: ".") + "/share/mechanisms/" + string(mech) + ".c";

    kernel_properties["includes"].asArray();
    kernel_properties["includes"] += mechFile;
    kernel_properties["defines/float"] = "double";
    kernel_properties["defines/fg_exp2"] = "exp2"; // "__expf"
    kernel_properties["defines/p_blockSize"] = to_string(group_size);
    kernel_properties["flags"].asObject(); // ?
    //kernel_properties["compiler_flags"] += " --prec-div=false --prec-sqrt=false";
    //kernel_properties["compiler_flags"] += " --use_fast_math";
    //kernel_properties["okl/enabled"] = false;

    string okl_path = string(getenv("NEKRK_PATH") ?: ".")+"/okl/fuego_wrapper.okl";

    // ?
    int rank;
    MPI_Comm_rank(comm, &rank);
    for (int r = 0; r < 2; r++) {
      if ((r == 0 && rank == 0) || (r == 1 && rank > 0)) {
            if (transport) {
                kernel_properties["defines/CFG_FEATURE_TRANSPORT"] = "1";
                transportCoeffs_kernel           = device.buildKernel(okl_path.c_str(), "transport", kernel_properties);
            }
            kernel_properties["defines/CFG_FEATURE_TRANSPORT"] = "0";
            production_rates_kernel           = device.buildKernel(okl_path.c_str(), "production_rates", kernel_properties);
            number_of_species_kernel          = device.buildKernel(okl_path.c_str(), "number_of_species", kernel_properties);
            mean_specific_heat_at_CP_R_kernel = device.buildKernel(okl_path.c_str(), "mean_specific_heat_at_CP_R", kernel_properties); // FIXME: Should always be host CPU
            molar_mass_kernel                 = device.buildKernel(okl_path.c_str(), "molar_mass", kernel_properties); // FIXME: Should always be host CPU, used once as the interface to get molar_mass
      }
      MPI_Barrier(comm);
    }
    set_number_of_species();
    set_molar_mass();

    if(rank==0 && verbose) {
      std::cout << "nekRK initialized successfully\n";
        std::cout << "mechanism file: "<< mechFile <<"\n";
        std::cout << "nSpecies: "<< n_species <<"\n";
    }
}

#include "nekrk.h"

/* API */

void nekRK::init(const char* model_path, occa::device device,
      occa::properties kernel_properties, int group_size, MPI_Comm comm, bool transport)
{
  setup(model_path, device, kernel_properties, group_size, comm, transport, false);
}

double nekRK::mean_specific_heat_at_CP_R(double T, double* mole_fractions)
{
    auto mcp = new double[1];
    auto o_mcp = device.malloc<double>(1);
  auto o_mole_fractions = device.malloc<double>(n_species, mole_fractions);
    // This is not a kernel, just to interface a single call to fg_molar_heat_capacity_at_constant_pressure_R on CPU
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
  return n_species;
}

const double* nekRK::molar_mass()
{
  return (const double*) m_molar;
}

void nekRK::transportCoeffs(int nStates, double p, occa::memory T, occa::memory Yi, occa::memory mue, occa::memory lambda, occa::memory rho_Di, double reference_temperature)
{
    transportCoeffs_kernel(
        nStates,
        p,
        T,
        Yi,
        mue,
        lambda,
        rho_Di,
        reference_temperature
    );
}
