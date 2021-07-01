#ifndef NEKRK_H
#define NEKRK_H

#include <occa.hpp>
#include "mpi.h"

namespace nekRK {

    void init(
        const char* mechanism,
        occa::device device,
        occa::properties kernel_properties,
        int blockSize,
        MPI_Comm comm,
        bool transport,
        bool verbose
    );

    int number_of_species();
    int number_of_active_species();
    const std::vector<double> species_molar_mass();
    const std::vector<std::string> species_names();

    double mean_specific_heat_at_CP_R(
        double T,
    vector<double> mole_fractions
    );

    void set_reference_parameters(
        double reference_pressure,
        double reference_temperature,
        double reference_velocity,
        double reference_length,
        double reference_mass_fractions[]
    );

    void production_rates(
        int n_states,
        double pressure,
        occa::memory temperature,
        occa::memory mass_fractions,
        occa::memory mass_rates,
        occa::memory energy_rate
    );

    void transportCoeffs(int nStates, double pressure_Pa, occa::memory T, occa::memory Yi, occa::memory mue, occa::memory lambda, occa::memory rho_Di, const double reference_temperature);

}

#endif
