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
        MPI_Comm comm	
    );

    int number_of_species();
    const double* molar_mass();

    double mean_specific_heat_at_CP_R(
        double T, 
	double* mole_fractions
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
}

#endif
