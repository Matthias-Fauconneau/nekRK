#include <iostream>
#include <fstream>
#include <cstddef>
#include <string>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <unistd.h>
#include <getopt.h>

#include "mpi.h"

#include "occa.hpp"
#include <nekrk.h>
#include "refBK1.h"

using namespace std;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // parse command line options
    int err = 0;
    std::string threadModel;
    int n_states = 1;
    int mode = 0;
    int blockSize = 1;
    int nRep = 0; // repetitions
    int fp32 = 0;
    int ci = 0;
    int deviceId = 0;
    int deviceIdFlag = 0;
    std::string mech("GRIMech-3.0");

    while(1) {
      static struct option long_options[] =
      {
        {"mode", required_argument, 0, 'e'},
        {"backend", required_argument, 0, 'd'},
        {"n-states", required_argument, 0, 'n'},
        {"block-size", required_argument, 0, 'b'},
        {"repetitions", required_argument, 0, 'r'},
        {"fp32", no_argument, 0, 'p'},
        {"ci", no_argument, 0, 'c'},
        {"mechanism-name", required_argument, 0, 'f'},
        {"device-id", required_argument, 0, 'i'},
        {0, 0, 0, 0}
      };
      int option_index = 0;
      int c = getopt_long (argc, argv, "s:", long_options, &option_index);

      if (c == -1)
        break;

      switch(c) {
        case 'e':
          mode = std::stoi(optarg);
          break;
        case 'd':
          threadModel.assign(strdup(optarg));
          break;
        case 'n':
          n_states = std::stoi(optarg)/size;
          break;
        case 'b':
          blockSize = std::stoi(optarg);
          break;
        case 'r':
          nRep = std::stoi(optarg);
          break;
        case 'p':
          fp32 = 1;
          break;
        case 'c':
          ci = 1;
          break;
        case 'f':
          mech.assign(optarg);
          break;
        case 'i':
          deviceId = std::stoi(optarg);
                    deviceIdFlag = 1;
          break;
    default:
          err++;
      }
    }

    if(mode < 1 ||  threadModel.size() < 1) err++;

    if(err > 0) {
        if(rank == 0)
                    fprintf(stderr, "Usage: ./bk --mode 1|2 --backend SERIAL|CUDA|HIP|OPENCL --n-states n "
                 "[--repetitions n] [--fp32] [--ci] [--mechanism-name] [--block-size  n] [--device-id  n]\n");
        exit(EXIT_FAILURE);
    }

    // setup device
    if(!deviceIdFlag) {
      deviceId = 0;
      long int hostId = gethostid();
      long int* hostIds = (long int*) std::calloc(size, sizeof(long int));
      MPI_Allgather(&hostId, 1, MPI_LONG, hostIds,1, MPI_LONG, MPI_COMM_WORLD);
      for (int r = 0; r < rank; r++)
        if (hostIds[r] == hostId) deviceId++;
    }

    char deviceConfig[BUFSIZ];
    const int platformId = 0;
    if(strstr(threadModel.c_str(), "CUDA")) {
        sprintf(deviceConfig, "{mode: 'CUDA', device_id: %d}",deviceId);
    }else if(strstr(threadModel.c_str(),  "HIP")) {
        sprintf(deviceConfig, "{mode: 'HIP', device_id: %d}",deviceId);
    }else if(strstr(threadModel.c_str(),  "OpenCL")) {
        sprintf(deviceConfig, "{mode: 'OpenCL', device_id: %d, platform_id: %d}", deviceId, platformId);
    }else if(strstr(threadModel.c_str(),  "OpenMP")) {
        sprintf(deviceConfig, "{mode: 'OpenMP'}");
    }else {
        sprintf(deviceConfig, "{mode: 'Serial'}");
    }

    occa::device device;
    std::string deviceConfigString(deviceConfig);
    device.setup(deviceConfigString);

    if(rank == 0) {
      std::cerr << "number of states: " << n_states << '\n';
      std::cerr << "use fp32: " << fp32 << '\n';
      std::cerr << "number of repetitions: " << nRep << '\n';
    }

    nekRK::init(mech.c_str(), device, {}, blockSize, MPI_COMM_WORLD/*,*/ /*transport*//*mode==2*/);
    const int n_species = nekRK::number_of_species();

        for (int k=0; k<n_species; k++) { cout << nekRK::species_names()[k]; if (k<n_species-1) { cout << ' '; } }
        cout << '\n';

    // setup reference quantities
    double pressure_Pa   = 101325.;
    double temperature_K = 1000.;
    auto mole_fractions = new double[n_species];
    for (int i=0; i<n_species; i++) mole_fractions[i] = 1./(double)n_species;

    auto molar_mass_species = nekRK::molar_mass();

    double molar_mass = 0.;
    for(int k=0; k<n_species; k++) molar_mass += mole_fractions[k] * molar_mass_species[k];

    auto reference_mass_fractions = new double[n_species];
    for(int k=0; k<n_species; k++) {
      reference_mass_fractions[k] = mole_fractions[k] * molar_mass_species[k]  / molar_mass;
    }
    const double reference_length = 1.;
    const double reference_velocity = 1.;
    const double reference_pressure = pressure_Pa;
    const double reference_temperature = temperature_K;

    nekRK::set_reference_parameters(reference_pressure,
                            reference_temperature,
                    reference_length,
                    reference_velocity,
                    reference_mass_fractions);

    // populate state vector
    auto mass_fractions = new double[n_species*n_states];
    for (int i=0; i<n_species; i++) {
        for (int id=0; id<n_states; id++) {
      mass_fractions[i*n_states+id] = reference_mass_fractions[i];
    }
    }
    auto o_mass_fractions = device.malloc<double>(n_species*n_states, mass_fractions);

    auto temperatures = new double[n_states];
    for (int i=0; i<n_states; i++) temperatures[i] = temperature_K / reference_temperature;
        auto o_temperature_normalized = device.malloc<double>(n_states, temperatures);

    const double pressure = pressure_Pa / reference_pressure;

    // run benchmarks
    if(mode == 1 && n_states) {
      auto o_rates = device.malloc<double>(n_species*n_states);
      auto o_hrr = device.malloc<double>(n_states);

      // warm up
      nekRK::production_rates(n_states,
                              pressure,
                                                            o_temperature_normalized,
                          o_mass_fractions,
                          o_rates,
                          o_hrr,
                      fp32);

      // actual run
      device.finish();
      MPI_Barrier(MPI_COMM_WORLD);
      const auto startTime = MPI_Wtime();
      for(int i=0; i<nRep; i++) {
          nekRK::production_rates(n_states,
                        pressure,
                                            o_temperature_normalized,
                    o_mass_fractions,
                    o_rates,
                    o_hrr,
                    fp32);
      }
      device.finish();
      MPI_Barrier(MPI_COMM_WORLD);
      const auto elapsedTime = MPI_Wtime() - startTime;
      if(rank==0) {
                fprintf(stderr, "avg elapsed time: %.5f s\n", elapsedTime);
                fprintf(stderr, "avg aggregated throughput: %.3f GDOF/s\n",
          (size*(double)(n_states*(n_species+1))*nRep)/elapsedTime/1e9);
      }

      // validate results
      if(ci && rank == 0) {
        // get results from device
        auto rates = new double[n_species*n_states];
        o_rates.copyTo(rates);
        auto hrr = new double[n_states];
        o_hrr.copyTo(hrr);

    const double rtol = (fp32) ? 2e-5 : 1e-14;
    double errInf = fmax(abs((hrr[0] - refBK1Data[0])/(refBK1Data[0])), 0);
        for (int k=0; k < n_species; k++) {
      if(refBK1Data[k+1] > 1e-15)
        errInf = fmax(abs((rates[k*n_states+0] - refBK1Data[k+1])/refBK1Data[k+1]), errInf);
    }
    const int passed = (errInf < rtol);
    fprintf(stderr, "BK1 error_inf: %g (%s)\n", errInf, (passed) ? "passed" : "failed");
        if(!passed) (EXIT_FAILURE);
      }
    }

    if(mode == 2) {
            auto o_thermal_conductivity = device.malloc<double>(n_states);
            auto o_viscosity = device.malloc<double>(n_states);
            auto o_rho_Di = device.malloc<double>(n_species*n_states);

            nekRK::transportCoeffs(n_states,
                                                         pressure_Pa,
                                                    o_temperature_normalized,
                                                    o_mass_fractions,
                                                    o_thermal_conductivity,
                                                    o_viscosity,
                                                    o_rho_Di,
                                                    reference_temperature);

            auto conductivity = new double[n_states];
            o_thermal_conductivity.copyTo(conductivity);
            printf("%15e ", conductivity[0]);
            auto viscosity = new double[n_states];
            o_viscosity.copyTo(viscosity);
            //printf("μ: %1.3e, ", viscosity[0])	;
            printf("%.15e ", viscosity[0])	;
            auto rho_Di = new double[n_species*n_states];
            o_rho_Di.copyTo(rho_Di);
            // print results
            const double K = 1.380649e-23; // J/kelvin
            const double NA = 6.02214076e23; // /mole
            const double R = K*NA;
            const double concentration = reference_pressure / R / reference_temperature;
            const double reference_density = concentration * molar_mass;
            const double reference_time = reference_length / reference_velocity;

            double density = concentration * molar_mass;
            //printf("D: ");
            for (int k=0; k<n_species; k++) {
                double density_times_diffusion_coefficient = rho_Di[k*n_states+0];
                //if(rank==0 && argc > 5) printf("species %5zu density_times_diffusion_coefficient=%.15e\n", k+1, density_times_diffusion_coefficient);
                printf("%.15e ", density_times_diffusion_coefficient/density);
            }
    }

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
