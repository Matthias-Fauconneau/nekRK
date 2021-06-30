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
#include <fenv.h>
#include <cstddef>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <unistd.h>
#include "mpi.h"
#include "omp.h"

#include "occa.hpp"
#include <nekrk.h>

int main(int argc, char **argv) {
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    MPI_Init(&argc, &argv);
    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(argc < 5) {
        printf("Usage: ./bk1 SERIAL|CUDA|HIP nStates blockSize nRepetitions [mechanism]\n");
        return 1;
    }
    std::string threadModel; threadModel.assign(strdup(argv[1]));
    const int n_states = std::stoi(argv[2])/size;
    const int blockSize = std::stoi(argv[3]);
    const int nRep = std::stoi(argv[4]);
    std::string mech("LiDryer");
    if(argc > 5) mech.assign(argv[5]);
    const bool verbose = argc < 7;

    char deviceConfig[BUFSIZ];
    const int deviceId = 0;
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
    if(rank == 0 && verbose) {
        std::cerr << "active occa mode: " << device.mode() << '\n';
    }

    nekRK::init(mech.c_str(), device, {}, blockSize, MPI_COMM_WORLD, /*transport:*/false, verbose);
    const int n_species = nekRK::species_names().size();
    const int n_active_species = nekRK::number_of_active_species();
    // setup reference quantities
    double pressure_Pa   = 101325.;
    double temperature_K = 1000.;
    auto mole_fractions = new double[n_species];
    for (int i=0; i<n_species; i++) mole_fractions[i] = 0;
    if (argc>=7) {
        auto values = split(argv[6], " ");
        vector<double> amount_proportions;
        for(auto&& value: values) amount_proportions.push_back(std::stof(value));
        assert(amount_proportions.size() == n_species);
        double sum = 0.;
        for (int i=0; i<n_species; i++) sum += amount_proportions[i];
        for (int i=0; i<n_species; i++) mole_fractions[i] = amount_proportions[i]/sum;
    } else {
        for (int i=0; i<n_species; i++) mole_fractions[i] = 1./(double)n_species;
    }
    auto species_molar_mass = nekRK::species_molar_mass();

    double molar_mass = 0.;
    for(int k=0; k<n_species; k++) molar_mass += mole_fractions[k] * species_molar_mass[k];

    auto reference_mass_fractions = new double[n_species];
    for(int k=0; k<n_species; k++) {
      reference_mass_fractions[k] = mole_fractions[k] * species_molar_mass[k]  / molar_mass;
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

    // populate states
    auto mass_fractions = new double[n_species*n_states];
    for (int i=0; i<n_species; i++) {
        for (int id=0; id<n_states; id++) {
            mass_fractions[i*n_states+id] = reference_mass_fractions[i];
        }
    }
    auto o_mass_fractions = device.malloc<double>(n_species*n_states, mass_fractions);

    auto temperatures = new double[n_states];
    for (int i=0; i<n_states; i++) temperatures[i] = temperature_K / reference_temperature;
    auto o_temperature = device.malloc<double>(n_states, temperatures);

    const double pressure = pressure_Pa / reference_pressure;

    auto o_rates = device.malloc<double>(n_active_species*n_states);
    auto o_heat_release_rate = device.malloc<double>(n_states);

    // warm up
    nekRK::production_rates(n_states,
                            pressure,
                    o_temperature,
                    o_mass_fractions,
                    o_rates,
                    o_heat_release_rate);

    device.finish();
    MPI_Barrier(MPI_COMM_WORLD);
    auto startTime = MPI_Wtime();
    for(int i=0; i<nRep; i++) {
        nekRK::production_rates(n_states,
                    pressure,
                o_temperature,
                o_mass_fractions,
                o_rates,
                o_heat_release_rate);
    }
    device.finish();
    MPI_Barrier(MPI_COMM_WORLD);
    auto elapsedTime = MPI_Wtime() - startTime;
    if(rank==0 && nRep>0) {
        cerr<<"average throughput: "<<((size*(double)(n_states*(n_species+1))*nRep)/elapsedTime/1e9)<<" GDOF/s\n";
    }

    // get results from device
    auto rates = new double[n_active_species*n_states];
    o_rates.copyTo(rates);
    auto heat_release_rate = new double[n_states];
    o_heat_release_rate.copyTo(heat_release_rate);

    double K = 1.380649e-23; // J/kelvin
    double NA = 6.02214076e23; // /mole
    double R = K*NA;
    double concentration = reference_pressure / R / reference_temperature;
    double reference_density = concentration * molar_mass;
    const double reference_time = reference_length / reference_velocity;

    // Prints results
    for (int k=0; k<n_active_species; k++) {
        double mass_production_rate = rates[k*n_states+0];
        // mass_rates[k*n_states + id] = rcp_mass_rate * fg_molar_mass[k] * molar_rates[k];
        double reference_mass_rate = reference_density / reference_time;
        double rcp_mass_rate = 1./reference_mass_rate;
        double molar_rate = mass_production_rate / (rcp_mass_rate * species_molar_mass[k]);
        if(rank==0) printf("%s: %.0f, ", nekRK::species_names()[k].c_str(), molar_rate);
    }
    double molar_heat_capacity_R = nekRK::mean_specific_heat_at_CP_R(reference_temperature, mole_fractions);
    const double energy_rate = (molar_heat_capacity_R * reference_pressure) / reference_time;
    printf("HRR: %.3e\n", heat_release_rate[0] * energy_rate);

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}