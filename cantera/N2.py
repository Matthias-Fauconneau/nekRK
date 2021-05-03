#!/usr/bin/env python
print("with the more realistic N2, reactions with N2 as either reactant or product are activated")
from cantera import Solution, gas_constant, one_atm
gas= Solution('gri30.cti')
temperature= 1000
pressure= one_atm
gas.TPX= temperature, pressure, {'O2':0.4, 'CH4':0.2, 'N2':0.4}
time= 1.
density= gas.density
print("%+e" % (gas.cp_mole/1e3))
molar_heat_capacity_R= gas.cp_mole / gas_constant
molar_mass= gas.mean_molecular_weight
mass_rate= density / time;
energy_rate= (molar_heat_capacity_R * pressure) / time;
for i in range(gas.n_total_species):
    if gas.net_production_rates[i] != 0:
        print("%8s %+e" % (gas.species_names[i], gas.net_production_rates[i]*gas.molecular_weights[i]/mass_rate))
print("%8s %+e" % ("HRR", gas.heat_release_rate))
print("%8s %+e" % ("HRR0", energy_rate))
print("%8s %+e" % ("HRR/HRR0", gas.heat_release_rate/energy_rate))
