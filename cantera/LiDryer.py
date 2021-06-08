#!/usr/bin/env python
from cantera import Solution, gas_constant, one_atm
gas= Solution('LiDryer.yaml')
temperature= 1000
pressure= 101325
print(one_atm)
gas.TPX= temperature, pressure, {"H2": 1, "O2": 1, "H2O": 1, "H": 1, "O": 1, "OH": 1, "HO2": 1, "H2O2": 1, "N2": 1}
time= 1.
density= gas.density
print("%+e" % (gas.cp_mole/1e3))
molar_heat_capacity_R= gas.cp_mole / gas_constant
molar_mass= gas.mean_molecular_weight
mass_rate= density / time;
energy_rate= (molar_heat_capacity_R * pressure) / time;
for i in range(gas.n_total_species):
    if gas.net_production_rates[i] != 0:
        print("%8s %+e" % (gas.species_names[i], gas.net_production_rates[i]))#*gas.molecular_weights[i]/mass_rate))
print("%8s %+e" % ("HRR", gas.heat_release_rate))
