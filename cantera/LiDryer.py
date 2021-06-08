#!/usr/bin/env python
from cantera import Solution
gas= Solution('LiDryer.yaml')
temperature= 1000
pressure= 101325
gas.TPX= temperature, pressure, {"H2": 1, "O2": 1, "H2O": 1, "H": 1, "O": 1, "OH": 1, "HO2": 1, "H2O2": 1, "N2": 1}
for i in range(gas.n_total_species): if gas.net_production_rates[i] != 0: print("%8s %+e" % (gas.species_names[i], gas.net_production_rates[i]))
