#!/usr/bin/env python
from cantera import Solution
gas= Solution('LiDryer.yaml')
temperature= 1000
pressure= 101325
gas.TPX= temperature, pressure, {"H2": 1, "O2": 1, "H2O": 1, "H": 1, "O": 1, "OH": 1, "HO2": 1, "H2O2": 1, "N2": 1}
for (i, rate) in enumerate(gas.net_rates_of_progress):
    if rate != 0: print("%4d %+e" % (i, rate))
for (name, rate) in zip(gas.species_names, gas.net_production_rates):
    if rate != 0: print("%4s %+e" % (name, rate))
