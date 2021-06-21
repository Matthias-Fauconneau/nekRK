#!/bin/env python3.9
#
import sys, cantera
gas= cantera.Solution(sys.argv[1])
temperature= 1000
pressure= 101325
gas.TPX= (temperature, pressure, [1 for specie in gas.species_names])
print(', '.join([f'{name}: {rate:+.3f}' for (name, rate) in zip(gas.species_names, gas.net_production_rates) if rate != 0]+[f'HRR: {gas.heat_release_rate:.3e}']))
