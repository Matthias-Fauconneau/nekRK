#!/bin/env python3.9
#ck2yaml --input=mechanism.inp --output=mechanisms/$name.yaml #--thermo=therm.dat --transport=tran.dat
import sys, cantera
gas= cantera.Solution(sys.argv[1])
temperature= 1000
pressure= 101325
gas.TPX= (temperature, pressure, [1 for specie in gas.species_names])
print(' '.join(gas.species_names))
print(' '.join(f'{r}' for r in gas.net_production_rates))
print(f"{gas.viscosity} {gas.thermal_conductivity} {' '.join([f'{x:e}' for x in gas.mix_diff_coeffs])}")
#print(', '.join([f'{name}: {rate:+.3f}' for (name, rate) in zip(gas.species_names, gas.net_production_rates) if rate != 0]+[f'HRR: {gas.heat_release_rate:.3e}']))
#print(f"μ: {gas.viscosity:.4}, λ: {gas.thermal_conductivity:.4}, D: {' '.join([f'{x:.3e}' for x in gas.mix_diff_coeffs])}")
