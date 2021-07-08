#!/bin/env python3.9
#ck2yaml --input=mechanism.inp --output=mechanisms/$name.yaml #--thermo=therm.dat --transport=tran.dat
import sys, cantera
gas= cantera.Solution(sys.argv[1])
temperature= 1000
pressure= 101325
gas.TPX= (temperature, pressure, [1 for specie in gas.species_names])
print(f"{gas.viscosity} {gas.thermal_conductivity} {' '.join([f'{x:e}' for x in gas.mix_diff_coeffs])}")
#print(f"μ: {gas.viscosity:.4}, λ: {gas.thermal_conductivity:.4}, D: {' '.join([f'{x:.3e}' for x in gas.mix_diff_coeffs])}")
