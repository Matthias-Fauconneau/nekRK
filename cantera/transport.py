#!/bin/env python3.9
#PelePhysics/Support/Fuego/Mechanism/Models/LiDryer
#ck2yaml --input=mechanism.inp --output=mechanisms/LiDryer.yaml #--thermo=therm.dat --transport=tran.dat
from cantera import Solution
gas= Solution('LiDryer.yaml', thermo='IdealGas')
temperature= 1000
pressure= 101325
gas.TPX = temperature, pressure, {'H2': 1, 'O2': 1, 'H2O': 1, 'H': 1, 'O': 1, 'OH': 1, 'HO2': 1, 'H2O2': 1, 'N2': 1}
print(f"μ: {gas.viscosity:.4}, λ: {gas.thermal_conductivity:.4}, D: {' '.join([f'{x:.3e}' for x in gas.mix_diff_coeffs])}")
