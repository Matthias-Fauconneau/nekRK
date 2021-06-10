#!/usr/bin/env python
from cantera import Solution
#PelePhysics/Support/Fuego/Mechanism/Models/LiDryer
#ck2yaml --input=mechanism.inp --output=/usr/share/cantera/data/LiDryer.yaml #--thermo=therm.dat --transport=tran.dat
gas= Solution('LiDryer.yaml', thermo='IdealGas')
temperature= 1000
pressure= 101325
gas.TPX = temperature, pressure, {'H2': 1, 'O2': 1, 'H2O': 1, 'H': 1, 'O': 1, 'OH': 1, 'HO2': 1, 'H2O2': 1, 'N2': 1}
gas()
print(gas.mix_diff_coeffs)
