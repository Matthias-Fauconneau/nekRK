#!/usr/bin/env python

import cantera as ct
import numpy   as np

#PelePhysics/Support/Fuego/Mechanism/Models/LiDryer
#ck2yaml --input=mechanism.inp --thermo=therm.dat --transport=tran.dat --permissive
#mv mechanism.yaml /usr/share/cantera/data/LiDryer.yaml
gas=ct.Solution('LiDryer.yaml', thermo='IdealGas')
gas.transport_model = 'Mix'
gas.TPX = 1000, ct.one_atm, {'H2': 1./9., 'O2': 1./9., 'H2O': 1./9., 'H': 1./9., 'O': 1./9., 'OH': 1./9., 'HO2': 1./9., 'H2O2': 1./9., 'N2': 1./9.}
gas()
print(ct.Species.listFromFile('gri30.yaml'))
print(gas.mix_diff_coeffs)
