#!/usr/bin/env python

import cantera as ct
import numpy   as np

gas=ct.Solution('gri30.yaml', thermo='IdealGas')
gas.transport_model = 'Mix'
gas.TPX = 1000, ct.one_atm, {'O2':0.40, 'CH4':0.20, 'AR':0.40}
#gas.TPX = 1000, ct.one_atm, {'O2':0.40, 'CH4':0.20, 'N2':0.40}
gas()

print(ct.Species.listFromFile('gri30.yaml'))
print(gas.mix_diff_coeffs)
