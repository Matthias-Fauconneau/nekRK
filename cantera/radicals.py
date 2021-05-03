#!/usr/bin/env python
print("with all radicals (except H) to see how large part of the mechanism becomes active")
import cantera as ct
import numpy   as np
gas=ct.Solution('gri30.cti')
Tref = 1000     # [K]
pref = 1.0      # [atm]
gas.TPX = Tref, pref*ct.one_atm, {'O2':0.4, 'CH4':0.2, 'O': 1e-5, 'OH': 1e-5, 'HO2': 1e-5, 'H2O2': 1e-5, 'CH3': 1e-5, 'CH': 1e-5, 'N2':0.4-6*1e-5}
tRef   = 1.0                    # [s]
rhoRef = gas.density*1e3/1e6    # [gr/cm3]
cpRef  = gas.cp_mass*1e7/1e3    # [erg/(gr K)]
factorW = tRef/rhoRef
factorT = Tref*cpRef
for i in range(gas.n_total_species):
    if (gas.net_production_rates[i] !=0):
        print("%6s %+.7e" % ( gas.species_names[i], gas.net_production_rates[i]*1e3/1e6*gas.molecular_weights[i]*factorW ) )
