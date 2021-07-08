#!/bin/env python
from sys import argv
from os import getenv
import subprocess
# TODO: display parameters
# TODO: automatic code regeneration when needed
species, rates = subprocess.run(['build/main', argv[2]], capture_output=True).stdout.decode().splitlines()
species = species.split()
rates = [float(s) for s in rates.split()]

cantera_species, cantera = subprocess.run(['cantera/reaction.py', argv[1]], capture_output=True).stdout.decode().splitlines();
cantera_species = cantera_species.split()
cantera = [float(cantera.split()[cantera_species.index(s)])*1000 for s in species]

#ln -s ../okl # only supports GRIMech-3.0
stgeke_species, stgeke = 'H2      H       O       O2      OH      H2O     HO2     H2O2 C       CH      CH2     CH2(S)  CH3     CH4     CO      CO2 HCO     CH2O    CH2OH   CH3O    CH3OH   C2H     C2H2    C2H3 C2H4    C2H5    C2H6    HCCO    CH2CO   HCCOH   N       NH NH2     NH3     NNH     NO      NO2     N2O     HNO     CN HCN     H2CN    HCNN    HCNO    HOCN    HNCO    NCO     N2 AR      C3H7    C3H8    CH2CHO  CH3CHO', subprocess.run(cwd=getenv('HOME')+'/stgeke/build',args=['bk1','Serial','1','1','0',argv[3]], capture_output=True).stdout.decode() # FIXME: parse Chemkin model
stgeke_species = stgeke_species.split()
stgeke = [float(stgeke.split()[stgeke_species.index(s)]) for s in species]

error = lambda a,b: abs(a-b)/abs(b) if min(abs(a),abs(b)) > 0 else 0 #min(abs(a),abs(b))
print(f"{'':6}: {'NekRK':7} {'Cantera':7} {'Stgeke':7} new abs rel old abs rel [kmol/m³/s])")
print('\n'.join([f'{name:6}: {a/1e3:+7.0f} {b/1e3:+7.0f} {c/1e3:+7.0f} {abs(a-b):.0e} {error(a,b):.0e} {abs(c-b):.0e} {error(c,b):.0e}' for (name, a, b, c) in zip(species, rates, cantera, stgeke) if a != 0 and b != 0]))
