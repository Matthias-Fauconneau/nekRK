#!/bin/env python
from sys import argv
import subprocess
# TODO: automatic code regeneration when needed
species, rates = subprocess.run(["build/main", argv[2]], capture_output=True).stdout.decode().splitlines()
species = species.split()
rates = [float(s)/1000 for s in rates.split()]

cantera_species, cantera = subprocess.run(["cantera/reaction.py", argv[1]], capture_output=True).stdout.decode().splitlines();
cantera_species = cantera_species.split()
cantera = [float(cantera.split()[cantera_species.index(s)]) for s in species]

error = lambda a,b: abs(a-b)/min(abs(a),abs(b)) if min(abs(a),abs(b)) > 0 else 0
print(f"{'':6}: {'NekRK':7} {'Cantera':7} Relative Error [mol/m³/s])")
print('\n'.join([f'{name:6}: {a:+7.0f} {b:+7.0f} {error(a,b):.0e}' for (name, a, b) in zip(species, rates, cantera) if a != 0 and b != 0]))
