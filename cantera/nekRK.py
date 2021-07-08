#!/bin/env python
from sys import argv
from os import environ
from os.path import dirname, realpath
import subprocess
# TODO: display parameters
# TODO: automatic code regeneration when needed
nekRK = dirname(realpath(argv[0]))+'/..'
run = subprocess.run(cwd=nekRK, args=['build/main', f'mechanisms/{argv[1]}.c'], capture_output=True)
run.check_returncode() #assert len(run.stderr)==0, run.stderr
species, rates = run.stdout.decode().splitlines()
species = species.split()
rates = [float(s) for s in rates.split()]

yaml = f'/usr/share/cantera/data/{argv[1]}.yaml'
cantera_species, cantera = subprocess.run([nekRK+'/cantera/reaction.py', yaml], capture_output=True).stdout.decode().splitlines();
cantera_species = cantera_species.split()
cantera = [float(cantera.split()[cantera_species.index(s)])*1000 for s in species]

fuego_species = cantera_species #s['name'] for s in ruamel.yaml.YAML().load(open(yaml))['species'] # FIXME: Assumes Fuego and Cantera use same specie indices
env = environ
env["LD_LIBRARY_PATH"] = 'build'
run = subprocess.run(cwd='fuego', env=env, args=['build/bk1','Serial','1','1','0',argv[1]], capture_output=True);
assert len(run.stderr)==0, run.stderr
fuego = [float(run.stdout.decode().split()[fuego_species.index(s)]) for s in species]

error = lambda a,b: abs(a-b)/abs(b) if min(abs(a),abs(b)) > 0 else 0 #min(abs(a),abs(b))
print(f"{'':6}: {'NekRK':7} {'Cantera':7} {'Fuego':7} new abs rel old abs rel [kmol/m³/s])")
print('\n'.join([f'{name:6}: {a/1e3:+7.0f} {b/1e3:+7.0f} {c/1e3:+7.0f} {abs(a-b):.0e} {error(a,b):.0e} {abs(c-b):.0e} {error(c,b):.0e}' for (name, a, b, c) in zip(species, rates, cantera, stgeke) if a != 0 and b != 0]))
