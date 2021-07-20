#!/bin/env python
# TODO: display parameters
# TODO: automatic code regeneration when needed
from sys import argv
from os import environ
from os.path import dirname, realpath
import subprocess

print('NekRK')
nekRK = dirname(realpath(argv[0]))

build = subprocess.run(cwd=nekRK, args=['make',"-Cbuild",'install'])
build.check_returncode()

#$NEKRK_PATH/bin
args = f'build/bk --mechanism-name {argv[1]} --backend Serial --n-states 1 --repetitions 0 --mode 1 --ci'
print(args)
run = subprocess.run(cwd=nekRK, args=args.split(), capture_output=True)
run.check_returncode()
print(run.stderr.decode())
print(run.stdout.decode())
species, rates = run.stdout.decode().splitlines()
species = species.split()
rates = [float(s) for s in rates.split()]

args = f'build/bk --mechanism-name {argv[1]} --backend Serial --n-states 1 --repetitions 0 --mode 2'
print(args)
run = subprocess.run(cwd=nekRK, args=args.split(), capture_output=True)
run.check_returncode() #assert len(run.stderr)==0, run.stderr
print(run.stderr.decode())
print(run.stdout.decode())
_species, nekRK_transport = run.stdout.decode().splitlines()
nekRK_transport = [float(s) for s in nekRK_transport.split()]

print('Cantera')
yaml = f'/usr/share/cantera/data/{argv[1]}.yaml'
cantera_species, cantera, cantera_transport = subprocess.run([nekRK+'/gas.py', yaml], capture_output=True).stdout.decode().splitlines();
cantera_species = cantera_species.split()
cantera = cantera.split()
cantera = [float(cantera[cantera_species.index(s)])*1000 for s in species]
conductivity, viscosity, *diffusion = [float(s) for s in cantera_transport.split()]
diffusion = [diffusion[cantera_species.index(s)] for s in species]
cantera_transport = conductivity, viscosity, *diffusion

error = lambda a,b: abs(a-b)/abs(b) if min(abs(a),abs(b)) > 0 else 0 #min(abs(a),abs(b))

print(f"{'':6}: {'NekRK':8} {'Cantera':8} {'abs':6} rel")
print('\n'.join([f'{name:6}: {a/1e6:+8.2f} {b/1e6:+8.2f} {abs(a-b)/1e6:.0e} {error(a,b):.0e} ' for (name, a, b) in zip(species, rates, cantera) if a != 0 and b != 0]))
print('\n'.join([f'{name:6}: {a:.2e} {b:.2e} {abs(a-b):.0e} {error(a,b):.0e}' for (name, a, b) in zip(['λ','μ']+species, nekRK_transport, cantera_transport) if a != 0 and b != 0]))
