#!/bin/fish
#echo Cantera
#cantera/transport.py
echo Rust
run ~/combustion transport,itertools,float-pretty-print,program > share/mechanisms/$argv[1].viscosity.program || exit 1
#run ~/combustion transport,itertools,float-pretty-print,reaction || exit 1
./main.py "/usr/share/cantera/data/$argv[1].yaml" > /var/tmp/$argv[1].c || exit 1
mv /var/tmp/$argv[1].c share/mechanisms
fish -c 'make -Cbuild -j --no-print-directory | rg -v Built | rg -v Consolidate | rg -v Linking | rg -v Building; exit $pipestatus[1]' || exit 1
#build/bk1 Serial 1 1 0 $argv[1] || exit 1
echo NekRK
build/bk2 Serial 1 1 0 $argv[1] quiet
