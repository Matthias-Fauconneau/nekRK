#!/bin/fish
cd (dirname $argv[1])
echo -ne 'Cantera\t| '
cantera/transport.py
echo -ne 'Rust\t| '
#run combustion transport,itertools || exit 1
run combustion transport,debug > share/mechanisms/$argv[1].diffusion.program || exit 1
echo -ne 'NekRK\t| '
./main.py "/usr/share/cantera/data/$argv[1].yaml" > /var/tmp/$argv[1].c || exit 1
mv /var/tmp/$argv[1].c share/mechanisms
fish -c 'make -Cbuild -j --no-print-directory | rg -v Built | rg -v Consolidate | rg -v Linking | rg -v Building; exit $pipestatus[1]' || exit 1
#build/bk1 Serial 1 1 0 $argv[1] || exit 1
build/bk2 Serial 1 1 0 $argv[1] quiet
