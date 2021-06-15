#!/bin/fish
run ~/combustion transport,itertools,float-pretty-print,reaction || exit 1
./main.py "/usr/share/cantera/data/$argv[1].yaml" > /var/tmp/$argv[1].c || exit 1
mv /var/tmp/$argv[1].c share/mechanisms
fish -c 'make -Cbuild -j --no-print-directory | rg -v Built; exit $pipestatus[1]' || exit 1
build/bk1 Serial 1 1 0 $argv[1] || exit 1
#build/bk2 Serial 1 1 0 $argv[1]
