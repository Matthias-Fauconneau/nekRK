#!/bin/fish
run ~/combustion transport,itertools,float-pretty-print || exit 1
./main.py '/usr/share/cantera/data/LiDryer.yaml'
#fish -c 'make -Cbuild -j --no-print-directory | rg -v Built; exit $pipestatus[1]' && build/bk2 Serial 1 1 0 $argv[1] quiet
