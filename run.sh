#!/bin/env fish
cd (dirname (status --current-filename))
set arg $argv[1]
set file (test -e $arg && echo $arg || echo "/usr/share/cantera/data/$arg.yaml")
echo -e 'Cantera\t| '
cantera/reaction.py $file
cantera/transport.py $file
#echo -e 'Rust\t| '
#run combustion $file || exit 1
echo -e 'NekRK\t| '
set name (basename $file)
./main.py $file > /var/tmp/$name.c || exit 1
mv /var/tmp/$name.c share/mechanisms
fish -c 'make -Cbuild -j --no-print-directory | rg -v Built | rg -v Consolidate | rg -v Linking | rg -v Building; exit $pipestatus[1]' || exit 1
build/bk1 Serial 1 1 0 $name quiet || exit 1
build/bk2 Serial 1 1 0 $name quiet
