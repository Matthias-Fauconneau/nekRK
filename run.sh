#!/bin/env fish
cd (dirname (status --current-filename))
set arg $argv[1]
set file (test -e $arg && echo $arg || echo "/usr/share/cantera/data/$arg.yaml")
echo -e 'Cantera\t| '
cantera/reaction.py $file
cantera/transport.py $file
#echo -e 'Rust\t| '
#run combustion $file || exit 1
set name (basename $file)
echo -e "NekRK $name\t| "
set target share/mechanisms/$name.c
command test $target -nt $file -a $target -nt main.py || begin; echo 'Fit binary thermal diffusion coefficients using Python (slow version)' && ./main.py $file > /var/tmp/$name.c && mv /var/tmp/$name.c share/mechanisms || exit 1; end
#cmake -Bdebug -DOCCA_DIR=../occa/debug -DCMAKE_BUILD_TYPE=Debug
echo "Build tests"
fish -c 'make -Cdebug -j --no-print-directory | rg -v Built | rg -v Consolidate | rg -v Linking | rg -v Building; exit $pipestatus[1]' || exit 1
test -e share/mechanisms/$name.c || exit 1
echo "Rates"
PATH=/usr/lib/llvm/12/bin:/usr/bin debug/bk1 Serial 1 1 0 $name || exit 1
echo "Transport"
PATH=/usr/lib/llvm/12/bin:/usr/bin debug/bk2 Serial 1 1 0 $name
