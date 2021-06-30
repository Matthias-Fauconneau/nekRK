#!/bin/env fish
cd (dirname (status --current-filename))
set arg $argv[1]
set file (test -e $arg && echo $arg || echo "/usr/share/cantera/data/$arg.yaml")
echo -e 'Cantera\t| '
cantera/reaction.py $file
cantera/transport.py $file
echo -e 'Rust\t| '
run combustion transport $file || exit 1
set name (basename $file .yaml)
echo -e "NekRK $name\t| "
set target mechanisms/$name.c
command test $target -nt $file -a $target -nt main.py || begin; echo 'Fit binary thermal diffusion coefficients using Python (slow version)' && ./main.py $file > /var/tmp/$name.c && mv /var/tmp/$name.c mechanisms || exit 1; end
#set build build; cmake -B$build -DOCCA_DIR=../occa/build -DCMAKE_BUILD_TYPE=Release
set build debug; cmake -B$build -DOCCA_DIR=../occa/build -DCMAKE_BUILD_TYPE=Debug
fish -c "make -C$build -j --no-print-directory | rg -v Built | rg -v Consolidate | rg -v Linking | rg -v Building; exit \$pipestatus[1]" || exit 1
test -e mechanisms/$name.c || exit 1
PATH=/usr/lib/llvm/12/bin:/usr/bin $build/main Serial 1 1 0 $name || exit 1
