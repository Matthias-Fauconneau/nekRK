#!/bin/fish
#cmake -B build OCCA_DIR=~/occa/gentoo
run ~/combustion transport,itertools,float-pretty-print || exit 1
#run ~/combustion reaction 1>LiDryer.c || exit 1
~/.local/share/junest/bin/junest -b "--bind $HOME/nekRK /mnt" fish -c "cd /mnt && fuego/run.sh etc/fuego/mechanisms/$argv[1]" || exit 1
fish -c 'make -Cbuild -j --no-print-directory | rg -v Built; exit $pipestatus[1]' && build/bk2 Serial 1 1 0 $argv[1] quiet
