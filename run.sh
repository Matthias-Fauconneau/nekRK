#!/bin/fish
#cmake -B build OCCA_DIR=~/occa/gentoo
cargo run --manifest-path ~/combustion/Cargo.toml 1>LiDryer.c || exit 1
~/.local/share/junest/bin/junest -b "--bind $HOME/nekRK /mnt" fish -c 'cd /mnt && fuego/run.sh etc/fuego/mechanisms/LiDryer'
fish -c 'make -Cbuild -j --no-print-directory | rg -v Built; exit $pipestatus[1]' && build/bk1 Serial 1 1 0 LiDryer quiet
