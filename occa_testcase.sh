#!/bin/sh
./run.sh LiDryer
./run.sh grimech30
cat share/mechanisms/LiDryer.c okl/fuego_wrapper.okl > short.okl.txt
cat share/mechanisms/grimech30.c okl/fuego_wrapper.okl > long.okl.txt
~/occa/gentoo/bin/occa translate short.okl.txt
~/occa/gentoo/bin/occa translate long.okl.txt
