#!/bin/bash
mkdir -p share/mechanisms
output=`pwd`/share/mechanisms/`echo $1 | sed 's/[^.]*$/c/'`
cd fuego
./fmc.py --mechanism=$1 --name=$output &&\
sed -i 's/+ -/-/g' $output
