#!/bin/bash
fuego=$(pwd)/$(dirname "$0")
name=$(basename $1)
mkdir -p share/mechanisms
output=$(pwd)/share/mechanisms/$name.c
echo $1 "->" $output
#/!\ fmc/fuego accepts a --trans parameter which is dead code, transport properties can only be loaded from a modified reaction mechanism
# There MUST be a therm.dat file in the same folder (It is not optional) (/!\ The flag is thermo but the file is therm)
cd $1
PYTHONPATH=$fuego:$fuego/fuego:$fuego/pyre:$fuego/journal:$fuego/weaver python2 $fuego/applications/fmc.py -mechanism=mechanism.inp -thermo=therm.dat -name=$output &&\
sed -i 's/+ -/-/g' $output
cat $fuego/verbatim.c >> $output
