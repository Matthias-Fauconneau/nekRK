#!/bin/bash
die() { echo "$*" 1>&2 ; exit 1; }
cd "$(dirname "$0")"
./run.py `echo $1 | sed 's/[^.]*$//'` || die
mkdir -p share/mechanisms
output=`pwd`/share/mechanisms/`echo $1 | sed 's/[^.]*$/c/'`
PYTHONPATH=.:fuego:pyre:journal:weaver python2 applications/fmc.py --mechanism=$1 --trans=$1.transport.dat --name=$output &&\
sed -i 's/+ -/-/g' $output
