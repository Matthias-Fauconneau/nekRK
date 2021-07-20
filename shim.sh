#!/bin/env fish
./main.py $argv
string replace float dfloat
cat shim.okl
#ifdef __FG_ENABLE__ #endif
__FG_CONST__ #__constant__
__FG_DEVICE__ #__device__
