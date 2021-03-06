# nekRK: reaction kinetics

nekRK is a software package for [nekRS](https://github.com/Nek5000/nekRS) to compute species production rates and transport properties for chemically reactive flows. The repository contains a modified version of the source code generation tool FUEGO (originally created by Michael Aivazis at CalTech) in order to support the inclusion of chemical models specified in the CHEMKIN-II format.

## How to compile
Dependencies:
- C++
- Python3
- CMake-3.11
- [OCCA](https://github.com/libocca/occa) (/!\ The occa source directory needs to stay available at the same path as it is used by the OCCA installation)
- MPI

```sh
git clone --recurse-submodules https://github.com/Matthias-Fauconneau/nekRK
export OCCA_DIR=$HOME/occa
export NEKRK_PATH=$HOME/.local/nekRK
./main.py mechanisms/grimech30.yaml > share/mechanisms/grimech30.c (only required if mechanisms does not exist in share/mechanism)
cmake -B build -DOCCA_DIR=$OCCA_DIR -DCMAKE_INSTALL_PREFIX=$NEKRK_PATH
cd build && make -j && make install
```
Please ensure that the env-var `$OCCA_DIR` points to your OCCA installation.

## Benchmark Kernels

### BK1: Species production rates

```sh
Usage: ./bk1 mechanism SERIAL|CUDA|HIP nStates blockSize nRepetitions
>cd $NEKRK_PATH; bin/main LiDryer CUDA 100000 256 1000
mechanism file: ./share/mechanisms/GRIMech-3.0.c
nSpecies: 53
avg throughput: 2.609 GDOF/s
```

| CPU/GPU           | MECH        | GDOF/s | GRXN/s |
| ----------------- | ----------- | ------ | ------ |
| Nvidia V100       | GRIMech-3.0 |  2.61  |  15.7  |
| Nvidia A100       | GRIMech-3.0 |  3.98  |  23.9  |
| AMD MI100         | GRIMech-3.0 |  2.00  |  12.0  |
| 2xAMD EPYC 7742   | GRIMech-3.0 |  0.63  |  3.79  |
| 2xIntel XEON 6252 | GRIMech-3.0 |  0.17  |  1.02  |
|                   |             |        |        |
| Nvidia V100       | LiH2        |  21.6  |  45.4  |
| Nvidia A100       | LiH2        |  28.7  |  60.3  |
| AMD MI100         | LiH2        |  14.3  |  30.3  |
| 2xAMD EPYC 7742   | LiH2        |  1.68  |  3.53  |
| 2xIntel XEON 6252 | LiH2        |  0.46  |  0.97  |

### BK2: Properties

```sh
Usage: ./bk1 SERIAL|CUDA|HIP nStates blockSize nRepetitions [mechanism]
>cd $NEKRK_PATH; bin/bk2 CUDA 100000 256 1000
```

## Benchmark Problems

### BP1: Ignition 0D-reactor

TODO

