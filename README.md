# nekRK: nek reaction kinetics

nekRK is a software package for nekRS to compute species production rates and transport properties for chemically reactive flows.

## How to compile
Dependencies:
- CMake-3.11
- MPI
- [OCCA](https://github.com/libocca/occa) 
- Python

```sh
export NEKRK_PATH=$HOME/.local/nekRK
fuego/run.sh GRIMech-3.0.ck2 (only required if mechanisms does not exist in share/mechanism)
cmake -B build -DOCCA_DIR=$OCCA_DIR -DCMAKE_INSTALL_PREFIX=$NEKRK_PATH 
cd build && make && make install
```
Please ensure that the env-var `$OCCA_DIR` points to your OCCA installation. 

## Benchmark Kernels

### BK1: Species production rates

```sh
Usage: ./bk1 SERIAL|CUDA|HIP nStates blockSize nRepetitions [mechanism]
>cd $NEKRK_PATH; bin/bk1 CUDA 100000 256 1000
mechanism file: ./share/mechanisms/GRIMech-3.0.c
nSpecies: 53
avg throughput: 2.609 GDOF/s
```

| CPU/GPU           | MECH    | GDOF/s |
| ----------------- | ------- | ------ |
| Nvidia V100       | GRI 3.0 |  2.61  | 
| Nvidia A100       | GRI 3.0 |  3.98  |
| AMD MI100         | GRI 3.0 |  ?.??  |
| 2xAMD EPYC 7742   | GRI 3.0 |  0.63  |
| 2xIntel XEON 6252 | GRI 3.0 |  0.17  |
|                   |         |        | 
| Nvidia V100       | LiH2    | 21.57  |
| Nvidia A100       | LiH2    | 28.67  | 
| AMD MI100         | LiH2    | ??.??  |
| 2xAMD EPYC 7742   | LiH2    |  1.68  |
| 2xIntel XEON 6252 | LiH2    |  0.46  |

### BK2: Properties

TODO

## Benchmark Problems 

### BP1: Ignition 0D-reactor

TODO
