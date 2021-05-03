# nekRK: nek reaction kinetics

nekRK is a software package for nekRS to compute species production rates and transport properties for chemically reactive flows.

## How to compile
Dependencies:
- CMake-3.11
- MPI
- OCCA
- Python

```sh
export NEKRK_PATH=$HOME/.local/nekRK
fuego/run.sh GRIMech-3.0.ck2 (only required if mechanisms does not exist in share/mechanism)
cmake -B build -DOCCA_DIR=$OCCA_DIR -DCMAKE_INSTALL_PREFIX=$NEKRK_PATH 
cd build && make && make install
```
Please ensure that env-var `$OCCA_DIR` points to your OCCA installation. 

## Benchmark Kernels/Problems
* BK1: Species production rates
* BK2: Properties (not available yet)

* BP1: Ignition of 0D-reactor (not available yet)

## Running BK1 on NVIDIA V100
```sh
>cd $NEKRK_PATH
>bin/bk1 CUDA 100000 256 100
active occa mode: CUDA
mechanism: GRIMech-3.0
nSpecies: 53
throughput: 48.0111 MStates/s
```

## Running BK1 on 2xIntel 6252 
```sh
>cd $NEKRK_PATH
>OCCA_CXXFLAGS="-O3 -ffast-math" mpirun -np 48 bin/bk1 SERIAL 100000 256 20
active occa mode: CUDA
mechanism: GRIMech-3.0
nSpecies: 53
throughput: 3.11479 MStates/s
```
