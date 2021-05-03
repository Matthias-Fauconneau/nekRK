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

## Results BK1

### NVidia A100
```sh
>cd $NEKRK_PATH
>bin/bk1 CUDA 100000 256 1000
active occa mode: CUDA
mechanism file: GRIMech-3.0.c
nSpecies: 53
avg throughput: 73.0634 MStates/s
```

### 2 x AMD EPYC 7742 
```sh
>cd $NEKRK_PATH
>OCCA_CXXFLAGS="-O3 -ffast-math" mpirun -np 48 bin/bk1 SERIAL 100000 1 100 
active occa mode: SERIAL
mechanism file: GRIMech-3.0.c
nSpecies: 53
avg throughput: 11.6161 MStates/s
```
