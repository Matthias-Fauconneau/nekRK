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
#fuego/run.sh GRIMech-3.0.ck2
cmake -B build -DOCCA_DIR=$OCCA_DIR -DCMAKE_INSTALL_PREFIX=$NEKRK_PATH 
cd build && make && make install
```

## Running BP1 on NVIDIA V100
```sh
>cd $NEKRK_PATH
>bin/bk1 CUDA 100000 256 1000
active occa mode: CUDA
throughput: 48.0111 MStates/s
```
