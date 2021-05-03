Dependencies:
- CMake-3.11
- MPI
- OCCA
- Python

## How to compile
```sh
export NEKRK_PATH=$HOME/.local/nekRK
#fuego/run.sh GRIMech-3.0.ck2
cmake -B build -DOCCA_DIR=$OCCA_DIR -DCMAKE_INSTALL_PREFIX=$NEKRK_PATH 
cd build && make && make install
```

# Running BP1 on NVIDIA V100
```sh
>cd $NEKRK_PATH
>bin/bk1 CUDA 100000 256 1000
active occa mode: CUDA
throughput: 48.0111 MStates/s
```
