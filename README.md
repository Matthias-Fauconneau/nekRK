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

### BK1: Species production rates

| CPU/GPU         | MECH    | GDOF/s |
| --------------- | ------- | ------ |
| Nvidia V100     | GRI 3.0 |  2.61  | 
| Nvidia A100     | GRI 3.0 |  3.98  |
| AMD MI100       | GRI 3.0 |  ???   |
| 2xAMD EPYC 7742 | GRI 3.0 |  0.63  |
| 2xIntel 6252    | GRI 3.0 |  0.17  |
|                 |         |        | 
| Nvidia V100     | LiH2    | 21.57  |
| Nvidia A100     | LiH2    | 28.67  | 
| AMD MI100       | LiH2    |  ???   |
| 2xAMD EPYC 7742 | LiH2    |  1.68  |
| 2xIntel 6252    | LiH2    |  0.46  |

### BK2: Properties

TODO

### BP1: Ignition of 0D-reactor

TODO
