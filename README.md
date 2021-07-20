# nekRK: nek reaction kinetics

nekRK is a software package for [nekRS](https://github.com/Nek5000/nekRS) to compute species production rates and transport properties for chemically reactive flows. The repository contains a modified version of the source code generation tool FUEGO (originally created by Michael Aivazis at CalTech) in order to support the inclusion of chemical models specified in the CHEMKIN-II format.

## How to compile
Dependencies:
- CMake-3.11
- MPI
- [OCCA](https://github.com/libocca/occa)
- Python

```sh
export NEKRK_PATH=$HOME/.local/nekRK
./main.py /usr/share/cantera/data/gri30.yaml > share/mechanisms/gri30.c (only required if not commited)
cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$NEKRK_PATH ..
module load spectrum_mpi
make -Cbuild -j install
```

## Benchmark Kernels

### BK1: Source Terms / Production Rates

```sh
Usage: module load spectrum_mpi cuda && rm ~/.occa -R && mpirun -np 1 $NEKRK_PATH/bin/bk --mode 1|2 --backend SERIAL|CUDA|HIP --n-states n [--block-size n] [--repetitions n] [--fp32] [--ci] [--mechanism-name]
>cd $NEKRK_PATH; mpirun -np 1 ./bin/bk --mode 1 --backend CUDA --n-states 1000000
number of states: 1000000
use fp32: 0
number of repetitions: 1000
nekRK initialized successfully
mechanism file: GRIMech-3.0.c
nSpecies: 53
active occa mode: CUDA
blockSize: 512
avg elapsed time: 8.04059 s
avg aggregated throughput: 6.716 GDOF/s
```

| CPU/GPU           | MECH        | GDOF/s | GRXN/s |
| ----------------- | ----------- | ------ | ------ |
| Nvidia V100       | GRIMech-3.0 |  3.59  |  21.4  |
| Nvidia A100       | GRIMech-3.0 |  6.71  |  40.1  |
| AMD MI100         | GRIMech-3.0 |  tbd   |  tbd   |
| 2xAMD EPYC 7742   | GRIMech-3.0 |  0.63  |  3.79  |
| 2xIntel XEON 6252 | GRIMech-3.0 |  0.17  |  1.02  |
|                   |             |        |        |
| Nvidia V100       | LiH2        |  tbd   |  tbd   |
| Nvidia A100       | LiH2        | 31.9   |  66.9  |
| AMD MI100         | LiH2        |  tbd   |  tbd   |
| 2xAMD EPYC 7742   | LiH2        |  tbd   |  tbd   |
| 2xIntel XEON 6252 | LiH2        |  tbd   |  tbd   |

### BK2: Mixture-Averaged Transport Properties

TODO

## Benchmark Problems

### BP1: Ignition 0D-reactor

TODO

## Related Work
 * https://github.com/SLACKHA/pyJac
 * https://github.com/AMReX-Combustion/PelePhysics
 * https://dl.acm.org/doi/10.1145/2692916.2555258
 * https://opencommons.uconn.edu/dissertations/2044
 * https://link.springer.com/chapter/10.1007%2F978-3-319-68394-2_11
