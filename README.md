# Princeton GPU Hackathon 2019
This repository describes how to build/run kernels from FV3 on the OLCF Ascent system.

# Requirements
The GNU autotools build system and a fortran compiler.

# Configuring the environment
```
$ module swap xl pgi/19.5
```

# Building
```
$ autoreconf --install
$ ./configure CC=mpicc FC=mpifort FCFLAGS="-g -acc -ta=nvidia:cc70 -Minfo=accel -Mcuda=lineinf"
$ make
```

# Running
A bash script is provided that will submit the job to the scheduler.  It can be run:
```
./submit.sh <path to executable> <walltime in minutes> <# of MPI ranks> <# of gpus>
```

# Checking For Correctness
TBD
