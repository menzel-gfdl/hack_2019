# Princeton GPU Hackathon 2019
This repository describes how to build/run kernels from FV3 on the OLCF Ascent system.

# Requirements
The GNU autotools build system and a fortran compiler.

# Configuring the environment
Since this project is targeting GPUs via OpenACC, we recommend configuring your environment
to use the PGI compiler.
```
$ source $MODULESHOME/init/bash
$ module swap xl pgi/19.5
```

# Building
This repository has been set up to use the GNU autotools build system.  Running
```
$ autoreconf --install
$ ./configure CC=mpicc FC=mpifort FCFLAGS="-g -acc -ta=nvidia:cc70 -Minfo=accel -Mcuda=lineinf"
$ make
```
should produce an OpenACC-enabled executable gfdl_fyppm/src/fyppm.

# Running
A bash script is provided that will submit the job to the scheduler.  It can be run:
```
$ ./submit.sh gfdl_fyppm/src/fyppm <walltime in minutes> <# of MPI ranks> <# of gpus>
```

# Checking For Correctness
TBD

# Debugging
Start ddt in the background on the login node, then use the submission script with the -d
option:
```
$ module load forge/19.0.2
$ ddt &
$ ./submit.sh -d gfdl_fyppm/src/fyppm <walltime in minutes> <# of MPI ranks> <# of gpus>
```
Once the job starts running, a dialog box should appear in ddt asking if you'd like to
accept a "reverse connect" from you job.  Click accept, then proceed to use ddt as
normal.

# Profiling
The provided bash submission script also includes a profiling option (-p or --profile), which
will allow the run the be profiled by nvprof:
```
$ ./submit.sh -p <name of profile output file> gfdl_fyppm/src/fyppm <walltime in minutes> <# of MPI ranks> <# of gpus>
```
Once the run completes, the output file containing the profile can be run through the
visual profiler (nvvp):
```
$ module load cuda/9.2.148
$ nvvp <path to profile output file>
```
