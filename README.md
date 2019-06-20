# Princeton GPU Hackathon 2019
This repository describes how to build/run 3 kernels (fyppm, fv_mapz, nh_core) from FV3
on the OLCF Ascent system.

# Requirements
The GNU autotools build system and a fortran compiler.  The fv_mapz and nh_core kernels
also require the netCDF fortran library.

# Configuring the environment
Since this project is targeting GPUs via OpenACC, we recommend configuring your environment
to use the PGI compiler.
```
$ source $MODULESHOME/init/bash
$ module swap xl pgi/19.5
$ module load netcdf/4.6.1
$ module load netcdf-fortran/4.4.4
```

# Building
This repository has been set up to use the GNU autotools build system.  Running
```
$ autoreconf --install
$ ./configure CC=mpicc CPPFLAGS="`pkg-config --cflags netcdf-fortran`" \
              FC=mpifort FCFLAGS="-g -acc -ta=nvidia:cc70 -Minfo=accel -Mcuda=lineinf" \
              LDFLAGS="`pkg-config --libs-only-L netcdf` `pkg-config --libs-only-L netcdf-fortran`" \
              LIBS="`pkg-config --libs-only-l netcdf` `pkg-config --libs-only-l netcdf-fortran`"
$ make
```
should produce OpenACC-enabled executables:
* gfdl_fyppm/src/fyppm.
* gfdl_mapz/src/fv_mapz
* gfdl_nh_core/src/nh_core

# Running
A bash script is provided that will submit the job to the scheduler.  It can be run:
```
$ ./submit.sh gfdl_fyppm/src/fyppm <walltime in minutes> <# of MPI ranks> <# of gpus>
```

# Checking For Correctness
### fyppm
A global sum is calculated at the end of the program and compared to a benchmark value.
If the two values agree to within a specified tolerance, the test is considered to have
run successfully and the program returns 0.  If not, then the program returns a non-zero
value.

### fv_mapz
TBD

### nh_core
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
