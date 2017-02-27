# Building and running Test cases for MFIX-Exa

## Directory overview

| File       | Description                                         |
| ---------  | --------------------------------------------------- |
| benchmarks | UC Benchmark cases (see benchmark/README.md)        |
| model      | Fortran source files                                |
| tests      | Regression tests (see tests/README.md)              |



## Prerequisite: Environment Dependencies on Joule (Joule specific)
For the Joule environment, load the gnu module and set environment variables first. If not on Joule, skip this step.
```shell
> module load gnu/6.1.0
> export CC=/nfs/apps/Compilers/GNU/6.1.0/bin/gcc
> export CXX=/nfs/apps/Compilers/GNU/6.1.0/bin/g++
> export F77=/nfs/apps/Compilers/GNU/6.1.0/bin/gfortran
> export FC=/nfs/apps/Compilers/GNU/6.1.0/bin/gfortran
```

## Prerequisite: Build and Install AMReX

Clone AMReX repo and set environment variable BOXLIB_HOME to where your want to
install AMReX. If that location is not the AMReX git repo directory, manually
copy the Tools/ directory to BOXLIB_HOME.

```shell
> git clone https://bitbucket.org/berkeleylab/amrex.git
> cd amrex
> git checkout master
> export BOXLIB_HOME=$HOME/local
> cmake -DENABLE_MPI=0 -DBL_USE_PARTICLES=1 -DCMAKE_INSTALL_PREFIX:PATH=$BOXLIB_HOME .
> cp -r Tools $BOXLIB_HOME
> make -j -k install
```
The make command may fail with an error involving mempool; if so rerun ```make -j -k install``` until it succeeds.


Go to the MFIX directory, run cmake and run make (make sure BOXLIB_HOME is still set)
```shell
> cd <MFIX source directory>
> cmake .
> make
```

---------------------------------------------------------------------

## Building MFIX-Exa

To build MFIX-Exa, run cmake then make.

### Building with defaults
```shell
> cd <path to>/mfix
> cmake .
> make -j
```

### Building with user-defined files (UDFs)

To build MFIX-Exa for case with case-specific user-defined source files (UDFs),
run cmake from the directory for that case. The following examples use DEM01 as the example
case.

```shell
> cd <path to>/mfix
> cd tests/DEM01
> ls
mfix.dat usr.f90 usr_mof.f90
> cmake ../..
> make -j
```

### Building with compiler flags

You can specify compiler options with CMAKE_Fortran_FLAGS and CMAKE_CXX_FLAGS,
either specified on the command line or by editing them in CMakeCache.txt.

```shell
> cmake -DCMAKE_Fortran_FLAGS="-O2 -g -ffpe-trap=invalid -fimplicit-none" -DCMAKE_CXX_FLAGS="-std=c++11" ..
```

### Building with SMP or DMP

Specify DMP or SMP to building with MPI or OpenMP support.

```shell
> cmake -DDMP=1 ../..   # for distributed memory (MPI) executable
> cmake -DSMP=1 ../..   # for shared memory (OpenMP) executable
```

## Running MFIX Test Suite

Running tests requires the `numdiff` command, which can be installed with `apt
install numdiff` on Ubuntu.

### Running all tests
```shell
> ctest
```

### Listing all tests (without running them)
```shell
> ctest -N
Test project ~/mfixexa
  Test #1: tests_FLD01
  Test #2: tests_FLD02
  Test #3: tests_DEM01
  Test #4: tests_DEM02
  Test #5: tests_DEM03
  Test #6: tests_DEM04
  Test #7: tests_DEM05
  Test #8: tests_DEM06

Total Tests: 8
```

### Running a particular test by the index listed in ctest -N
```shell
> ctest -I 3             # run the third test
Test project ~/mfixexa
    Start 3: tests_DEM01
```

### Running a particular test by name
```shell
> ctest -R DEM01  # running all tests with "DEM01" in the test name
```


--------------------------------------------------------------------

## Notice
Neither the United States Government nor any agency thereof, nor any
of their employees, makes any warranty, expressed or implied, or
assumes any legal liability or responsibility for the accuracy,
completeness, or usefulness of any information, apparatus, product,
or process disclosed or represents that its use would not infringe
privately owned rights.

* MFIX is provided without any user support for applications in the
  user's immediate organization. It should not be redistributed in
  whole or in part.

* The use of MFIX is to be acknowledged in any published paper based
  on computations using this software by citing the MFIX theory
  manual. Some of the submodels are being developed by researchers
  outside of NETL. The use of such submodels is to be acknowledged
  by citing the appropriate papers of the developers of the submodels.

* The authors would appreciate receiving any reports of bugs or other
  difficulties with the software, enhancements to the software, and
  accounts of practical applications of this software.

# Disclaimer
This report was prepared as an account of work sponsored by an agency
of the United States Government. Neither the United States Government
nor any agency thereof, nor any of their employees, makes any
warranty, express or implied, or assumes any legal liability or
responsibility for the accuracy, completeness, or usefulness of any
information, apparatus, product, or process disclosed, or represents
that its use would not infringe privately owned rights. Reference
herein to any specific commercial product, process, or service by
trade name, trademark, manufacturer, or otherwise does not
necessarily constitute or imply its endorsement, recommendation, or
favoring by the United States Government or any agency thereof. The
views and opinions of authors expressed herein do not necessarily
state or reflect those of the United States Government or any
agency thereof.
