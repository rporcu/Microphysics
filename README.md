# Directory overview

| File       | Description                                         |
| ---------  | --------------------------------------------------- |
| benchmarks | UC Benchmark cases (see benchmark/README.md)        |
| src        | Fortran source files                                |
| tests      | Regression tests (see tests/README.md)              |
| ThirdParty | External libraries sources                          |
| Tools      | CMake configuration files                           |


# Building MFIX-Exa

Refer to the **[MFiX-Exa user guide](https://amrex-codes.github.io/MFIX-Exa/docs_html/GettingStarted.html)** for instructions on how to build MFiX-Exa.


# Running MFIX Test Suite
MFIX-Exa comes with several tests aimed at evaluating software functionalities.
The source files as well as the required input files for each test are located in
the `tests` directory. The `tests` directory is copied to the build directory
during MFIX-Exa configuration process. When a test is run (see below), output
files are stored in `build_dir/tests/test-name`.

There are various dependencies for comparing test results.

o To compare results to archived flow slices stored in `AUTOTEST`
  directories with the case files, the environment variable `FEXTRACT`
  must point to the location of the AMReX `fextract` utility located
  in the directory, `amrex/Tools/PostProcessing/F_Src`. Additionally,
  `numdiff` must be installed for comparing results.

## Run all tests
```shell
> cd to mfix-build-dir
> ctest
```
## List all tests (without running them)
```shell
> cd to mfix-build-dir
> ctest -N
```
## Run a particular test by the index listed in ctest -N
```shell
> cd to mfix-build-dir
> ctest -I 3,3             # run the third test
```
## Run a particular test by name
```shell
> cd to mfix-build-dir
> ctest -R DEM01  # running all tests with "DEM01" in the test name
```
## Run a particular test via make
```shell
> cd to mfix-build-dir
> make run_DEM01-x  # running "DEM01-x" and output to the screen
```

## Run specific

If the environment variable GRID is defined, it specifies which grid types to run for the test(s).
If GRID variable is not defined, the default is to run the tests for all grid types.
> env GRID="tiled" ctest -R DEM01  # running all tests with "DEM01" for tiled grid
> env GRID="single multiple" ctest -R DEM01  # running all tests with "DEM01" for single grid and multiple grid
> ctest -R DEM01  # running all tests with "DEM01" for all grid types (single, multiple, tiled)

## Run a user-defined case
```shell
> ./mfix inputs-myrun
```
_inputs-myrun_ is a text file containing the AMReX input parameters; this can be named anything
as long as it is the __first__ argument following the executable.  Note that many of the problem
parameters are still defined in _mfix.dat_.

# See the User's Guide for more about MFIX-Exa

To build the User's Guide,

```shell
> cd doc/UsersGuide
> make
```

This will build a PDF of the MFIX-Exa User's Guide, which contains information
about the equations being solved, run-time parameters, checkpoint/restart
capability, options for visualization and more.

# Regression testing

## Generating local regression test data

Developers are encouraged to create local benchmark data for
regression testing prior to committing code the GitLab repository.

1. Create a location to store benchmark data and clone the MFIX
   and AMReX repositories.
``` shell
mkdir /home/user/exa-rt
mkdir /home/user/exa-rt/benchmark
cd /home/user/exa-rt
git clone http://mfix.netl.doe.gov/gitlab/exa/mfix.git
git clone https://github.com/AMReX-Codes/amrex.git
git clone https://github.com/AMReX-Codes/regression_testing.git
```

2. Create a local copy the regression test setup file from the
MFIX repository.
```shell
cp mfix/tools/MFIX-tests.ini MFIX-local.ini
```

3. Edit the local setup file.
Specify the top level directories for test and web output
under the `[main]` heading.
```shell
[main]
testTopDir = /home/user/exa-rt/benchmark
webTopDir =  /home/user/exa-rt/web
```
Specify the AMReX source directory location under the `[AMReX]` heading.
```shell
[AMReX]
dir = /home/user/exa-rt/amrex
branch = development
cmakeSetupOpts = -DENABLE_AMRDATA=ON -DENABLE_MPI=ON -DENABLE_OMP=ON -DENABLE_PARTICLES=ON -DENABLE_EB=ON 
```
Specify the MFIX-Exa source directory location under the `[source]` heading.
```shell
[source]
dir = /home/user/exa-rt/mfix
branch = develop
cmakeSetupOpts = -DAMREX_INSTALL_DIR=/home/user/exa-rt/amrex/installdir
```
4. Run the AMReX regression test tool. The second argument is a user
supplied comment.
```shell
cd /home/user/exa-rt
regression_testing/regtest.py --make_benchmarks "MFIX" MFIX-local.ini
```


--------------------------------------------------------------------
## Prerequisite: Environment Dependencies on Joule (Joule specific)
For the Joule environment, load the gnu module and set environment variables first. If not on Joule, skip this step.
```shell
> module load gnu/6.1.0
> export CC=/nfs/apps/Compilers/GNU/6.1.0/bin/gcc
> export CXX=/nfs/apps/Compilers/GNU/6.1.0/bin/g++
> export F77=/nfs/apps/Compilers/GNU/6.1.0/bin/gfortran
> export FC=/nfs/apps/Compilers/GNU/6.1.0/bin/gfortran
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
