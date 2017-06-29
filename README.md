# Directory overview

| File       | Description                                         |
| ---------  | --------------------------------------------------- |
| benchmarks | UC Benchmark cases (see benchmark/README.md)        |
| src        | Fortran source files                                |
| tests      | Regression tests (see tests/README.md)              |
| ThirdParty | External libraries sources                          |
| Tools      | CMake configuration files                           |


# Building MFIX-Exa

There are two options to building MFIX-Exa:

o __SUPERBUILD (default):__  Utilities download and build AMReX as
   part of the MFIX-Exa build process. This method is strongly
   encouraged as it ensures configure options for MFIX-Exa and
   AMReX are consistent.

o __Using an existing AMReX Library:__  MFIX-Exa is linked to an
   existing AMReX installation. This is ideal for continuous
   integration severs (CI) and regression testing applications.
   AMReX library version must meet MFIX-Exa requirements.

## SUPERBUILD Instructions (recommended)

When building in __SUPERBUILD__ mode, MFIX-Exa build system will take 
care of downloading, configuring and installing AMReX as part of the
MFIX-Exa build. The following commands build MFIX-Exa and AMReX in a single step.

```shell
> git clone http://mfix.netl.doe.gov/gitlab/exa/mfix.git
> cd mfix
> mkdir build
> cd build
> cmake CONFIG_OPTIONS ..
> make -j
```
AMReX is built the first time the `make` command is issued.
No external installation of AMReX is needed. However, internet access
to the AMReX github repository is required.
The optional string `CONFIG_OPTIONS` allows to customize the build of both AMReX and MFIX-Exa.
`CONFIG_OPTIONS` is a list of one or more configuration options given in the form
__-D__*OPTION=VALUE*`.
 
The table below lists configuration options, possible values, and their effect on the build.
Options prefixed by `AMREX_` are specific to the build of AMReX.

| Option name                  |  Description                                 | Possible values              | Default value       |
| -----------------------------|----------------------------------------------|------------------------------|---------------------|
| CMAKE_BUILD_TYPE             | Configuration of the build                   |   Debug/Release              |   Debug             |
| MFIX_FFLAGS_OVERRIDES        | User-defined Fortran flags                   | valid F90 compiler flags     |   None              |
| MFIX_CXXLAGS_OVERRIDES       | User-defined C++ flags                       | valid C++ compiler flags     |   None              |
| AMREX_ENABLE_MPI             | Enable build with MPI                        |   0/1                        |   1                 |
| AMREX_ENABLE_OMP             | Enable build with OpenMP                     |   0/1                        |   0                 |
| AMREX_ENABLE_DP              | Enable double precision                      |   0/1                        |   1                 |
| AMREX_ENABLE_DP_PARTICLES    | Enable double precision in particles classes |   0/1                        |   1                 |
| AMREX_ENABLE_PROFILING       | Include profiling info                       |   0/1                        |   0                 |
| AMREX_ENABLE_TINY_PROFILING  | Include tiny profiling info                  |   0/1                        |   0                 |
| AMREX_ENABLE_BACKTRACE       | Include backtrace info                       |   0/1                        |   1                 |
| AMREX_ENABLE_PIC             | Build position-independent code              |   0/1                        |   0                 |
| AMREX_GIT_COMMIT             | AMReX commit to be used in the build         | valid git commit id/branch   |   None              |
| AMREX_INSTALL_DIR            | Global path to AMReX install directory       | valid global path            |   None (superbuild) |
 
`__SUPERBUILD__ mode is enabled automatically when _AMREX_INSTALL_DIR_ is not given.`

Example: build mfix with custom fortran flags, AMReX profiling enabled and single precision particles:

```
cmake -DMFIX_FFLAGS_OVERRIDES="custom flags" -DAMREX_ENABLE_PROFILING=1 -DAMREX_ENABLE_DP_PARTICLES=0 ..
```  

## Building MFIX-Exa using a separate AMReX installation (no superbuild)

### Build AMReX Library

Clone AMReX from the official Git repository and checkout the _development_ branch.
```shell
> git clone https://github.com/AMReX-Codes/amrex.git
> cd amrex
> git checkout development
```
Next, configure, build and install AMReX as follows: 
```shell
> cmake AMREX_CONFIG_OPTIONS -DENABLE_PARTICLES=1 -DCMAKE_INSTALL_PREFIX:PATH=/absolute/path/to/installdir .
> make install
```
Here,`AMREX_CONFIG_OPTIONS` are optional configure options for AMReX. Please refer to the AMReX user guide
for a list of all the possible configuration options. The only option required is __-DENABLE_PARTICLES=1__.

### Build MFIX-Exa
Clone and build MFIX-Exa.
```shell
> git clone http://mfix.netl.doe.gov/gitlab/exa/mfix.git
> mkdir build
> cd build
> cmake CONFIG_OPTIONS -DAMREX_INSTALL_DIR=/absolute/path/to/amrex/installdir ..
> make -j
```
Here, `CONFIG_OPTIONS` are MFIX-Exa specific configuration options, that is, any option prefixed by `MFIX_`
in the table above. Options prefixed by `AMREX_` are always ignored 
when using an external AMReX installation.


## Custom configurations (`CMAKE_CONFIG_OPTIONS`)

The system defaults compilers can be overwritten as well by setting the flags
`FC` and `CXX` before invoking the cmake command.
```shell
> FC=fortran-compiler CXX=c++-compiler cmake CMAKE_CONFIG_OPTIONS  ..
```


# Running MFIX Test Suite
MFIX-Exa comes  with several tests aimed at evaluating software functionalities.
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

o To compare point-by-point field data, the environment variable
  `FCOMPARE` must point the AMReX utility `plt_compare_diff_grids`
  found in the directory, `amrex/Tools/PostProcessing/F_Src`.
  Additionally, the environment variable `MFIX_BENCHMARKS_HOME`
  must point to the location of a local regression test data set.
  See _Generating local regression test data_ for instructions on
  creating a local regression test data set.

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
> ./mfix inputs  mfix.input_file=<user_file_name>
```
_inputs_ is a text file containing the AMReX input parameters.
_inputs_  __has to be provided and cannot be renamed__.
_user_file_name_ is the name of a user-defined text file containing the MFIX input parameters.
If _mfix.input_file=input_file_name_ is not given, MFIX will try to read the file
_mfix.dat_. MFIX __requires__ either _user_file_name_ or _mfix.dat_.


# Writing plotfiles
In order to write out plotfiles, add the following to the _inputs_ file:
```shell
amr.plot_int=N
```
N needs to be > 1 for the plotfiles to be written out. For transient solves,
N indicates the number of time steps between two consecutive writes.
For steady state solve, N does not have any meaning: a plotfile will be written
after the steady state is reached, as long as N > 0. To specify the name of the
plotfiles directories, add the following to the _inputs_ file:
```shell
amr.plot_file=<plotfile_name>
```
If the name of the plotfile is not provided, MFIX will default to _plt_.


# Writing checkfiles
To dump a checkfile every N time steps, add the following to the _inputs_
file:
```shell
amr.check_int=N
```
The name of the checkfile can be specified by adding
```shell
amr.check_file=<checkfile_name>
```
to the _inputs_ file. In order to restart a calculation from a checkpoint,
add add the following to the _inputs_
file:
```shell
amr.restart_chkfile=<checkfile_name>
```

# Writing particles ASCII files
To write an ASCII file with particles data every N time steps, add the
following to the _inputs_ file:
```shell
amr.par_ascii_file=<ascii_file_name>
amr.par_ascii_int=N
```
N needs to be > 1 for the ASCII files to be written out. For transient solves,
N indicates the number of time steps between two consecutive writes.
For steady state solve, N does not have any meaning: an ASCII file will be written
after the steady state is reached, as long as N > 0. Specifying the
ASCII file name is optional.


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
```

2. Create a local copy the regression test setup file from the
MFIX repository.
```shell
cp mfix/RegressionTesting/MFIX-tests.ini MFIX-local.ini
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
```
Specify the MFIX-Exa source directory location under the `[source]` heading.
```shell
[source]
dir = /home/user/exa-rt/mfix
branch = develop
```
4. Run the AMReX regression test tool. The second argument is a user
supplied comment.
```shell
cd /home/user/exa-rt
amrex/Tools/RegressionTesting/regtest.py --make_benchmarks "MFIX" MFIX-local.ini
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
## Prerequisite: Building and installing AMReX (not needed if SUPERBUILD is enabled: see below)
Clone AMReX from the official Git repository and checkout the _development_ branch.
```shell
> git clone https://bitbucket.org/berkeleylab/amrex.git
> cd amrex
> git checkout development
```
Set the environment variable `AMREX_HOME` to point to *installdir*, an installation
directory of choice. If *installdir* does not exist, the build system will
create it for you. __It is strongly recommended that *installdir* not
be placed in the same folder as the AMReX source files.__
```shell
> export AMREX_HOME=absolute-path-to-installdir
```
From the AMReX source directory, run CMake to build and install AMReX.
```shell
> cmake CMAKE_CONFIG_OPTIONS -DBL_USE_PARTICLES=1 -DCMAKE_INSTALL_PREFIX:PATH=$AMREX_HOME .
> make install
```
CMAKE_CONFIG_OPTIONS represents a string of CMake configuration options as explained below.

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
