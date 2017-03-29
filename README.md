# Building and running Test cases for MFIX-Exa

## Directory overview

| File       | Description                                         |
| ---------  | --------------------------------------------------- |
| benchmarks | UC Benchmark cases (see benchmark/README.md)        |
| src        | Fortran source files                                |
| tests      | Regression tests (see tests/README.md)              |
| ThirdParty | External libraries sources                          |
| Tools      | CMake configuration files                           |

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
Set the environment variable AMREX_HOME to point to *installdir*, an installation directory of choice. If *installdir* does not exist, the build system will create it for you. __It is strongly recommended that *installdir* not be placed in the same folder as the AMReX source files.__ 
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

## Building MFIX-Exa

To build MFIX-Exa, you have two options:

1. Use an already existing AMReX installation
1. Build and install AMReX as part of the MFIX build (SUPERBUILD,default)

### 1. Building MFIX against a pre-installed version of AMReX 
Clone the mfix git repo
```shell
> git clone http://mfix.netl.doe.gov/gitlab/exa/mfix.git
> cd mfix
``` 
Create a build directory
```shell
> mkdir build
> cd build
``` 
Make sure that AMREX_HOME points to the AMReX installation directory (see above)
```shell
> export AMREX_HOME=absolute-path-to-amrex_installdir
```
Run CMake and build
```shell
> cmake CMAKE_CONFIG_OPTIONS -DENABLE_SUPERBUILD=0 ..
> make -j
```
CMAKE_CONFIG_OPTIONS must be the same string used for the AMReX configuration.

### 2. Building MFIX via SUPERBUILD
In order to avoid a separate AMReX installation, you can take advantage of the default MFIX build configuration
by doing as follows (starting from the clone phase)
```shell
> git clone http://mfix.netl.doe.gov/gitlab/exa/mfix.git
> cd mfix
> mkdir build
> cd build
> cmake CMAKE_CONFIG_OPTIONS  ..
> make -j
```  
This will build AMReX internally the first time the *make* command is issued. No external installation of AMReX will be needed. 

## Customizing MFIX configuration
The options string CMAKE_CONFIG_OPTIONS allows to enable/disable features of both AMReX and MFIX. If AMReX is built separatly, MFIX must be configured with the same CMAKE_CONFIG_OPTIONS used for the configuration of AMReX. If instead SUPERBUILD is used, MFIX build system will take care of applying CMAKE_CONFIG_OPTIONS to the AMReX configuration.

CMAKE_CONFIG_OPTIONS consists of a series of options of the form -D _option-name_=*option_value*. The table below lists all the possible options, their possible values and a description of their effect on the build.
 
| Option name          |  Description                                 | Possible values              | Default value  |
| ---------------------|----------------------------------------------|------------------------------|----------------| 
| ENABLE_MPI           | Enable build with MPI                        |   0/1                        |   0            |
| ENABLE_OpenMP        | Enable build with OpenMP                     |   0/1                        |   0            |
| ENABLE_PROFILING     | Include profiling information in AMReX build |   0/1                        |   0            |
| ENABLE_BACKTRACE     | Include backtrace information in AMReX build |   0/1                        |   1            |
| FFLAGS               | User-defined Fortran flags                   | all compiler-supported flags |   None         |
| CXXFLAGS             | User-defined C++ flags                       | all compiler-supported flags |   None         |

For example, invoking cmake as follows (_CMAKE_CONFIG_OPTIONS="-DFFLAGS=-fcray-pointer -DENABLE_MPI=1"_ )
```shell
> cmake -DFFLAGS=-fcray-pointer -DENABLE_MPI=1 ..
```
adds the flag *-fcray-pointer* to the Fortran compilation command and enable MPI subroutines.  
The system defaults compilers can be overwritten as well by setting the flags FC and CXX before invoking the cmake command.
```shell
> FC=fortran-compiler CXX=c++-compiler cmake CMAKE_CONFIG_OPTIONS  ..
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

## Running MFIX Test Suite
MFIX comes  with a number of test cases aimed at testing the functionalities of the software package.
The source files as well as the required input files for each test are located in _MFIX-git-repo_/_tests_/_test-name_. 
During the MFIX configuration process, the directory *tests* is copied to the build directory.
When a test is run (see below), the output files are stored in _build-dir_/_tests_/_test-name_.

**Running tests requires the `numdiff` command, which can be installed with `apt install numdiff` on Ubuntu.**

### Running all tests
```shell
> cd to mfix-build-dir
> ctest
```

### Listing all tests (without running them)
```shell
> cd to mfix-build-dir
> ctest -N
Test project path-to-build-dir
  Test  #1: FLD01-x
  Test  #2: FLD01-y
  Test  #3: FLD01-z
  Test  #4: FLD02-x
  Test  #5: FLD02-y
  Test  #6: FLD02-z
  Test  #7: DEM01-x
  Test  #8: DEM01-y
  Test  #9: DEM01-z
  Test #10: DEM02-x
  Test #11: DEM02-y
  Test #12: DEM02-z
  Test #13: DEM03-x
  Test #14: DEM03-y
  Test #15: DEM03-z
  Test #16: DEM04-x
  Test #17: DEM04-y
  Test #18: DEM04-z
  Test #19: DEM05-x
  Test #20: DEM05-y
  Test #21: DEM05-z
  Test #22: DEM06-y
  Test #23: DEM06-z
  Test #24: DEM06-x
```

### Running a particular test by the index listed in ctest -N
```shell
> cd to mfix-build-dir
> ctest -I 3,3             # run the third test
```
### Running a particular test by name
```shell
> cd to mfix-build-dir
> ctest -R DEM01  # running all tests with "DEM01" in the test name
```
### Running a particular test via make
```shell
> cd to mfix-build-dir
> make run_DEM01-x  # running "DEM01-x" and output to the screen
```
### Running a user-defined case
```shell
> ./mfix inputs  mfix.input_file=<user_file_name>
```
_inputs_ is a text file containing the BoxLib input parameters.
_inputs_  __has to be provided and cannot be renamed__.
_user_file_name_ is the name of a user-defined text file containing the MFIX input parameters.
If _mfix.input_file=input_file_name_ is not given, MFIX will try to read the file
_mfix.dat_. MFIX __requires__ either _user_file_name_ or _mfix.dat_.

### Writing plotfiles
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

### Writing checkfiles
To dump a checkfile every N time steps, add the following to the _inputs_
file:
```shell
mfix.check_int=N
```
The name of the checkfile can be specified by adding
```shell
mfix.check_file=<checkfile_name>
```
to the _inputs_ file. In order to restart a calculation from a checkpoint,
add add the following to the _inputs_
file:
```shell
mfix.restart_chkfile=<checkfile_name>
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
