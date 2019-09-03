
namespace amrex {

const char* buildInfoGetBuildDate() {

  static const char BUILD_DATE[] = "2019-09-03 11:19:41.347530";
  return BUILD_DATE;
}

const char* buildInfoGetBuildDir() {

  static const char BUILD_DIR[] = "/.nfs/home/2/jmusser/mfix/exa/exec";
  return BUILD_DIR;
}

const char* buildInfoGetBuildMachine() {

  static const char BUILD_MACHINE[] = "Linux login02 3.10.0-693.el7.x86_64 #1 SMP Tue Aug 22 21:09:27 UTC 2017 x86_64 x86_64 x86_64 GNU/Linux";
  return BUILD_MACHINE;
}

const char* buildInfoGetAMReXDir() {

  static const char AMREX_DIR[] = "../../amrex";
  return AMREX_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "gnu";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "6.5.0";
  return COMP_VERSION;
}

// deprecated
const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "";
  return FCOMP;
}

// deprecated
const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "";
  return FCOMP_VERSION;
}

const char* buildInfoGetCXXName() {

  static const char CXX_comp_name[] = "mpicxx";
  return CXX_comp_name;
}

const char* buildInfoGetFName() {

  static const char F_comp_name[] = "mpif90";
  return F_comp_name;
}

const char* buildInfoGetCXXFlags() {

  static const char CXX_flags[] = " -Werror=return-type -g -O3  -pthread  -DNDEBUG -DBL_USE_MPI -DAMREX_USE_MPI -DAMREX_GIT_VERSION=\"19.08-343-gdcfcc8913e79\" -DBL_GCC_VERSION=6.5.0 -DBL_GCC_MAJOR_VERSION=6 -DBL_GCC_MINOR_VERSION=5 -DAMREX_LAUNCH= -DAMREX_DEVICE= -DAMREX_CUDA_FORT_GLOBAL= -DAMREX_CUDA_FORT_DEVICE= -DAMREX_CUDA_FORT_HOST= -DAMREX_CUDA_FORT_HOST_DEVICE= -DBL_SPACEDIM=3 -DAMREX_SPACEDIM=3 -DBL_FORT_USE_UNDERSCORE -DAMREX_FORT_USE_UNDERSCORE -DBL_Linux -DAMREX_Linux -DAMREX_PARTICLES -DAMREX_USE_EB -I. -I../src -I../src/des_fluid -I../src/io -I../src/convection -I../src/diffusion -I../src/projection -I../src/setup -I../src/util -I../src/check_data -I../src/mods -I../src/src_des -I../src/eb -I../src/usr -I../src/include -I../../amrex/Src/Base -I../../amrex/Src/AmrCore -I../../amrex/Src/Boundary -I../../amrex/Src/Particle -I../../amrex/Src/EB -I../../amrex/Src/Base -I../../amrex/Src/AmrCore -I../../amrex/Src/Boundary -I../../amrex/Src/Particle -I../../amrex/Src/EB -I../../amrex/Src/LinearSolvers/MLMG -I../../amrex/Src/LinearSolvers/MLMG -I../../amrex/Tools/C_scripts";
  return CXX_flags;
}

const char* buildInfoGetFFlags() {

  static const char F_flags[] = " -g -O3 -ffree-line-length-none -fno-range-check -fno-second-underscore -fimplicit-none  ";
  return F_flags;
}

const char* buildInfoGetLinkFlags() {

  static const char link_flags[] = "-L. -L/.nfs/misc/apps/Compilers/GNU/6.5.0/bin/../lib/gcc/x86_64-pc-linux-gnu/6.5.0/../../../../lib64/";
  return link_flags;
}

const char* buildInfoGetLibraries() {

  static const char libraries[] = "-pthread -I/nfs/home/2/jmusser/packages/openmpi/3.1.4_gnu6.5-cuda10/lib -Wl,-rpath -Wl,/nfs/home/2/jmusser/packages/openmpi/3.1.4_gnu6.5-cuda10/lib -Wl,--enable-new-dtags -L/nfs/home/2/jmusser/packages/openmpi/3.1.4_gnu6.5-cuda10/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lquadmath";
  return libraries;
}

const char* buildInfoGetAux(int i) {

  //static const char AUX1[] = "${AUX[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

int buildInfoGetNumModules() {
  // int num_modules = X;
  int num_modules = 0;

  return num_modules;
}

const char* buildInfoGetModuleName(int i) {

  //static const char MNAME1[] = "${MNAME[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetModuleVal(int i) {

  //static const char MVAL1[] = "${MVAL[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetGitHash(int i) {

  //static const char HASH1[] = "${GIT[1]}";
  static const char HASH1[] = "18.12.1-1086-gc2c5786";
  static const char HASH2[] = "19.08-343-gdcfcc89";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return HASH1;
    case 2: return HASH2;

    default: return EMPT;
  }
}

const char* buildInfoGetBuildGitHash() {

  //static const char HASH[] = "${GIT}";
  static const char HASH[] = "";


  return HASH;
}

const char* buildInfoGetBuildGitName() {

  //static const char NAME[] = "";
  static const char NAME[] = "";


  return NAME;
}

}
