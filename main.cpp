#include <iostream>
#include <cstring>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
using namespace amrex;

#include <extern_parameters.H>

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  amrex::Finalize();
}
