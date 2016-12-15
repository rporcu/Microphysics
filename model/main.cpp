#include <fstream>
#include <iomanip>

#include <ParmParse.H>
#include <Geometry.H>
#include <VisMF.H>
#include <iMultiFab.H>

#include <mfix_F.H>

int main (int argc, char* argv[])
{

  // Copy arguments into MFIX
  for(int i=1; i < argc; i++) {
    int nlen = strlen(argv[i]);
    mfix_add_argument(argv[i], &nlen);
  }

  BoxLib::Initialize(argc,argv);

  Real strt_time = ParallelDescriptor::second();

  // Size of the entire domain
  int imax,jmax,kmax;
  mfix_get_data(&imax,&jmax,&kmax);

  IntVect dom_lo(IntVect(D_DECL(0,0,0)));
  IntVect dom_hi(IntVect(D_DECL(imax-1, jmax-1, kmax-1)));

  Box domain(dom_lo, dom_hi);

  // Initialize the boxarray "ba" from the single box "bx"
  BoxArray ba;
  ba.define(domain);

  int max_grid_size = 1024;

  // Break up boxarray "ba" into chunks no larger than
  // "max_grid_size" along a direction.
  ba.maxSize(max_grid_size);

  // This defines the physical size of the box.
  // Right now the box is [0,1] in each direction.
  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    real_box.setLo(n, 0.0);
    real_box.setHi(n, 1.0);
  }

  // This says we are using Cartesian coordinates
  int coord = 0;

  // This sets the boundary conditions to be doubly or triply periodic
  int is_periodic[BL_SPACEDIM];
  for (int i = 0; i < BL_SPACEDIM; i++) {
    is_periodic[i] = 0;
  }

  // This defines a Geometry object
  Geometry geom;
  geom.define(domain,&real_box,coord,is_periodic);

  // define dx[]
  const Real* dx = geom.CellSize();

  int nghost;
  if (ParallelDescriptor::NProcs() == 1) {
     nghost = 1;
  } else {
     nghost = 2;
  }

  // Define and allocate the integer MultiFab on BoxArray ba with 4 components and nghost ghost cells.
  iMultiFab flag(ba,4,nghost);
  flag.setVal(0);

  // Call set_domain for each subdomain
  // Read input data, check data, do computations for IC and BC locations
  // and flows, and set geometry parameters such as X, X_E, DToDX, etc.
  for (MFIter mfi(flag); mfi.isValid(); ++mfi)
     mfix_set_domain(flag[mfi].dataPtr());

  // Total volume of cell's DES stencil neighbors
  MultiFab vol_surr(ba,1,nghost);
  vol_surr.setVal(0.);

  // Matrix and rhs vector
  MultiFab A_m(ba,7,nghost);
  MultiFab b_m(ba,1,nghost);
  A_m.setVal(0.);
  b_m.setVal(0.);

  for (MFIter mfi(flag); mfi.isValid(); ++mfi)
     mfix_MAIN(flag[mfi].dataPtr(),vol_surr[mfi].dataPtr(),
               A_m[mfi].dataPtr(),      b_m[mfi].dataPtr());

  Real end_time = ParallelDescriptor::second() - strt_time;

  if (ParallelDescriptor::IOProcessor())
    std::cout << "Time spent in main " << end_time << std::endl;

  BoxLib::Finalize();
  return 0;
}
