#include <fstream>
#include <iomanip>

#include <ParmParse.H>
#include <Geometry.H>
#include <VisMF.H>
#include <iMultiFab.H>

#include <mfix_level.H>
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
  int max_nit;
  int solve_fluid;
  int solve_dem;
  int steady_state;
  Real dt, dt_min, dt_max, tstop, time;
  int nstep=0;  // Number of time steps
  Real normg;
  int set_normg;
  int cyclic_x, cyclic_y, cyclic_z, cyclic_mf;

  mfix_get_data(
    &imax, &jmax, &kmax,
    &solve_fluid, 
    &solve_dem,
    &steady_state,
    &dt, &dt_min, &dt_max, &tstop, &time, &max_nit,
    &normg, &set_normg,
    &cyclic_x, &cyclic_y, &cyclic_z, &cyclic_mf);

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

  int max_level = 0;
  int lev = 0;
  Array<int> n_cell(3);
  n_cell[0] = imax;
  n_cell[1] = jmax;
  n_cell[2] = kmax;

  const RealBox* rb_ptr = &real_box;

  // Note that the constructor constructs the Geometry object now.
  mfix_level my_mfix(rb_ptr,max_level,n_cell,coord);

  my_mfix.Init(solve_fluid,solve_dem,steady_state,cyclic_mf,max_nit);

  // define dx[]
  // const Real* dx = mfix_level.geom.CellSize();

  my_mfix.call_main(lev,nstep,dt,time);

  int finish, estatus;
  finish = 0;
  estatus = 0;

  Real prev_dt; // Actual dt used to solve fluid

  my_mfix.output(lev,estatus,finish,nstep,dt,time);

  my_mfix.evolve(0,nstep,estatus,finish,set_normg,
                 dt,dt_min,dt_max,tstop,time,normg);

  my_mfix.usr3(0);

  Real end_time = ParallelDescriptor::second() - strt_time;

  if (ParallelDescriptor::IOProcessor())
    std::cout << "Time spent in main " << end_time << std::endl;

  BoxLib::Finalize();
  return 0;
}
