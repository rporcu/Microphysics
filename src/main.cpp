#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include <mfix_level.H>
#include <mfix_F.H>

int main (int argc, char* argv[])
{
  // Issue an error if AMR input file is not given
  if ( argc < 2 )
    amrex::Abort("AMReX input file missing");

  // BoxLib will now read the inputs file and the command line arguments, but the
  //        command line arguments are in mfix-format so it will just ignore them.
  amrex::Initialize(argc,argv);

  // Copy arguments into MFIX -- note that the first argument is now the name of the
  //      inputs file to be read by BoxLib, so we only pass the arguments after that
  for(int i=2; i < argc; i++) {
    int nlen = strlen(argv[i]);
    // If-statement avoids passing the name of the mfix input file if it is
    // specified on the command line or any AMReX command.
    if ( strstr(argv[i], "input_file") == NULL && strstr(argv[i], "amr") == NULL) {
      mfix_add_argument(argv[i], &nlen);
    }
  }

  Real strt_time = ParallelDescriptor::second();

  int max_nit;
  int solve_fluid;
  int solve_dem;
  int steady_state;
  int call_udf;
  Real dt, dt_min, dt_max, tstop;
  Real time=0.0L;
  Real xlength, ylength, zlength;
  int nstep=0;  // Number of time steps
  Real normg;
  int set_normg;
  int coord;
  int cyclic_x, cyclic_y, cyclic_z, cyclic_mf;

  mfix_get_data( &solve_fluid,
                 &solve_dem,
                 &steady_state,
                 &dt, &dt_min, &dt_max, &tstop, &max_nit,
                 &normg, &set_normg, &call_udf,
                 &cyclic_x, &cyclic_y, &cyclic_z, &cyclic_mf,
                 &xlength, &ylength, &zlength, &coord);

#if 0
  // This sets the boundary conditions to be doubly or triply periodic
  Array<int> is_per(3);
  is_per[0] = cyclic_x;
  is_per[1] = cyclic_y;
  is_per[2] = cyclic_z;

  // This is equivalent to reading "geometry.is_per = ..." from the inputs file,
  // but here we read cyclic_* from the mfix input file and initialize the ParmParse
  // table here so when we call Geometry::Setup from the AmrCore constructor it will see it.
  ParmParse pp("geometry");
  pp.addarr("is_periodic",is_per);
#endif

  int lev = 0;

  // Note that the constructor constructs the Geometry object now.
  mfix_level my_mfix;
//mfix_level my_mfix(rb_ptr,max_level,n_cell,coord);

  my_mfix.InitParams(solve_fluid,solve_dem,cyclic_mf,max_nit,call_udf);

  my_mfix.Init(lev,dt,time);

  // Restart from checkpoint if needed
  if ( my_mfix.IsRestartEnabled() ) {
    if ( steady_state ) {
      amrex::Warning("Restart from checkpoint enabled for "
                     "steady state solve: ignoring");
    }
    else {
      my_mfix.InitFromCheckpoint( &nstep, &dt, &time);
    }
  }

  int finish  = 0;
  int estatus = 0;

  // Call to output before entering time march loop
  if (solve_fluid && ParallelDescriptor::IOProcessor())
    my_mfix.output(lev,estatus,finish,nstep,dt,time);

  // Initialize prev_dt here; it will be re-defined by call to evolve_fluid but
  // only if solve_fluid = T
  Real prev_dt = dt;

  while (finish == 0)
  {
    mfix_usr1();

    if (solve_fluid)
      my_mfix.evolve_fluid(lev,nstep,set_normg,dt,prev_dt,time,normg);

    if (solve_dem)
      my_mfix.evolve_dem(lev,nstep,dt,time);

    if (!steady_state)  {
      time += prev_dt;
      nstep++;
      my_mfix.WritePlotFile( nstep, dt, time );
      my_mfix.WriteCheckPointFile( nstep, dt, time );
    }

    if (ParallelDescriptor::IOProcessor())
       my_mfix.output(lev,estatus,finish,nstep,dt,time);

    // Mechanism to terminate MFIX normally.
    if (steady_state || (time + 0.1*dt >= tstop) || (solve_dem && !solve_fluid)) finish = 1;
  }

  // Dump plotfile if enabled for steady state
  if (steady_state) my_mfix.WritePlotFile();

  my_mfix.usr3(0);

  Real end_time = ParallelDescriptor::second() - strt_time;

  if (ParallelDescriptor::IOProcessor())
    std::cout << "Time spent in main " << end_time << std::endl;

  amrex::Finalize();
  return 0;
}
