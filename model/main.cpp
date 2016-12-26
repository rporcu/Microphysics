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
  int imax,jmax,kmax,dem_solids;
  mfix_get_data(&imax,&jmax,&kmax,&dem_solids);

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

  // Matrix and rhs vector
  MultiFab A_m(ba,7,nghost);
  MultiFab b_m(ba,1,nghost);
  A_m.setVal(0.);
  b_m.setVal(0.);

  // Void fraction
  MultiFab ep_g (ba,1,nghost);
  MultiFab ep_go(ba,1,nghost);
  ep_g.setVal(0.);
  ep_go.setVal(0.);

  // Gas pressure fraction
  MultiFab p_g (ba,1,nghost);
  MultiFab p_go(ba,1,nghost);
  p_g.setVal(0.);
  p_go.setVal(0.);

  // Gas density
  MultiFab ro_g (ba,1,nghost);
  MultiFab ro_go(ba,1,nghost);
  ro_g.setVal(0.);
  ro_go.setVal(0.);

  // Gas bulk density
  MultiFab rop_g (ba,1,nghost);
  MultiFab rop_go(ba,1,nghost);
  rop_g.setVal(0.);
  rop_go.setVal(0.);

  // X-axis gas velocity
  MultiFab u_g (ba,1,nghost);
  MultiFab u_go(ba,1,nghost);
  u_g.setVal(0.);
  u_go.setVal(0.);

  // Y-axis gas velocity
  MultiFab v_g (ba,1,nghost);
  MultiFab v_go(ba,1,nghost);
  v_g.setVal(0.);
  v_go.setVal(0.);

  // Z-axis gas velocity
  MultiFab w_g (ba,1,nghost);
  MultiFab w_go(ba,1,nghost);
  w_g.setVal(0.);
  w_go.setVal(0.);

  // Pressure correction equation
  MultiFab pp_g(ba,1,nghost);
  MultiFab d_e (ba,1,nghost);
  MultiFab d_n (ba,1,nghost);
  MultiFab d_t (ba,1,nghost);
  pp_g.setVal(0.);
  d_e.setVal(0.);
  d_n.setVal(0.);
  d_t.setVal(0.);

  // Molecular viscosity
  MultiFab mu_g (ba,1,nghost);
  mu_g.setVal(0.);

  //
  MultiFab lambda_g (ba,1,nghost);
  MultiFab trD_g (ba,1,nghost);
  lambda_g.setVal(0.);
  trD_g.setVal(0.);

  //
  MultiFab tau_u_g (ba,1,nghost);
  MultiFab tau_v_g (ba,1,nghost);
  MultiFab tau_w_g (ba,1,nghost);
  tau_u_g.setVal(0.);
  tau_v_g.setVal(0.);
  tau_w_g.setVal(0.);

  MultiFab flux_gE (ba,1,nghost);
  MultiFab flux_gN (ba,1,nghost);
  MultiFab flux_gT (ba,1,nghost);
  flux_gE.setVal(0.);
  flux_gN.setVal(0.);
  flux_gT.setVal(0.);

  MultiFab rop_gE (ba,1,nghost);
  MultiFab rop_gN (ba,1,nghost);
  MultiFab rop_gT (ba,1,nghost);
  rop_gE.setVal(0.);
  rop_gN.setVal(0.);
  rop_gT.setVal(0.);

  //
  MultiFab f_gds  (ba,1,nghost);
  MultiFab drag_bm(ba,3,nghost);
  f_gds.setVal(0.);
  drag_bm.setVal(0.);
  int nparticles = 5000;

  int mmax = 1;

  Array<int> particle_state (  nparticles);
  Array<int> particle_phase (  nparticles);

  Array<Real> des_radius    (  nparticles);
  Array<Real> ro_sol        (  nparticles);
  Array<Real> pvol          (  nparticles);
  Array<Real> pmass         (  nparticles);
  Array<Real> omoi          (  nparticles);
  Array<Real> des_pos_new   (3*nparticles);
  Array<Real> des_vel_new   (3*nparticles);
  Array<Real> des_usr_var   (  nparticles);
  Array<Real> omega_new     (3*nparticles);
  Array<Real> des_acc_old   (3*nparticles);
  Array<Real> rot_acc_old   (3*nparticles);
  Array<Real> drag_fc       (3*nparticles);
  Array<Real> fc            (3*nparticles);
  Array<Real> tow           (3*nparticles);
  Array<int> pairs          (12*nparticles);

  for (MFIter mfi(flag); mfi.isValid(); ++mfi)
     mfix_MAIN(
               u_g[mfi].dataPtr(),     v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
               u_go[mfi].dataPtr(),    v_go[mfi].dataPtr(),     w_go[mfi].dataPtr(),
               p_g[mfi].dataPtr(),     p_go[mfi].dataPtr(),     pp_g[mfi].dataPtr(),
               ep_g[mfi].dataPtr(),    ep_go[mfi].dataPtr(),
               ro_g[mfi].dataPtr(),    ro_go[mfi].dataPtr(),
               rop_g[mfi].dataPtr(),   rop_go[mfi].dataPtr(),
               rop_gE[mfi].dataPtr(),  rop_gN[mfi].dataPtr(),   rop_gT[mfi].dataPtr(),
               d_e[mfi].dataPtr(),     d_n[mfi].dataPtr(),      d_t[mfi].dataPtr(),
               tau_u_g[mfi].dataPtr(), tau_v_g[mfi].dataPtr(),  tau_w_g[mfi].dataPtr(),
               flux_gE[mfi].dataPtr(), flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
               trD_g[mfi].dataPtr(),   lambda_g[mfi].dataPtr(), mu_g[mfi].dataPtr(),
               f_gds[mfi].dataPtr(),   A_m[mfi].dataPtr(),      b_m[mfi].dataPtr(),
               drag_bm[mfi].dataPtr(),
               flag[mfi].dataPtr(),
               particle_state.dataPtr(),
               particle_phase.dataPtr(), des_radius.dataPtr(), ro_sol.dataPtr(),
               pvol.dataPtr(), pmass.dataPtr(), omoi.dataPtr(),
               des_pos_new.dataPtr(), des_vel_new.dataPtr(),
               des_usr_var.dataPtr(), omega_new.dataPtr(), des_acc_old.dataPtr(),
               rot_acc_old.dataPtr(), drag_fc.dataPtr(), fc.dataPtr(), tow.dataPtr());

  for (MFIter mfi(flag); mfi.isValid(); ++mfi)
     mfix_time_march(
               u_g[mfi].dataPtr(),     v_g[mfi].dataPtr(),      w_g[mfi].dataPtr(),
               u_go[mfi].dataPtr(),    v_go[mfi].dataPtr(),     w_go[mfi].dataPtr(),
               p_g[mfi].dataPtr(),     p_go[mfi].dataPtr(),     pp_g[mfi].dataPtr(),
               ep_g[mfi].dataPtr(),    ep_go[mfi].dataPtr(),
               ro_g[mfi].dataPtr(),    ro_go[mfi].dataPtr(),
               rop_g[mfi].dataPtr(),   rop_go[mfi].dataPtr(),
               rop_gE[mfi].dataPtr(),  rop_gN[mfi].dataPtr(),   rop_gT[mfi].dataPtr(),
               d_e[mfi].dataPtr(),     d_n[mfi].dataPtr(),      d_t[mfi].dataPtr(),
               tau_u_g[mfi].dataPtr(), tau_v_g[mfi].dataPtr(),  tau_w_g[mfi].dataPtr(),
               flux_gE[mfi].dataPtr(), flux_gN[mfi].dataPtr(),  flux_gT[mfi].dataPtr(),
               trD_g[mfi].dataPtr(),   lambda_g[mfi].dataPtr(), mu_g[mfi].dataPtr(),
               f_gds[mfi].dataPtr(),   A_m[mfi].dataPtr(),      b_m[mfi].dataPtr(),
               drag_bm[mfi].dataPtr(),
               flag[mfi].dataPtr(),
               particle_state.dataPtr(), particle_phase.dataPtr(),
               des_radius.dataPtr(), ro_sol.dataPtr(),
               pvol.dataPtr(), pmass.dataPtr(), omoi.dataPtr(),
               des_pos_new.dataPtr(), des_vel_new.dataPtr(),
               des_usr_var.dataPtr(), omega_new.dataPtr(), des_acc_old.dataPtr(),
               rot_acc_old.dataPtr(), drag_fc.dataPtr(), fc.dataPtr(), tow.dataPtr(), pairs.dataPtr());

  for (MFIter mfi(flag); mfi.isValid(); ++mfi)
     mfix_usr3(u_g[mfi].dataPtr(),    v_g[mfi].dataPtr(),
               w_g[mfi].dataPtr(),    p_g[mfi].dataPtr());

  Real end_time = ParallelDescriptor::second() - strt_time;

  if (ParallelDescriptor::IOProcessor())
    std::cout << "Time spent in main " << end_time << std::endl;

  BoxLib::Finalize();
  return 0;
}
