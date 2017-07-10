#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

// This subroutine is the driver for the whole time stepping (fluid + particles )
void
mfix_level::Evolve(int lev, int nstep, int set_normg, Real dt, Real& prev_dt,
                   Real time, Real normg)
{

    if (solve_fluid)
  EvolveFluid(lev,nstep,set_normg,dt,prev_dt,time,normg);

    if (solve_dem)
    {
  if (solve_fluid)
      mfix_calc_drag_particle(lev);

  pc ->  EvolveParticles( lev, nstep, dt, time);
    }
}

void
mfix_level::EvolveFluid(int lev, int nstep, int set_normg,
                        Real dt, Real& prev_dt, Real time, Real normg)
{
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

  for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
    {
      Box domain(geom[lev].Domain());
      const Box& sbx = (*ep_g[lev])[mfi].box();
      Box ubx((*u_g[lev])[mfi].box());
      Box vbx((*v_g[lev])[mfi].box());
      Box wbx((*w_g[lev])[mfi].box());

      set_bc1(sbx.loVect(), sbx.hiVect(),
                ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
              (*u_g[lev])[mfi].dataPtr(),     (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
              bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
              bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect());
    }

  Real sum_vol;
  if (solve_dem)
    mfix_calc_volume_fraction(lev,sum_vol);

  if (nstep == 0) 
  {
     sum_vol_orig = sum_vol;
     Print() << "Setting original sum_vol to " << sum_vol_orig << std::endl;
  } else {
     Print() << "Testing new sum_vol " << sum_vol << " against original sum_vol " << sum_vol_orig << std::endl;
     if (abs(sum_vol_orig - sum_vol) > 1.e-12 * sum_vol_orig) amrex::Abort("Volume fraction in domain has changed!");
  }
  

  // Calculate transport coefficients
  int calc_flag = 2;
  mfix_calc_coeffs(lev,calc_flag);

  // Calculate the stress tensor trace and cross terms for all phases.
  mfix_calc_trd_and_tau(lev);

  // Backup field variable to old
  int nghost = ep_go[lev]->nGrow();
  MultiFab::Copy(*ep_go[lev],  *ep_g[lev],  0, 0, 1, nghost);
  MultiFab::Copy(*p_go[lev],   *p_g[lev],   0, 0, 1, nghost);
  MultiFab::Copy(*ro_go[lev],  *ro_g[lev],  0, 0, 1, nghost);
  MultiFab::Copy(*rop_go[lev], *rop_g[lev], 0, 0, 1, nghost);
  MultiFab::Copy(*u_go[lev],   *u_g[lev],   0, 0, 1, nghost);
  MultiFab::Copy(*v_go[lev],   *v_g[lev],   0, 0, 1, nghost);
  MultiFab::Copy(*w_go[lev],   *w_g[lev],   0, 0, 1, nghost);

  // Loop over iterate for auto time-step size adjustment
  int reiterate;
  do {
    prev_dt = dt;

    // Calculate bulk density (epg*ro_g) at cell faces
    mfix_conv_rop(lev,dt);

    // Calculate face mass fluxes
    mfix_calc_mflux(lev);

    int converged=0;
    int nit=0;          // number of iterations
    int gsmf=0;         // number of outer iterations for goal seek mass flux (GSMF)
    Real delP_MF=0.0L;  // actual GSMF pressure drop
    Real lMFlux=0.0L;   // actual GSMF mass flux
    Real resg=0.0L;     // fluid pressure residual

    // int lset_normg=1-set_normg;
    Real lnormg=normg;

    ///////////////// ---- call to iterate -------- /////////////////
    do {
      nit++;

      Real residuals[2*8];
      for (int i=0; i<2*8; ++i)
        residuals[i] = 0.0L;

      // User hooks
      for (MFIter mfi(*ep_g[lev]); mfi.isValid(); ++mfi)
        mfix_usr2();

      // Calculate transport coefficients
      calc_flag = 1;
      mfix_calc_coeffs(lev,calc_flag);

      // Calculate drag coefficient
      if (solve_dem)
        mfix_calc_drag_fluid(lev);

      // Solve momentum equations
      mfix_solve_for_vels(lev, dt, residuals);

      // Calculate transport coefficients
      mfix_physical_prop(lev,0);

      // Calculate bulk density (epg*ro_g) at cell faces
      mfix_conv_rop(lev,dt);

      // Solve the pressure correction equation
      mfix_solve_for_pp(lev,dt,lnormg,resg, residuals);

      // Apply pressure correction to all Pg, Ug, Vg, Wg
      mfix_correct_0(lev);

      // Update fluid density
      mfix_physical_prop(lev,0);

      // Calculate face mass fluxes
      mfix_calc_mflux(lev);

      // Check for convergence
      ParallelDescriptor::ReduceRealSum(residuals,16);
      converged = check_convergence(&nit, residuals);

      // Display current iteration residuals
      if ( ParallelDescriptor::IOProcessor() )
        display_resid(&time, &dt, &nit, residuals);

    } while(converged==0 && nit<max_nit);

    // Adjust time step if iteration failed.
    reiterate = adjustdt(&converged, &nit, &dt);
    if(reiterate == 1) {

      // Reset the field variables
      MultiFab::Copy(*ep_g[lev],  *ep_go[lev],  0, 0, 1, nghost);
      MultiFab::Copy(*p_g[lev],   *p_go[lev],   0, 0, 1, nghost);
      MultiFab::Copy(*ro_g[lev],  *ro_go[lev],  0, 0, 1, nghost);
      MultiFab::Copy(*rop_g[lev], *rop_go[lev], 0, 0, 1, nghost);
      MultiFab::Copy(*u_g[lev],   *u_go[lev],   0, 0, 1, nghost);
      MultiFab::Copy(*v_g[lev],   *v_go[lev],   0, 0, 1, nghost);
      MultiFab::Copy(*w_g[lev],   *w_go[lev],   0, 0, 1, nghost);

    }
  } while (reiterate==1);
}
