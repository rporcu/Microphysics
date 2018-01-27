#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

void
mfix_level::EvolveFluidSimple(int lev, int nstep, int set_normg,
                              Real dt, Real& prev_dt, Real time, Real normg)
{
  BL_PROFILE_REGION_START("mfix::EvolveFluidSimple");

  // Reimpose boundary conditions -- make sure to do this before we compute tau
  mfix_set_bc1(lev);

    // Calculate transport coefficients
    int calc_flag = 2;
    mfix_calc_coeffs(lev,calc_flag);

    // Calculate the stress tensor trace and cross terms for all phases.
    mfix_calc_trd_and_tau(lev);

    // Backup field variable to old
    int nghost = ep_go[lev]->nGrow();

    MultiFab::Copy(*ep_go[lev],  *ep_g[lev],  0, 0, 1, nghost);
    MultiFab::Copy( *p_go[lev],   *p_g[lev],  0, 0, 1, nghost);
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
	mfix_conv_rop(lev);

	// Calculate face mass fluxes
	mfix_calc_mflux(lev);

	int converged=0;
	int nit=0;          // number of iterations

	///////////////// ---- call to iterate -------- /////////////////
	do {
	    nit++;

	    Real residuals[2*8];
	    for (int i=0; i<2*8; ++i)
		residuals[i] = 0.0L;

	    // User hooks
#ifdef _OPENMP
#pragma omp parallel
#endif
	    for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
		mfix_usr2();

	    // Calculate drag coefficient
	    if (solve_dem)
		mfix_calc_drag_fluid(lev);

	    // Solve momentum equations in a thread safe way.
	    Real num_u = 0.0L;
	    Real num_v = 0.0L;
	    Real num_w = 0.0L;

	    Real denom_u = 0.0L;
	    Real denom_v = 0.0L;
	    Real denom_w = 0.0L;

	    mfix_solve_for_u(lev, dt, num_u, denom_u);
	    mfix_solve_for_v(lev, dt, num_v, denom_v);
	    mfix_solve_for_w(lev, dt, num_w, denom_w);

	    residuals[1] = num_u;
	    residuals[2] = num_v;
	    residuals[3] = num_w;

	    residuals[9]  = denom_u;
	    residuals[10] = denom_v;
	    residuals[11] = denom_w;

	    //Called after each momentum equation variable solve
	    MultiFab::Copy(*u_g[lev], *u_gt[lev], 0, 0, 1, u_g[lev]->nGrow());
	    MultiFab::Copy(*v_g[lev], *v_gt[lev], 0, 0, 1, v_g[lev]->nGrow());
	    MultiFab::Copy(*w_g[lev], *w_gt[lev], 0, 0, 1, w_g[lev]->nGrow());

	    u_g[lev]->FillBoundary(geom[lev].periodicity());
	    v_g[lev]->FillBoundary(geom[lev].periodicity());
	    w_g[lev]->FillBoundary(geom[lev].periodicity());

	    // Calculate transport coefficients
	    mfix_physical_prop(lev,0);

	    // Reimpose boundary conditions
	    mfix_set_bc1(lev);

	    // Calculate bulk density (epg*ro_g) at cell faces
	    mfix_conv_rop(lev);

	    // Solve the pressure correction equation
	    Real   num_p = 0.0L;
	    Real denom_p = 0.0L;
	    mfix_solve_for_pp(lev,dt,num_p,denom_p);
	    residuals[0] = num_p;
	    residuals[8] = denom_p;

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

  BL_PROFILE_REGION_STOP("mfix::EvolveFluidSimple");
}

void
mfix_level::mfix_calc_trd_and_tau(int lev)
{
  BL_PROFILE("mfix_level::mfix_calc_trd_and_tau()");
  Real dx = geom[lev].CellSize(0);
  Real dy = geom[lev].CellSize(1);
  Real dz = geom[lev].CellSize(2);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(*ep_g[lev],true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      const Box& sbx = (*ep_g[lev])[mfi].box();

      Box ubx((*u_g[lev])[mfi].box());
      Box vbx((*v_g[lev])[mfi].box());
      Box wbx((*w_g[lev])[mfi].box());

      calc_trd_g(sbx.loVect(), sbx.hiVect(),
                 ubx.loVect(), ubx.hiVect(),
                 vbx.loVect(), vbx.hiVect(),
                 wbx.loVect(), wbx.hiVect(),
                 bx.loVect(),  bx.hiVect(),
                 (*trD_g[lev])[mfi].dataPtr(),
                 (*u_g[lev])[mfi].dataPtr(),
                 (*v_g[lev])[mfi].dataPtr(),
                 (*w_g[lev])[mfi].dataPtr(),
                 &dx, &dy, &dz);
    }
  fill_mf_bc(lev,*trD_g[lev]);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(*ep_g[lev],true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      const Box& sbx = (*ep_g[lev])[mfi].box();

      Box ubx((*u_g[lev])[mfi].box());
      Box vbx((*v_g[lev])[mfi].box());
      Box wbx((*w_g[lev])[mfi].box());

      calc_tau_g(sbx.loVect(), sbx.hiVect(),
                 ubx.loVect(), ubx.hiVect(),
                 vbx.loVect(), vbx.hiVect(),
                 wbx.loVect(), wbx.hiVect(),
                 bx.loVect(),  bx.hiVect(),
                 (*tau_u_g[lev])[mfi].dataPtr(),
                 (*tau_v_g[lev])[mfi].dataPtr(),
                 (*tau_w_g[lev])[mfi].dataPtr(),
                 (*u_g[lev])[mfi].dataPtr(),
                 (*v_g[lev])[mfi].dataPtr(),
                 (*w_g[lev])[mfi].dataPtr(),
                 (*trD_g[lev])[mfi].dataPtr(),
                 (*lambda_g[lev])[mfi].dataPtr(),
                 (*mu_g[lev])[mfi].dataPtr(),
                 &dx, &dy, &dz);
    }
  tau_u_g[lev]->FillBoundary(geom[lev].periodicity());
  tau_v_g[lev]->FillBoundary(geom[lev].periodicity());
  tau_w_g[lev]->FillBoundary(geom[lev].periodicity());
}

void
mfix_level::mfix_calc_mflux(int lev)
{
  BL_PROFILE("mfix_level::mfix_calc_mflux()");
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*u_g[lev], true); mfi.isValid(); ++mfi)
    {
      Box ubx((*u_g[lev])[mfi].box());
      Box vbx((*v_g[lev])[mfi].box());
      Box wbx((*w_g[lev])[mfi].box());

  calc_mflux(
      ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
      (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
      (*ropX[lev])[mfi].dataPtr(),   (*ropY[lev])[mfi].dataPtr(),   (*ropZ[lev])[mfi].dataPtr(),
      (*fluxX[lev])[mfi].dataPtr(),  (*fluxY[lev])[mfi].dataPtr(),  (*fluxZ[lev])[mfi].dataPtr(),
      &dx, &dy, &dz);
    }

    // Impose periodic bc's at domain boundaries and fine-fine copies in the interio
    fluxX[lev]->FillBoundary(geom[lev].periodicity());
    fluxY[lev]->FillBoundary(geom[lev].periodicity());
    fluxZ[lev]->FillBoundary(geom[lev].periodicity());
}

void
mfix_level::mfix_conv_rop(int lev)
{
  BL_PROFILE("mfix_level::mfix_conv_rop()");
    Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*rop_g[lev], true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        Box  sbx((*rop_g[lev])[mfi].box());
        Box  ubx((  *u_g[lev])[mfi].box());
        Box  vbx((  *v_g[lev])[mfi].box());
        Box  wbx((  *w_g[lev])[mfi].box());
        Box rxbx(( *ropX[lev])[mfi].box());
        Box rybx(( *ropY[lev])[mfi].box());
        Box rzbx(( *ropZ[lev])[mfi].box());

        conv_rop( bx.loVect(),  bx.hiVect(),
                  (*rop_g[lev])[mfi].dataPtr(),  sbx.loVect(),  sbx.hiVect(),
                  (  *u_g[lev])[mfi].dataPtr(),  ubx.loVect(),  ubx.hiVect(),
                  (  *v_g[lev])[mfi].dataPtr(),  vbx.loVect(),  vbx.hiVect(),
                  (  *w_g[lev])[mfi].dataPtr(),  wbx.loVect(),  wbx.hiVect(),
                  ( *ropX[lev])[mfi].dataPtr(), rxbx.loVect(), rxbx.hiVect(),
                  ( *ropY[lev])[mfi].dataPtr(), rybx.loVect(), rybx.hiVect(),
                  ( *ropZ[lev])[mfi].dataPtr(), rzbx.loVect(), rzbx.hiVect());
    }

    ropX[lev]->FillBoundary(geom[lev].periodicity());
    ropY[lev]->FillBoundary(geom[lev].periodicity());
    ropZ[lev]->FillBoundary(geom[lev].periodicity());
}

void
mfix_level::mfix_solve_for_u(int lev, Real dt, Real& num_u, Real& denom_u)
{
    BL_PROFILE("mfix_level::solve_for_u()");
    Box domain(geom[lev].Domain());

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // Matrix and rhs vector
    BoxArray x_edge_ba = grids[lev];
    x_edge_ba.surroundingNodes(0);
    A_m[lev].reset(new MultiFab(x_edge_ba,dmap[lev],7,0));
    b_m[lev].reset(new MultiFab(x_edge_ba,dmap[lev],1,0));

    // Solve U-Momentum equation
    MultiFab::Copy(*u_gt[lev], *u_g[lev], 0, 0, 1, u_g[lev]->nGrow());

    // Initialize d_e to 0
    d_e[lev]->setVal(0.);

    auto mask = u_g[lev]->OverlapMask(geom[lev].periodicity());

    // We make temporaries because the omp parallel reduction won't take reference values.
    Real temp_num = num_u;
    Real temp_denom = denom_u;


#ifdef _OPENMP
#pragma omp parallel reduction(+:temp_num,temp_denom)
#endif
    for (MFIter mfi(*u_g[lev],true); mfi.isValid(); ++mfi)
    {

       Real wt = ParallelDescriptor::second();

       const Box& bx = mfi.tilebox();
       const Box& sbx = (*ep_g[lev])[mfi].box();
       Box abx((*A_m[lev])[mfi].box());

       Box ubx((*u_g[lev])[mfi].box());
       Box vbx((*v_g[lev])[mfi].box());
       Box wbx((*w_g[lev])[mfi].box());

       Box dbx((*drag_u[lev])[mfi].box());

       solve_u_g_star(sbx.loVect(), sbx.hiVect(),
           ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(),
           wbx.loVect(), wbx.hiVect(), abx.loVect(), abx.hiVect(), 
           dbx.loVect(), dbx.hiVect(), bx.loVect(),  bx.hiVect(),
           (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
           (*u_go[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),      (*ro_g[lev])[mfi].dataPtr(),
           (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),   (*ep_g[lev])[mfi].dataPtr(),
           (*tau_u_g[lev])[mfi].dataPtr(),  (*d_e[lev])[mfi].dataPtr(),
           (*fluxX[lev])[mfi].dataPtr(),  (*fluxY[lev])[mfi].dataPtr(),  (*fluxZ[lev])[mfi].dataPtr(),
           (*mu_g[lev])[mfi].dataPtr(),     (*f_gds_u[lev])[mfi].dataPtr(), (*drag_u[lev])[mfi].dataPtr(),   
           (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),      (*mask)[mfi].dataPtr(),
           bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
           bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect(),
           &dt, &dx, &dy, &dz, &temp_num, &temp_denom);

       if (fluid_cost[lev]) {
	 const Box& tbx = mfi.tilebox(IntVect::TheZeroVector());
	 wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
	 (*fluid_cost[lev])[mfi].plus(wt, tbx);
       }
    }
    
    num_u = temp_num;
    denom_u = temp_denom;

    int eq_id=2;

    mfix_solve_linear_equation(eq_id,lev,(*u_gt[lev]),(*A_m[lev]),(*b_m[lev]));

    if (u_gt[lev]->contains_nan())
    {
        std::cout << "U_GT HAS NANS AFTER SOLVE" << std::endl;
        exit(0);
    }
}

void
mfix_level::mfix_solve_for_v(int lev, Real dt, Real& num_v, Real& denom_v)
{
    BL_PROFILE("mfix_level::solve_for_v()");
    Box domain(geom[lev].Domain());

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // Solve V-Momentum equation

    // Matrix and rhs vector
    BoxArray y_edge_ba = grids[lev];
    y_edge_ba.surroundingNodes(1);
    A_m[lev].reset(new MultiFab(y_edge_ba,dmap[lev],7,0));
    b_m[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,0));

    MultiFab::Copy(*v_gt[lev], *v_g[lev], 0, 0, 1, v_g[lev]->nGrow());

    // Initialize d_n to 0
    d_n[lev]->setVal(0.);

    auto mask = v_g[lev]->OverlapMask(geom[lev].periodicity());

    // We make temporaries because the omp parallel reduction won't take reference values.
    Real temp_num = num_v;
    Real temp_denom = denom_v;

#ifdef _OPENMP
#pragma omp parallel reduction(+:temp_num,temp_denom)
#endif
    for (MFIter mfi(*v_g[lev],true); mfi.isValid(); ++mfi)
    {

       Real wt = ParallelDescriptor::second();

       const Box& bx = mfi.tilebox();
       const Box& sbx = (*ep_g[lev])[mfi].box();
       Box abx((*A_m[lev])[mfi].box());

       Box ubx((*u_g[lev])[mfi].box());
       Box vbx((*v_g[lev])[mfi].box());
       Box wbx((*w_g[lev])[mfi].box());

       Box dbx((*drag_v[lev])[mfi].box());

       solve_v_g_star(sbx.loVect(), sbx.hiVect(),
           ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(),
           wbx.loVect(), wbx.hiVect(), abx.loVect(), abx.hiVect(),
           dbx.loVect(), dbx.hiVect(), bx.loVect(),  bx.hiVect(),
           (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
           (*v_go[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),      (*ro_g[lev])[mfi].dataPtr(),
           (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),   (*ep_g[lev])[mfi].dataPtr(),
           (*tau_v_g[lev])[mfi].dataPtr(),  (*d_n[lev])[mfi].dataPtr(),
           (*fluxX[lev])[mfi].dataPtr(),  (*fluxY[lev])[mfi].dataPtr(),  (*fluxZ[lev])[mfi].dataPtr(),
           (*mu_g[lev])[mfi].dataPtr(),     (*f_gds_v[lev])[mfi].dataPtr(), (*drag_v[lev])[mfi].dataPtr(),
           (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),      (*mask)[mfi].dataPtr(),
           bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
           bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect(),
           &dt, &dx, &dy, &dz, &temp_num, &temp_denom);

       if (fluid_cost[lev]) {
	 const Box& tbx = mfi.tilebox(IntVect::TheZeroVector());
	 wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
	 (*fluid_cost[lev])[mfi].plus(wt, tbx);
       }
    }

    num_v = temp_num;
    denom_v = temp_denom;

    int eq_id = 3;

    mfix_solve_linear_equation(eq_id,lev,(*v_gt[lev]),(*A_m[lev]),(*b_m[lev]));

    if (v_gt[lev]->contains_nan())
    {
        std::cout << "V_GT HAS NANS AFTER SOLVE" << std::endl;
        exit(0);
    }
}

void
mfix_level::mfix_solve_for_w(int lev, Real dt, Real& num_w, Real& denom_w)
{
    BL_PROFILE("mfix_level::solve_for_w()");
    Box domain(geom[lev].Domain());

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // Solve W-Momentum equation

    // Matrix and rhs vector
    BoxArray z_edge_ba = grids[lev];
    z_edge_ba.surroundingNodes(2);
    A_m[lev].reset(new MultiFab(z_edge_ba,dmap[lev],7,0));
    b_m[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,0));

    MultiFab::Copy(*w_gt[lev], *w_g[lev], 0, 0, 1, w_g[lev]->nGrow());

    // Initialize d_t to 0
    d_t[lev]->setVal(0.);
    auto mask = w_g[lev]->OverlapMask(geom[lev].periodicity());

    // We make temporaries because the omp parallel reduction won't take reference values.
    Real temp_num = num_w;
    Real temp_denom = denom_w;

#ifdef _OPENMP
#pragma omp parallel reduction(+:temp_num,temp_denom)
#endif
    for (MFIter mfi(*w_g[lev],true); mfi.isValid(); ++mfi)
    {

       Real wt = ParallelDescriptor::second();
      
       const Box& bx = mfi.tilebox();
       const Box& sbx = (*ep_g[lev])[mfi].box();
       Box abx((*A_m[lev])[mfi].box());

       Box ubx((*u_g[lev])[mfi].box());
       Box vbx((*v_g[lev])[mfi].box());
       Box wbx((*w_g[lev])[mfi].box());

       Box dbx((*drag_w[lev])[mfi].box());

       solve_w_g_star(sbx.loVect(), sbx.hiVect(),
           ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(),
           wbx.loVect(), wbx.hiVect(), abx.loVect(), abx.hiVect(),
           dbx.loVect(), dbx.hiVect(), bx.loVect(),  bx.hiVect(),
           (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
           (*w_go[lev])[mfi].dataPtr(),     (*p_g[lev])[mfi].dataPtr(),      (*ro_g[lev])[mfi].dataPtr(),
           (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),   (*ep_g[lev])[mfi].dataPtr(),
           (*tau_w_g[lev])[mfi].dataPtr(),  (*d_t[lev])[mfi].dataPtr(),
           (*fluxX[lev])[mfi].dataPtr(),  (*fluxY[lev])[mfi].dataPtr(),  (*fluxZ[lev])[mfi].dataPtr(),
           (*mu_g[lev])[mfi].dataPtr(),     (*f_gds_w[lev])[mfi].dataPtr(), (*drag_w[lev])[mfi].dataPtr(),
           (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),      (*mask)[mfi].dataPtr(),
           bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
           bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect(),
           &dt, &dx, &dy, &dz, &temp_num, &temp_denom);

       if (fluid_cost[lev]) {
	 const Box& tbx = mfi.tilebox(IntVect::TheZeroVector());
	 wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
	 (*fluid_cost[lev])[mfi].plus(wt, tbx);
       }
    }

    num_w = temp_num;
    denom_w = temp_denom;

    int eq_id = 4;

    mfix_solve_linear_equation(eq_id,lev,(*w_gt[lev]),(*A_m[lev]),(*b_m[lev]));

    if (w_gt[lev]->contains_nan())
    {
        std::cout << "W_GT HAS NANS AFTER SOLVE" << std::endl;
        exit(0);
    }
}

void
mfix_level::mfix_solve_for_pp(int lev, Real dt, Real& num_p, Real& denom_p)
{
    BL_PROFILE("mfix_level::solve_for_pp()");
    Box domain(geom[lev].Domain());

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // Matrix and rhs vector
    A_m[lev].reset(new MultiFab(grids[lev],dmap[lev],7,0));
    b_m[lev].reset(new MultiFab(grids[lev],dmap[lev],1,0));

    MultiFab b_mmax(b_m[lev]->boxArray(),dmap[lev],1,b_m[lev]->nGrow());
    b_mmax.setVal(0.);

    // We make temporaries because the omp parallel reduction won't take reference values.
    Real temp_num = num_p;
    Real temp_denom = denom_p;

#ifdef _OPENMP
#pragma omp parallel reduction(+:temp_num,temp_denom)
#endif
    for (MFIter mfi(*A_m[lev],true); mfi.isValid(); ++mfi)
    {

       Real wt = ParallelDescriptor::second();

       const Box& bx = mfi.tilebox();
       const Box& sbx = (*ep_g[lev])[mfi].box();

       Box abx((*A_m[lev])[mfi].box());
       Box ubx((*u_g[lev])[mfi].box());
       Box vbx((*v_g[lev])[mfi].box());
       Box wbx((*w_g[lev])[mfi].box());

       // Solve the pressure correction equation
       solve_pp_g(sbx.loVect(), sbx.hiVect(),
            ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
            abx.loVect(), abx.hiVect(), bx.loVect(),  bx.hiVect(),
            (*u_g[lev])[mfi].dataPtr(),      (*v_g[lev])[mfi].dataPtr(),      (*w_g[lev])[mfi].dataPtr(),
            (*p_g[lev])[mfi].dataPtr(),      (*ep_g[lev])[mfi].dataPtr(),
            (*rop_g[lev])[mfi].dataPtr(),    (*rop_go[lev])[mfi].dataPtr(),
            (     *ro_g[lev])[mfi].dataPtr(),
            (*ropX[lev])[mfi].dataPtr(),   (*ropY[lev])[mfi].dataPtr(),   (*ropZ[lev])[mfi].dataPtr(),
            (*d_e[lev])[mfi].dataPtr(),      (*d_n[lev])[mfi].dataPtr(),      (*d_t[lev])[mfi].dataPtr(),
            (*A_m[lev])[mfi].dataPtr(),      (*b_m[lev])[mfi].dataPtr(),           b_mmax[mfi].dataPtr(),
            bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
            bc_klo.dataPtr(), bc_khi.dataPtr(),
            &dt, &dx, &dy, &dz, domain.loVect(), domain.hiVect(), &temp_num, &temp_denom);

       if (fluid_cost[lev]) {
	 const Box& tbx = mfi.tilebox(IntVect::TheZeroVector());
	 wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
	 (*fluid_cost[lev])[mfi].plus(wt, tbx);
       }
    }
    num_p = temp_num;
    denom_p = temp_denom;

    pp_g[lev]->setVal(0.);

    int eq_id = 1;

    mfix_solve_linear_equation(eq_id,lev,(*pp_g[lev]),(*A_m[lev]),(*b_m[lev]));

    fill_mf_bc(lev,*pp_g[lev]);

    if (pp_g[lev]->contains_nan())
    {
        std::cout << "PP_G HAS NANS AFTER SOLVE" << std::endl;
        exit(0);
    }
}

void
mfix_level::mfix_correct_0(int lev)
{
  BL_PROFILE("mfix_level::correct0()");
    Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
    {
  const Box& bx = mfi.growntilebox();
  const Box& sbx = (*p_g[lev])[mfi].box();

  Box ubx((*u_g[lev])[mfi].box());
  Box vbx((*v_g[lev])[mfi].box());
  Box wbx((*w_g[lev])[mfi].box());

  correct_p_0( bx.loVect(),  bx.hiVect(),
         sbx.loVect(),  sbx.hiVect(),
         (*p_g[lev])[mfi].dataPtr(),
         (*pp_g[lev])[mfi].dataPtr());
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*u_g[lev],true); mfi.isValid(); ++mfi)
    {
  const Box& bx = mfi.tilebox();
  const Box& sbx = (*p_g[lev])[mfi].box();

  Box ubx((*u_g[lev])[mfi].box());
  Box vbx((*v_g[lev])[mfi].box());
  Box wbx((*w_g[lev])[mfi].box());

  correct_u_0( bx.loVect(),  bx.hiVect(),
         ubx.loVect(), ubx.hiVect(),
         sbx.loVect(),  sbx.hiVect(),
         (*pp_g[lev])[mfi].dataPtr(),
         (*u_g[lev])[mfi].dataPtr(),
         (*d_e[lev])[mfi].dataPtr());
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*v_g[lev],true); mfi.isValid(); ++mfi)
    {
  const Box& bx = mfi.tilebox();
  const Box& sbx = (*p_g[lev])[mfi].box();

  Box ubx((*u_g[lev])[mfi].box());
  Box vbx((*v_g[lev])[mfi].box());
  Box wbx((*w_g[lev])[mfi].box());

  correct_v_0( bx.loVect(),  bx.hiVect(),
         vbx.loVect(), vbx.hiVect(),
         sbx.loVect(),  sbx.hiVect(),
         (*pp_g[lev])[mfi].dataPtr(),
         (*v_g[lev])[mfi].dataPtr(),
         (*d_n[lev])[mfi].dataPtr());
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*w_g[lev],true); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.tilebox();
       const Box& sbx = (*p_g[lev])[mfi].box();

       Box ubx((*u_g[lev])[mfi].box());
       Box vbx((*v_g[lev])[mfi].box());
       Box wbx((*w_g[lev])[mfi].box());

       correct_w_0( bx.loVect(),  bx.hiVect(),
              wbx.loVect(), wbx.hiVect(),
              sbx.loVect(),  sbx.hiVect(),
              (*pp_g[lev])[mfi].dataPtr(),
              (*w_g[lev])[mfi].dataPtr(),
              (*d_t[lev])[mfi].dataPtr());
    }

    fill_mf_bc(lev,*p_g[lev]);

    u_g[lev]->FillBoundary(geom[lev].periodicity());
    v_g[lev]->FillBoundary(geom[lev].periodicity());
    w_g[lev]->FillBoundary(geom[lev].periodicity());
}
