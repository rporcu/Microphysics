#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

std::string mfix_level::particle_init_type = "AsciiFile";
std::string mfix_level::load_balance_type = "FixedSize";

mfix_level::~mfix_level ()
{};

mfix_level::mfix_level ()
{
    // Geometry on all levels has just been defined in the AmrCore constructor

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;
#if 0
    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= maxLevel(); ++lev) {
  nsubsteps[lev] = MaxRefRatio(lev-1);
    }
#endif

    // Particle Container
    pc = std::unique_ptr<MFIXParticleContainer> (new MFIXParticleContainer(this));

    A_m.resize(nlevs_max);
    b_m.resize(nlevs_max);

    ep_g.resize(nlevs_max);
    ep_go.resize(nlevs_max);

    p_g.resize(nlevs_max);
    p_go.resize(nlevs_max);

    pp_g.resize(nlevs_max);

    ro_g.resize(nlevs_max);
    ro_go.resize(nlevs_max);

    rop_g.resize(nlevs_max);
    rop_go.resize(nlevs_max);

    u_g.resize(nlevs_max);
    u_go.resize(nlevs_max);
    u_gt.resize(nlevs_max);

    v_g.resize(nlevs_max);
    v_go.resize(nlevs_max);
    v_gt.resize(nlevs_max);

    w_g.resize(nlevs_max);
    w_go.resize(nlevs_max);
    w_gt.resize(nlevs_max);

    d_e.resize(nlevs_max);
    d_n.resize(nlevs_max);
    d_t.resize(nlevs_max);

    mu_g.resize(nlevs_max);
    lambda_g.resize(nlevs_max);
    trD_g.resize(nlevs_max);

    fluxX.resize(nlevs_max);
    fluxY.resize(nlevs_max);
    fluxZ.resize(nlevs_max);

    ropX.resize(nlevs_max);
    ropY.resize(nlevs_max);
    ropZ.resize(nlevs_max);

    tau_u_g.resize(nlevs_max);
    tau_v_g.resize(nlevs_max);
    tau_w_g.resize(nlevs_max);

    f_gds_u.resize(nlevs_max);
    f_gds_v.resize(nlevs_max);
    f_gds_w.resize(nlevs_max);

    drag_u.resize(nlevs_max);
    drag_v.resize(nlevs_max);
    drag_w.resize(nlevs_max);
}

void mfix_level::mfix_calc_coeffs(int lev, int calc_flag)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*ep_g[lev],true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& sbx = (*ep_g[lev])[mfi].box();

        calc_coeff(sbx.loVect(), sbx.hiVect(), bx.loVect(),  bx.hiVect(), &calc_flag,
                   (*ro_g[lev])[mfi].dataPtr(), (*p_g[lev])[mfi].dataPtr(),
                   (*ep_g[lev])[mfi].dataPtr(), (*rop_g[lev])[mfi].dataPtr());
    }

    fill_mf_bc(lev,*ro_g[lev]);
    fill_mf_bc(lev,*rop_g[lev]);
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

void
mfix_level::mfix_physical_prop(int lev, int calc_flag)
{
  BL_PROFILE("mfix_level::mfix_physical_prop()");
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*p_g[lev],true); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.tilebox();
       const Box& sbx = (*p_g[lev])[mfi].box();

       physical_prop(sbx.loVect(), sbx.hiVect(), bx.loVect(), bx.hiVect(),&calc_flag,
               (*ro_g[lev])[mfi].dataPtr(), (*p_g[lev])[mfi].dataPtr(),
               (*ep_g[lev])[mfi].dataPtr(), (*rop_g[lev])[mfi].dataPtr());
    }
    fill_mf_bc(lev,*ro_g[lev]);
    fill_mf_bc(lev,*rop_g[lev]);
}

void
mfix_level::usr3(int lev)
{
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // We deliberately don't tile this loop since we will be looping
    //    over bc's on faces and it makes more sense to do this one grid at a time
    for (MFIter mfi(*p_g[lev]); mfi.isValid(); ++mfi)
    {
  const Box& sbx = (*p_g[lev])[mfi].box();
  Box ubx((*u_g[lev])[mfi].box());
  Box vbx((*v_g[lev])[mfi].box());
  Box wbx((*w_g[lev])[mfi].box());

  mfix_usr3((*u_g[lev])[mfi].dataPtr(), ubx.loVect(), ubx.hiVect(),
      (*v_g[lev])[mfi].dataPtr(), vbx.loVect(), vbx.hiVect(),
      (*w_g[lev])[mfi].dataPtr(), wbx.loVect(), wbx.hiVect(),
      (*p_g[lev])[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
      &dx, &dy, &dz);
    }
}

void
mfix_level::mfix_solve_linear_equation(int eq_id,int lev,MultiFab& sol, MultiFab& matrix, MultiFab& rhs)
{
    int sweep_type, precond_type, max_it;
    Real tol;

    BL_PROFILE("mfix_level::mfix_solve_linear_equation()");
    get_solver_params (&eq_id,&sweep_type,&precond_type,&max_it,&tol);

    solve_bicgstab(sol, rhs, matrix, sweep_type, precond_type, max_it, tol, lev);
}

void
mfix_level::mfix_set_bc_type(int lev)
{
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);
    Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
    Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
    Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);
    Box domain(geom[lev].Domain());

    set_bc_type(bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                bc_klo.dataPtr(), bc_khi.dataPtr(),
                domain.loVect(),domain.hiVect(),
                &dx, &dy, &dz, &xlen, &ylen, &zlen, &nghost_bc);
}

void
mfix_level::fill_mf_bc(int lev, MultiFab& mf)
{
    Box domain(geom[lev].Domain());

    if (!mf.boxArray().ixType().cellCentered())
  amrex::Error("fill_mf_bc only used for cell-centered arrays!");

    // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
    mf.FillBoundary(geom[lev].periodicity());

    // Fill all cell-centered arrays with first-order extrapolation at domain boundaries
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
  const Box& sbx = mf[mfi].box();
  fill_bc0(mf[mfi].dataPtr(),sbx.loVect(),sbx.hiVect(),
     bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
     bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect());
    }
}

void mfix_level::mfix_calc_volume_fraction(int lev, Real& sum_vol)
{
    BL_PROFILE("mfix_level::mfix_calc_volume_fraction()");

    Box domain(geom[lev].Domain());

    // This re-calculates the volume fraction within the domain
    // but does not change the values outside the domain

    // This call simply deposits the particle volume onto the grid in a PIC-like manner
    pc->CalcVolumeFraction(*ep_g[lev],bc_ilo,bc_ihi,bc_jlo,bc_jhi,bc_klo,bc_khi);

    // Now define rop_g = ro_g * ep_g
    rop_g[lev]->copy((*ro_g[lev]),0,0,1,rop_g[lev]->nGrow(),rop_g[lev]->nGrow());
    MultiFab::Multiply((*rop_g[lev]), (*ep_g[lev]), 0, 0, 1, rop_g[lev]->nGrow());

    // This sets the values outside walls or periodic boundaries
    fill_mf_bc(lev,*ep_g[lev]);
    fill_mf_bc(lev,*rop_g[lev]);

    // Sum up all the values of ep_g[lev] -- this value should never change!
    sum_vol = ep_g[lev]->sum();
}

void mfix_level::mfix_calc_drag_fluid(int lev)
{
  BL_PROFILE("mfix_level::mfix_calc_drag_fluid()");
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    if (OnSameGrids) 
    {
       // ************************************************************
       // First create the beta of individual particles 
       // ************************************************************
#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
       {
           const Box& sbx = (*ep_g[lev])[pti].box();
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();

           Box ubx((*u_g[lev])[pti].box());
           Box vbx((*v_g[lev])[pti].box());
           Box wbx((*w_g[lev])[pti].box());

           calc_particle_beta(
               sbx.loVect(), sbx.hiVect(),
               ubx.loVect(), ubx.hiVect(),
               vbx.loVect(), vbx.hiVect(),
               wbx.loVect(), wbx.hiVect(), &np,
               (*ep_g[lev])[pti].dataPtr(), (*ro_g[lev])[pti].dataPtr(),
               (*u_g[lev])[pti].dataPtr(),  (*v_g[lev])[pti].dataPtr(),
               (*w_g[lev])[pti].dataPtr(),  (*mu_g[lev])[pti].dataPtr(),
               particles.data(), &dx, &dy, &dz);
       }

       // ******************************************************************************
       // Now use the beta of individual particles to create the drag terms on the fluid 
       // ******************************************************************************

       f_gds_u[lev]->setVal(0.0L);
       f_gds_v[lev]->setVal(0.0L);
       f_gds_w[lev]->setVal(0.0L);
       drag_u[lev]->setVal(0.0L);
       drag_v[lev]->setVal(0.0L);
       drag_w[lev]->setVal(0.0L);

       pc -> CalcDragOnFluid(*f_gds_u[lev],*f_gds_v[lev],*f_gds_w[lev],
                             *drag_u[lev],*drag_v[lev],*drag_w[lev],
                              bc_ilo,bc_ihi,bc_jlo,bc_jhi,bc_klo,bc_khi);
    }
    else 
    {

       BoxArray            pba = pc->ParticleBoxArray(lev);
       DistributionMapping pdm = pc->ParticleDistributionMap(lev);

       // Temporary arrays
       int ng = ep_g[lev]->nGrow();
       std::unique_ptr<MultiFab> ep_g_pba(new MultiFab(pba,pdm,ep_g[lev]->nComp(),ng));
       ep_g_pba->copy(*ep_g[lev],0,0,1,ng,ng,geom[lev].periodicity());

       ng = ro_g[lev]->nGrow();
       std::unique_ptr<MultiFab> ro_g_pba(new MultiFab(pba,pdm,ro_g[lev]->nComp(),ro_g[lev]->nGrow()));
       ro_g_pba->copy(*ro_g[lev],0,0,1,ng,ng,geom[lev].periodicity());

       ng = mu_g[lev]->nGrow();
       std::unique_ptr<MultiFab> mu_g_pba(new MultiFab(pba,pdm,mu_g[lev]->nComp(),mu_g[lev]->nGrow()));
       mu_g_pba->copy(*mu_g[lev],0,0,1,ng,ng,geom[lev].periodicity());

       BoxArray x_face_ba = pba;
       x_face_ba.surroundingNodes(0);
       std::unique_ptr<MultiFab> u_g_pba(new MultiFab(x_face_ba,pdm,u_g[lev]->nComp(),u_g[lev]->nGrow()));
       u_g_pba->copy(*u_g[lev],0,0,1,ng,ng,geom[lev].periodicity());

       BoxArray y_face_ba = pba;
       y_face_ba.surroundingNodes(1);
       std::unique_ptr<MultiFab> v_g_pba(new MultiFab(y_face_ba,pdm,v_g[lev]->nComp(),v_g[lev]->nGrow()));
       v_g_pba->copy(*v_g[lev]);

       BoxArray z_face_ba = pba;
       z_face_ba.surroundingNodes(2);
       std::unique_ptr<MultiFab> w_g_pba(new MultiFab(z_face_ba,pdm,w_g[lev]->nComp(),w_g[lev]->nGrow()));
       w_g_pba->copy(*w_g[lev]);

       // ************************************************************
       // First create the beta of individual particles 
       // ************************************************************

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
       {
           const Box& sbx = (*ep_g_pba)[pti].box();
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();

           Box ubx((*u_g_pba)[pti].box());
           Box vbx((*v_g_pba)[pti].box());
           Box wbx((*w_g_pba)[pti].box());

           calc_particle_beta(
               sbx.loVect(), sbx.hiVect(),
               ubx.loVect(), ubx.hiVect(),
               vbx.loVect(), vbx.hiVect(),
               wbx.loVect(), wbx.hiVect(), &np,
               (*ep_g_pba)[pti].dataPtr(), (*ro_g_pba)[pti].dataPtr(),
               (*u_g_pba)[pti].dataPtr(),  (*v_g_pba)[pti].dataPtr(),
               (*w_g_pba)[pti].dataPtr(),  (*mu_g_pba)[pti].dataPtr(),
               particles.data(), &dx, &dy, &dz);
       }

       // ******************************************************************************
       // Now use the beta of individual particles to create the drag terms on the fluid 
       // ******************************************************************************

       std::unique_ptr<MultiFab> f_gds_v_pba(new MultiFab(y_face_ba,pdm,f_gds_v[lev]->nComp(),f_gds_v[lev]->nGrow()));
       std::unique_ptr<MultiFab> drag_v_pba(new MultiFab(y_face_ba,pdm,drag_v[lev]->nComp(),drag_v[lev]->nGrow()));

       std::unique_ptr<MultiFab> f_gds_u_pba(new MultiFab(x_face_ba,pdm,f_gds_u[lev]->nComp(),f_gds_u[lev]->nGrow()));
       std::unique_ptr<MultiFab> drag_u_pba(new MultiFab(x_face_ba,pdm,drag_u[lev]->nComp(),drag_u[lev]->nGrow()));

       std::unique_ptr<MultiFab> f_gds_w_pba(new MultiFab(z_face_ba,pdm,f_gds_w[lev]->nComp(),f_gds_w[lev]->nGrow()));
       std::unique_ptr<MultiFab> drag_w_pba(new MultiFab(z_face_ba,pdm,drag_w[lev]->nComp(),drag_w[lev]->nGrow()));

       f_gds_u_pba->setVal(0.0L);
       f_gds_v_pba->setVal(0.0L);
       f_gds_w_pba->setVal(0.0L);
       drag_u_pba->setVal(0.0L);
       drag_v_pba->setVal(0.0L);
       drag_w_pba->setVal(0.0L);

       pc -> CalcDragOnFluid(*f_gds_u_pba,*f_gds_v_pba,*f_gds_w_pba,
                             *drag_u_pba,*drag_v_pba,*drag_w_pba,
                              bc_ilo,bc_ihi,bc_jlo,bc_jhi,bc_klo,bc_khi);

       // Copy back from the dual grids.
       f_gds_u[lev] ->copy(*f_gds_u_pba);
       f_gds_v[lev] ->copy(*f_gds_v_pba);
       f_gds_w[lev] ->copy(*f_gds_w_pba);

       drag_u[lev] ->copy(*drag_u_pba);
       drag_v[lev] ->copy(*drag_v_pba);
       drag_w[lev] ->copy(*drag_w_pba);

    } // if not OnSameGrids

    // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
    f_gds_u[lev]->FillBoundary(geom[lev].periodicity());
    f_gds_v[lev]->FillBoundary(geom[lev].periodicity());
    f_gds_w[lev]->FillBoundary(geom[lev].periodicity());

    drag_u[lev]->FillBoundary(geom[lev].periodicity());
    drag_v[lev]->FillBoundary(geom[lev].periodicity());
    drag_w[lev]->FillBoundary(geom[lev].periodicity());
}

void
mfix_level::mfix_calc_drag_particle(int lev)
{
    BL_PROFILE("mfix_level::mfix_calc_drag_particle()");

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
    Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
    Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    if (OnSameGrids) 
    {

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
       {
           const Box& sbx = (*p_g[lev])[pti].box();
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();

           Box ubx((*u_g[lev])[pti].box());
           Box vbx((*v_g[lev])[pti].box());
           Box wbx((*w_g[lev])[pti].box());
   
           calc_drag_particle(
               sbx.loVect(), sbx.hiVect(),
               ubx.loVect(), ubx.hiVect(),
               vbx.loVect(), vbx.hiVect(),
               wbx.loVect(), wbx.hiVect(), &np,
               (*p_g[lev])[pti].dataPtr(), (*u_g[lev])[pti].dataPtr(),
               (*v_g[lev])[pti].dataPtr(), (*w_g[lev])[pti].dataPtr(),
               particles.data(), &dx, &dy, &dz, &xlen, &ylen, &zlen);
       }
    }
    else 
    {

       BoxArray            pba = pc->ParticleBoxArray(lev);
       DistributionMapping pdm = pc->ParticleDistributionMap(lev);

       // Temporary arrays
       int ng = p_g[lev]->nGrow();
       std::unique_ptr<MultiFab> p_g_pba(new MultiFab(pba,pdm,p_g[lev]->nComp(),ng));
       p_g_pba->copy(*p_g[lev],0,0,1,ng,ng,geom[lev].periodicity());

       BoxArray x_face_ba = pba;
       x_face_ba.surroundingNodes(0);
       ng = u_g[lev]->nGrow();
       std::unique_ptr<MultiFab> u_g_pba(new MultiFab(x_face_ba,pdm,u_g[lev]->nComp(),ng));
       u_g_pba->copy(*u_g[lev],0,0,1,ng,ng,geom[lev].periodicity());

       BoxArray y_face_ba = pba;
       y_face_ba.surroundingNodes(1);
       ng = v_g[lev]->nGrow();
       std::unique_ptr<MultiFab> v_g_pba(new MultiFab(y_face_ba,pdm,v_g[lev]->nComp(),ng));
       v_g_pba->copy(*v_g[lev],0,0,1,ng,ng,geom[lev].periodicity());

       BoxArray z_face_ba = pba;
       z_face_ba.surroundingNodes(2);
       ng = w_g[lev]->nGrow();
       std::unique_ptr<MultiFab> w_g_pba(new MultiFab(z_face_ba,pdm,w_g[lev]->nComp(),ng));
       w_g_pba->copy(*w_g[lev],0,0,1,ng,ng,geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
       {
           const Box& sbx = (*p_g_pba)[pti].box();
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();

           Box ubx((*u_g_pba)[pti].box());
           Box vbx((*v_g_pba)[pti].box());
           Box wbx((*w_g_pba)[pti].box());
   
           calc_drag_particle(
               sbx.loVect(), sbx.hiVect(),
               ubx.loVect(), ubx.hiVect(),
               vbx.loVect(), vbx.hiVect(),
               wbx.loVect(), wbx.hiVect(), &np,
               (*p_g_pba)[pti].dataPtr(), (*u_g_pba)[pti].dataPtr(),
               (*v_g_pba)[pti].dataPtr(), (*w_g_pba)[pti].dataPtr(),
               particles.data(), &dx, &dy, &dz, &xlen, &ylen, &zlen);
       }
    }
}

void
mfix_level::mfix_set_bc1(int lev)
{
  BL_PROFILE("mfix_level::mfix_set_bc1()");

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
    {
      Box domain(geom[lev].Domain());
      const Box& sbx = (*ep_g[lev])[mfi].box();
      Box ubx((*u_g[lev])[mfi].box());
      Box vbx((*v_g[lev])[mfi].box());
      Box wbx((*w_g[lev])[mfi].box());

      set_bc1(sbx.loVect(), sbx.hiVect(),
              ubx.loVect(), ubx.hiVect(), vbx.loVect(), vbx.hiVect(), wbx.loVect(), wbx.hiVect(),
              (*u_g[lev])[mfi].dataPtr(), (*v_g[lev])[mfi].dataPtr(), (*w_g[lev])[mfi].dataPtr(),
              (*p_g[lev])[mfi].dataPtr(),     (*ep_g[lev])[mfi].dataPtr(),
              (*ro_g[lev])[mfi].dataPtr(), (*rop_g[lev])[mfi].dataPtr(),
              (*mu_g[lev])[mfi].dataPtr(), (*lambda_g[lev])[mfi].dataPtr(),
              bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
              bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect());
    }
}
