#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

std::string mfix::particle_init_type   = "AsciiFile";
std::string mfix::load_balance_type    = "FixedSize";
std::string mfix::knapsack_weight_type = "RunTimeCosts";

// Define unit vectors for easily convert indeces
amrex::IntVect mfix::e_x(1,0,0);
amrex::IntVect mfix::e_y(0,1,0);
amrex::IntVect mfix::e_z(0,0,1);

int mfix::m_eb_basic_grow_cells = 2;
int mfix::m_eb_volume_grow_cells = 2;
int mfix::m_eb_full_grow_cells = 2;
EBSupport mfix::m_eb_support_level = EBSupport::full;


mfix::~mfix ()
{
};


mfix::mfix ()
{
    // Geometry on all levels has just been defined in the AmrCore constructor

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

#if 0
    int nlevs_max = maxLevel() + 1;
    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= maxLevel(); ++lev)
        nsubsteps[lev] = MaxRefRatio(lev-1);
#endif
}

void
mfix::ResizeArrays ()
{
    int nlevs_max = maxLevel() + 1;

    // Particle Container
    pc = std::unique_ptr<MFIXParticleContainer> (new MFIXParticleContainer(this));

    // HACK: temporary flag used to turn on legacy mode
    //   (used in evlove particles)
    pc -> legacy__eb_collisions = legacy__eb_collisions;

    A_m.resize(nlevs_max);
    b_m.resize(nlevs_max);

    ep_g.resize(nlevs_max);
    ep_go.resize(nlevs_max);

    p_g.resize(nlevs_max);
    p_go.resize(nlevs_max);

    p0_g.resize(nlevs_max);
    pp_g.resize(nlevs_max);

    ro_g.resize(nlevs_max);
    ro_go.resize(nlevs_max);

    rop_g.resize(nlevs_max);
    rop_go.resize(nlevs_max);

    phi.resize(nlevs_max);
    diveu.resize(nlevs_max);

    u_g.resize(nlevs_max);
    u_go.resize(nlevs_max);
    u_gt.resize(nlevs_max);

    v_g.resize(nlevs_max);
    v_go.resize(nlevs_max);
    v_gt.resize(nlevs_max);

    w_g.resize(nlevs_max);
    w_go.resize(nlevs_max);
    w_gt.resize(nlevs_max);

    fp_x.resize(nlevs_max);
    fp_y.resize(nlevs_max);
    fp_z.resize(nlevs_max);

    bcoeff.resize(nlevs_max);
    for(int i = 0; i < nlevs_max; ++i)
        bcoeff[i].resize(3);

    d_e.resize(nlevs_max);
    d_n.resize(nlevs_max);
    d_t.resize(nlevs_max);

    mu_g.resize(nlevs_max);
    lambda_g.resize(nlevs_max);
    trD_g.resize(nlevs_max);
    vort.resize(nlevs_max);

    // MultiFab storing level-set data
    ls.resize(nlevs_max);

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

    slopes_u.resize(nlevs_max);
    slopes_v.resize(nlevs_max);
    slopes_w.resize(nlevs_max);

    uacc.resize(nlevs_max);
    vacc.resize(nlevs_max);
    wacc.resize(nlevs_max);

    if(solve_dem)
       particle_cost.resize(nlevs_max);
    if(solve_fluid)
       fluid_cost.resize(nlevs_max);

    // EB factory
    ebfactory.resize(nlevs_max);
    particle_ebfactory.resize(nlevs_max);
}

void
mfix::usr3(int lev)
{
    if(solve_fluid)
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

          mfix_usr3( (*u_g[lev])[mfi].dataPtr(), ubx.loVect(), ubx.hiVect(),
                     (*v_g[lev])[mfi].dataPtr(), vbx.loVect(), vbx.hiVect(),
                     (*w_g[lev])[mfi].dataPtr(), wbx.loVect(), wbx.hiVect(),
                     (*p_g[lev])[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
                     & dx, & dy, & dz
                   );
       }
    }
}

void
mfix::mfix_set_bc_type(int lev)
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
                &dx, &dy, &dz, &xlen, &ylen, &zlen, &nghost);
}

void mfix::fill_mf_bc(int lev, MultiFab & mf) {
    Box domain(geom[lev].Domain());

    if(!mf.boxArray().ixType().cellCentered())
        amrex::Error("fill_mf_bc only used for cell-centered arrays!");

    // Impose periodic bc's at domain boundaries and fine-fine copies in the
    // interior It is essential that we do this before the call to fill_bc0
    // below since fill_bc0 can extrapolate out to fill ghost cells outside the
    // domain after we have filled ghost cells inside the domain, but doing
    // this call after fill_bc0 can't fill ghost cells from ghost cells.
    mf.FillBoundary(geom[lev].periodicity());

    // Fill all cell-centered arrays with first-order extrapolation at domain
    // boundaries
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
        const Box& sbx = mf[mfi].box();
        fill_bc0(mf[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
                 bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                 bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                 bc_klo.dataPtr(), bc_khi.dataPtr(),
                 domain.loVect(),  domain.hiVect(),
                 &nghost                                        );
    }

    // Impose periodic bc's at domain boundaries and fine-fine copies in the
    // interior It's not 100% clear whether we need this call or not.  Worth
    // testing.
    mf.FillBoundary(geom[lev].periodicity());
}

//
// Set the BCs for velocity only
//
void
mfix::mfix_set_velocity_bcs (int lev)
{
  BL_PROFILE("mfix::mfix_set_velocity_bcs()");

  u_g[lev] -> FillBoundary (geom[lev].periodicity());
  v_g[lev] -> FillBoundary (geom[lev].periodicity());
  w_g[lev] -> FillBoundary (geom[lev].periodicity());

  Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(*p_g[lev], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = (*p_g[lev])[mfi].box();
      set_mac_velocity_bcs ( bx.loVect(), bx.hiVect(),
			     BL_TO_FORTRAN_ANYD((*u_g[lev])[mfi]),
			     BL_TO_FORTRAN_ANYD((*v_g[lev])[mfi]),
			     BL_TO_FORTRAN_ANYD((*w_g[lev])[mfi]),
			     bc_ilo.dataPtr(), bc_ihi.dataPtr(),
			     bc_jlo.dataPtr(), bc_jhi.dataPtr(),
			     bc_klo.dataPtr(), bc_khi.dataPtr(),
			     domain.loVect(), domain.hiVect(),
			     &nghost );
    }
}

void mfix::mfix_calc_volume_fraction(int lev, Real & sum_vol) {
    BL_PROFILE("mfix::mfix_calc_volume_fraction()");

    if(solve_dem) {
       // This re-calculates the volume fraction within the domain
       // but does not change the values outside the domain

       // This call simply deposits the particle volume onto the grid in a
       // PIC-like manner
       pc->CalcVolumeFraction( * ep_g[lev],
                              bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo, bc_khi,
                              nghost                                         );
    } else {
       ep_g[lev]->setVal(1.);
    }

    // Now define rop_g = ro_g * ep_g
    MultiFab::Copy(*rop_g[lev], *ro_g[lev], 0, 0, 1, ro_g[lev]->nGrow());
    MultiFab::Multiply((*rop_g[lev]), (*ep_g[lev]), 0, 0, 1, rop_g[lev]->nGrow());

    // This sets the values outside walls or periodic boundaries
    fill_mf_bc(lev,*ep_g[lev]);
    fill_mf_bc(lev,*rop_g[lev]);

    // Sum up all the values of ep_g[lev] -- this value should never change!
    sum_vol = ep_g[lev]->sum();
}

void mfix::mfix_calc_drag_fluid(int lev)
{
  BL_PROFILE("mfix::mfix_calc_drag_fluid()");
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    // Make sure to fill tangential velocities outside the domain for use in
    // interpolation onto particle positions in the calc_particle_beta calculation
    mfix_set_velocity_bcs(lev);

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    if (OnSameGrids) {

       // ************************************************************
       // First create the beta of individual particles
       // ************************************************************
#ifdef _OPENMP
#pragma omp parallel
#endif
       for(MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {
           const Box& sbx = (*ep_g[lev])[pti].box();
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();

           Box ubx((*u_g[lev])[pti].box());
           Box vbx((*v_g[lev])[pti].box());
           Box wbx((*w_g[lev])[pti].box());

           calc_particle_beta( sbx.loVect(), sbx.hiVect(),
                               ubx.loVect(), ubx.hiVect(),
                               vbx.loVect(), vbx.hiVect(),
                               wbx.loVect(), wbx.hiVect(), &np,
                               (*ep_g[lev])[pti].dataPtr(), (*ro_g[lev])[pti].dataPtr(),
                               (*u_g[lev])[pti].dataPtr(),  (*v_g[lev])[pti].dataPtr(),
                               (*w_g[lev])[pti].dataPtr(),  (*mu_g[lev])[pti].dataPtr(),
                               particles.data(), &dx, &dy, &dz
                             );
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

       pc -> CalcDragOnFluid(* f_gds_u[lev], * f_gds_v[lev], * f_gds_w[lev],
                             * drag_u[lev],  * drag_v[lev],  * drag_w[lev],
                             bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo,bc_khi,
                             nghost                                         );
    } else {

       BoxArray            pba = pc->ParticleBoxArray(lev);
       DistributionMapping pdm = pc->ParticleDistributionMap(lev);

       // Temporary arrays
       int ng = ep_g[lev]->nGrow();
       std::unique_ptr<MultiFab> ep_g_pba(new MultiFab(pba,pdm,ep_g[lev]->nComp(),ng));
       ep_g_pba->copy(*ep_g[lev],0,0,1,ng,ng);
       ep_g_pba->FillBoundary(geom[lev].periodicity());

       ng = ro_g[lev]->nGrow();
       std::unique_ptr<MultiFab> ro_g_pba(new MultiFab(pba,pdm,ro_g[lev]->nComp(),ro_g[lev]->nGrow()));
       ro_g_pba->copy(*ro_g[lev],0,0,1,ng,ng);
       ro_g_pba->FillBoundary(geom[lev].periodicity());

       ng = mu_g[lev]->nGrow();
       std::unique_ptr<MultiFab> mu_g_pba(new MultiFab(pba,pdm,mu_g[lev]->nComp(),mu_g[lev]->nGrow()));
       mu_g_pba->copy(*mu_g[lev],0,0,1,ng,ng);
       mu_g_pba->FillBoundary(geom[lev].periodicity());

       BoxArray x_face_ba = pba;
       x_face_ba.surroundingNodes(0);
       std::unique_ptr<MultiFab> u_g_pba(new MultiFab(x_face_ba,pdm,u_g[lev]->nComp(),u_g[lev]->nGrow()));
       u_g_pba->copy(*u_g[lev],0,0,1,ng,ng);
       u_g_pba->FillBoundary(geom[lev].periodicity());

       BoxArray y_face_ba = pba;
       y_face_ba.surroundingNodes(1);
       std::unique_ptr<MultiFab> v_g_pba(new MultiFab(y_face_ba,pdm,v_g[lev]->nComp(),v_g[lev]->nGrow()));
       v_g_pba->copy(*v_g[lev],0,0,1,ng,ng);
       v_g_pba->FillBoundary(geom[lev].periodicity());

       BoxArray z_face_ba = pba;
       z_face_ba.surroundingNodes(2);
       std::unique_ptr<MultiFab> w_g_pba(new MultiFab(z_face_ba,pdm,w_g[lev]->nComp(),w_g[lev]->nGrow()));
       w_g_pba->copy(*w_g[lev],0,0,1,ng,ng);
       w_g_pba->FillBoundary(geom[lev].periodicity());

       // ************************************************************
       // First create the beta of individual particles
       // ************************************************************

#ifdef _OPENMP
#pragma omp parallel
#endif
       for(MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {
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

       pc -> CalcDragOnFluid(* f_gds_u_pba, * f_gds_v_pba, * f_gds_w_pba,
                             * drag_u_pba,  * drag_v_pba,  * drag_w_pba,
                             bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo, bc_khi,
                             nghost                                         );

       // Copy back from the dual grids.
       f_gds_u[lev] ->copy(*f_gds_u_pba);
       f_gds_v[lev] ->copy(*f_gds_v_pba);
       f_gds_w[lev] ->copy(*f_gds_w_pba);

       drag_u[lev] ->copy(*drag_u_pba);
       drag_v[lev] ->copy(*drag_v_pba);
       drag_w[lev] ->copy(*drag_w_pba);

    } // if not OnSameGrids

    // The projection method uses drag to update u, not (cell_vol * u), so we must divide by vol here
    //     and we will divide by density in the update.
    if (use_proj_method)
    {
        Real ovol = 1./(dx*dy*dz);
         drag_u[lev]->mult(ovol);
         drag_v[lev]->mult(ovol);
         drag_w[lev]->mult(ovol);
        f_gds_u[lev]->mult(ovol);
        f_gds_v[lev]->mult(ovol);
        f_gds_w[lev]->mult(ovol);
    }

    // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
    f_gds_u[lev]->FillBoundary(geom[lev].periodicity());
    f_gds_v[lev]->FillBoundary(geom[lev].periodicity());
    f_gds_w[lev]->FillBoundary(geom[lev].periodicity());

    drag_u[lev]->FillBoundary(geom[lev].periodicity());
    drag_v[lev]->FillBoundary(geom[lev].periodicity());
    drag_w[lev]->FillBoundary(geom[lev].periodicity());
}

void
mfix::mfix_calc_drag_particle(int lev)
{
    BL_PROFILE("mfix::mfix_calc_drag_particle()");

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    Box domain(geom[lev].Domain());
    if (OnSameGrids)
    {
       // Temporary arrays
       std::unique_ptr<MultiFab> gpx(new MultiFab(u_g[lev]->boxArray(),dmap[lev],1,1));
       std::unique_ptr<MultiFab> gpy(new MultiFab(v_g[lev]->boxArray(),dmap[lev],1,1));
       std::unique_ptr<MultiFab> gpz(new MultiFab(w_g[lev]->boxArray(),dmap[lev],1,1));

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIter mfi(*p_g[lev], true); mfi.isValid(); ++mfi)
       {
           const Box& sbx = (*p_g[lev])[mfi].box();

           construct_gradp( sbx.loVect(),   sbx.hiVect(),
                            (*p_g[lev])[mfi].dataPtr(), (*p0_g[lev])[mfi].dataPtr(),
                            BL_TO_FORTRAN_ANYD((*gpx)[mfi]),
                            BL_TO_FORTRAN_ANYD((*gpy)[mfi]),
                            BL_TO_FORTRAN_ANYD((*gpz)[mfi]),
                            &dx, &dy, &dz,
                            bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                            bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect(),
                            &nghost
                          );
       }

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIter mfi(*p_g[lev], true); mfi.isValid(); ++mfi)
       {
           const Box& sbx = (*p_g[lev])[mfi].box();
           set_gradp_bcs ( sbx.loVect(), sbx.hiVect(),
                           BL_TO_FORTRAN_ANYD((*gpx)[mfi]),
                           BL_TO_FORTRAN_ANYD((*gpy)[mfi]),
                           BL_TO_FORTRAN_ANYD((*gpz)[mfi]),
                           bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                           bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                           bc_klo.dataPtr(), bc_khi.dataPtr(),
                           domain.loVect(), domain.hiVect(),
                           &nghost );
       }

       gpx->FillBoundary(geom[lev].periodicity());
       gpy->FillBoundary(geom[lev].periodicity());
       gpz->FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
       {
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();

           calc_drag_particle( BL_TO_FORTRAN_ANYD((*gpx)[pti]),
                               BL_TO_FORTRAN_ANYD((*gpy)[pti]),
                               BL_TO_FORTRAN_ANYD((*gpz)[pti]),
                               BL_TO_FORTRAN_ANYD((*u_g[lev])[pti]),
                               BL_TO_FORTRAN_ANYD((*v_g[lev])[pti]),
                               BL_TO_FORTRAN_ANYD((*w_g[lev])[pti]),
                               &np, particles.data(), &dx, &dy, &dz
                             );
       }
    }
    else
    {
       BoxArray            pba = pc->ParticleBoxArray(lev);
       DistributionMapping pdm = pc->ParticleDistributionMap(lev);

       // Temporary arrays
       int ng = p_g[lev]->nGrow();
       std::unique_ptr<MultiFab> p_g_pba(new MultiFab(pba,pdm,p_g[lev]->nComp(),ng));
       p_g_pba->copy(*p_g[lev],0,0,1,ng,ng);
       p_g_pba->FillBoundary(geom[lev].periodicity());

       ng = p0_g[lev]->nGrow();
       std::unique_ptr<MultiFab> p0_g_pba(new MultiFab(pba,pdm,p0_g[lev]->nComp(),ng));
       p0_g_pba->copy(*p0_g[lev],0,0,1,ng,ng);
       p0_g_pba->FillBoundary(p0_periodicity);

       BoxArray x_face_ba = pba;
       x_face_ba.surroundingNodes(0);
       ng = u_g[lev]->nGrow();
       std::unique_ptr<MultiFab> u_g_pba(new MultiFab(x_face_ba,pdm,u_g[lev]->nComp(),ng));
       std::unique_ptr<MultiFab> gpx    (new MultiFab(x_face_ba,pdm,1,1));
       u_g_pba->copy(*u_g[lev],0,0,1,ng,ng);
       u_g_pba->FillBoundary(geom[lev].periodicity());

       BoxArray y_face_ba = pba;
       y_face_ba.surroundingNodes(1);
       ng = v_g[lev]->nGrow();
       std::unique_ptr<MultiFab> v_g_pba(new MultiFab(y_face_ba,pdm,v_g[lev]->nComp(),ng));
       std::unique_ptr<MultiFab> gpy    (new MultiFab(y_face_ba,pdm,1,1));
       v_g_pba->copy(*v_g[lev],0,0,1,ng,ng);
       v_g_pba->FillBoundary(geom[lev].periodicity());

       BoxArray z_face_ba = pba;
       z_face_ba.surroundingNodes(2);
       ng = w_g[lev]->nGrow();
       std::unique_ptr<MultiFab> w_g_pba(new MultiFab(z_face_ba,pdm,w_g[lev]->nComp(),ng));
       std::unique_ptr<MultiFab> gpz    (new MultiFab(z_face_ba,pdm,1,1));
       w_g_pba->copy(*w_g[lev],0,0,1,ng,ng);
       w_g_pba->FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIter mfi(*p_g_pba, true); mfi.isValid(); ++mfi)
       {
           const Box& sbx = (*p_g_pba)[mfi].box();

           construct_gradp( sbx.loVect(),   sbx.hiVect(),
                            (*p_g_pba)[mfi].dataPtr(), (*p0_g_pba)[mfi].dataPtr(),
                            BL_TO_FORTRAN_ANYD((*gpx)[mfi]),
   		            BL_TO_FORTRAN_ANYD((*gpy)[mfi]),
   			    BL_TO_FORTRAN_ANYD((*gpz)[mfi]),
                            &dx, &dy, &dz,
                            bc_ilo.dataPtr(), bc_ihi.dataPtr(), bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                            bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect(),
                            &nghost
                          );
       }

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIter mfi(*p_g_pba, true); mfi.isValid(); ++mfi)
       {
           const Box& sbx = (*p_g_pba)[mfi].box();
           set_gradp_bcs ( sbx.loVect(), sbx.hiVect(),
                         BL_TO_FORTRAN_ANYD((*gpx)[mfi]),
		         BL_TO_FORTRAN_ANYD((*gpy)[mfi]),
			 BL_TO_FORTRAN_ANYD((*gpz)[mfi]),
			 bc_ilo.dataPtr(), bc_ihi.dataPtr(),
			 bc_jlo.dataPtr(), bc_jhi.dataPtr(),
			 bc_klo.dataPtr(), bc_khi.dataPtr(),
			 domain.loVect(), domain.hiVect(),
			 &nghost );
       }

       gpx->FillBoundary(geom[lev].periodicity());
       gpy->FillBoundary(geom[lev].periodicity());
       gpz->FillBoundary(geom[lev].periodicity());

       for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
       {
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();

           calc_drag_particle( BL_TO_FORTRAN_ANYD((*gpx)[pti]),
      		               BL_TO_FORTRAN_ANYD((*gpy)[pti]),
      			       BL_TO_FORTRAN_ANYD((*gpz)[pti]),
                               BL_TO_FORTRAN_ANYD((*u_g_pba)[pti]),
      		               BL_TO_FORTRAN_ANYD((*v_g_pba)[pti]),
      			       BL_TO_FORTRAN_ANYD((*w_g_pba)[pti]),
                               &np, particles.data(), &dx, &dy, &dz
                             );
       }
    }
}
