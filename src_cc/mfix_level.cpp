#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>

std::string mfix_level::particle_init_type   = "AsciiFile";
std::string mfix_level::load_balance_type    = "FixedSize";
std::string mfix_level::knapsack_weight_type = "RunTimeCosts";

// Define unit vectors for easily convert indeces
amrex::IntVect mfix_level::e_x(1,0,0);
amrex::IntVect mfix_level::e_y(0,1,0);
amrex::IntVect mfix_level::e_z(0,0,1);

int mfix_level::m_eb_basic_grow_cells = 2;
int mfix_level::m_eb_volume_grow_cells = 2;
int mfix_level::m_eb_full_grow_cells = 2;
EBSupport mfix_level::m_eb_support_level = EBSupport::full;


mfix_level::~mfix_level ()
{
};


mfix_level::mfix_level ()
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
mfix_level::ResizeArrays ()
{
    int nlevs_max = maxLevel() + 1;

    // Particle Container
    pc = std::unique_ptr<MFIXParticleContainer> (new MFIXParticleContainer(this));

    // HACK: temporary flag used to turn on legacy mode
    //   (used in evlove particles)
    pc -> legacy__eb_collisions = legacy__eb_collisions;

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

    // RHS and solution arrays for diffusive solve
    rhs_diff.resize(nlevs_max);
    phi_diff.resize(nlevs_max);

    // Current (vel_g) and old (vel_go) velocities
    vel_g.resize(nlevs_max);
    vel_go.resize(nlevs_max);

    // Pressure gradients
    gp.resize(nlevs_max);
    gp0.resize(nlevs_max);

    f_gds.resize(nlevs_max);
    drag.resize(nlevs_max);

    mu_g.resize(nlevs_max);
    lambda_g.resize(nlevs_max);
    trD_g.resize(nlevs_max);

    // Vorticity
    vort.resize(nlevs_max);

    // MAC velocities used for defining convective term
    m_u_mac.resize(nlevs_max);
    m_v_mac.resize(nlevs_max);
    m_w_mac.resize(nlevs_max);

    // MultiFab storing level-set data
    ls.resize(nlevs_max);

    xslopes.resize(nlevs_max);
    yslopes.resize(nlevs_max);
    zslopes.resize(nlevs_max);

    bcoeff.resize(nlevs_max);
    for (int i = 0; i < nlevs_max; ++i ) {
        bcoeff[i].resize(3);
    }

    bcoeff_diff.resize(nlevs_max);
    for (int i = 0; i < nlevs_max; ++i ) {
        bcoeff_diff[i].resize(3);
    }

    if (solve_dem)
       particle_cost.resize(nlevs_max);
    if (solve_fluid)
       fluid_cost.resize(nlevs_max);

    // EB factory
    ebfactory.resize(nlevs_max);
    particle_ebfactory.resize(nlevs_max);
}

void
mfix_level::usr3(int lev)
{
    if (solve_fluid) 
    {
       Real dx = geom[lev].CellSize(0);
       Real dy = geom[lev].CellSize(1);
       Real dz = geom[lev].CellSize(2);

       // We deliberately don't tile this loop 
       for (MFIter mfi(*p_g[lev]); mfi.isValid(); ++mfi)
       {
          const Box& sbx = (*p_g[lev])[mfi].box();
          const Box& ubx = (*vel_g[lev])[mfi].box();
   
          mfix_usr3((*vel_g[lev])[mfi].dataPtr(0), ubx.loVect(), ubx.hiVect(),
                    (*vel_g[lev])[mfi].dataPtr(1), ubx.loVect(), ubx.hiVect(),
                    (*vel_g[lev])[mfi].dataPtr(2), ubx.loVect(), ubx.hiVect(),
                    (*p_g[lev])[mfi].dataPtr(), sbx.loVect(), sbx.hiVect(),
                     &dx, &dy, &dz);
       }
    }
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
                &dx, &dy, &dz, &xlen, &ylen, &zlen, &nghost);
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
		 bc_klo.dataPtr(), bc_khi.dataPtr(), domain.loVect(), domain.hiVect(), &nghost);
    }
}

void mfix_level::mfix_calc_volume_fraction(int lev, Real& sum_vol)
{
    BL_PROFILE("mfix_level::mfix_calc_volume_fraction()");

    if (solve_dem)
    {
       // This re-calculates the volume fraction within the domain
       // but does not change the values outside the domain

       // This call simply deposits the particle volume onto the grid in a PIC-like manner
       pc->CalcVolumeFraction(*ep_g[lev],bc_ilo,bc_ihi,bc_jlo,bc_jhi,bc_klo,bc_khi,nghost);
    }
    else
    {
       if (use_epg_hack)
       {
          Box domain(geom[lev].Domain());
          for (MFIter mfi(*ep_g[lev],false); mfi.isValid(); ++mfi)
          {
              set_epg( BL_TO_FORTRAN_ANYD((*ep_g[lev])[mfi]),
                       domain.loVect (), domain.hiVect () );
          }
       } else {
          ep_g[lev]->setVal(1.);
       }
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

           calc_particle_beta(
               sbx.loVect(), sbx.hiVect(),
               (*ep_g[lev])[pti].dataPtr() , (*ro_g[lev])[pti].dataPtr(),
               (*vel_g[lev])[pti].dataPtr(), (*mu_g[lev])[pti].dataPtr(),
               &np, particles.data(), &dx, &dy, &dz);
       }

       // ******************************************************************************
       // Now use the beta of individual particles to create the drag terms on the fluid 
       // ******************************************************************************

       drag[lev]->setVal(0.0L);
       f_gds[lev]->setVal(0.0L);

       pc -> CalcDragOnFluid(*f_gds[lev], *drag[lev],
                             bc_ilo,bc_ihi,bc_jlo,bc_jhi,bc_klo,bc_khi,nghost);
    }
    else 
    {

       BoxArray            pba = pc->ParticleBoxArray(lev);
       DistributionMapping pdm = pc->ParticleDistributionMap(lev);

       // Temporary arrays
       int ng = ep_g[lev]->nGrow();
       std::unique_ptr<MultiFab> ep_g_pba(new MultiFab(pba,pdm,ep_g[lev]->nComp(),ng));
       ep_g_pba->copy(*ep_g[lev],0,0,1,ng,ng);
       ep_g_pba->FillBoundary(geom[lev].periodicity());

       ng = ro_g[lev]->nGrow();
       std::unique_ptr<MultiFab> ro_g_pba(new MultiFab(pba,pdm,ro_g[lev]->nComp(),ng));
       ro_g_pba->copy(*ro_g[lev],0,0,1,ng,ng);
       ro_g_pba->FillBoundary(geom[lev].periodicity());

       ng = mu_g[lev]->nGrow();
       std::unique_ptr<MultiFab> mu_g_pba(new MultiFab(pba,pdm,mu_g[lev]->nComp(),ng));
       mu_g_pba->copy(*mu_g[lev],0,0,1,ng,ng);
       mu_g_pba->FillBoundary(geom[lev].periodicity());

       ng = vel_g[lev]->nGrow();
       std::unique_ptr<MultiFab> vel_g_pba(new MultiFab(pba,pdm,vel_g[lev]->nComp(),ng));
       vel_g_pba->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),ng,ng);
       vel_g_pba->FillBoundary(geom[lev].periodicity());

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

           calc_particle_beta(
               sbx.loVect(), sbx.hiVect(),
               (*ep_g_pba)[pti].dataPtr(), (*ro_g_pba)[pti].dataPtr(),
               (*vel_g_pba)[pti].dataPtr(),  (*mu_g_pba)[pti].dataPtr(),
               &np, particles.data(), &dx, &dy, &dz);
       }

       // ******************************************************************************
       // Now use the beta of individual particles to create the drag terms on the fluid 
       // ******************************************************************************

       std::unique_ptr<MultiFab>  drag_pba(new MultiFab(pba,pdm,drag[lev]->nComp(),drag[lev]->nGrow()));
       std::unique_ptr<MultiFab> f_gds_pba(new MultiFab(pba,pdm,f_gds[lev]->nComp(),f_gds[lev]->nGrow()));

       f_gds_pba->setVal(0.0L);
       drag_pba->setVal(0.0L);

       pc -> CalcDragOnFluid(*f_gds_pba,*drag_pba,
                             bc_ilo,bc_ihi,bc_jlo,bc_jhi,bc_klo,bc_khi,nghost);

       // Copy back from the dual grids.
        drag[lev] ->copy(*drag_pba);
       f_gds[lev] ->copy(*f_gds_pba);

    } // if not OnSameGrids

    // The projection method uses drag to update u, not (cell_vol * u), so we must divide by vol here
    //     and we will divide by density in the update.
    Real ovol = 1./(dx*dy*dz);
     drag[lev]->mult(ovol);
    f_gds[lev]->mult(ovol);

    // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
     drag[lev]->FillBoundary(geom[lev].periodicity());
    f_gds[lev]->FillBoundary(geom[lev].periodicity());
}

void
mfix_level::mfix_calc_drag_particle(int lev)
{
    BL_PROFILE("mfix_level::mfix_calc_drag_particle()");

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    Box domain(geom[lev].Domain());

    if (OnSameGrids) 
    {
       MultiFab gp_tmp, gp0_tmp;

       gp_tmp.define(grids[lev],dmap[lev],3,1);
       gp_tmp.copy(*gp[lev],0,0,3,1,1);
       gp_tmp.FillBoundary(geom[lev].periodicity());

       //
       // NOTE -- it is essential that we call set_gradp_bcs after calling FillBoundary
       //         because the set_gradp_bcs call hopefully sets the ghost cells exterior 
       //         to the domain from ghost cells interior to the domain
       //
#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIter mfi(gp_tmp, true); mfi.isValid(); ++mfi)
       {
           const Box& bx = mfi.tilebox();
           set_gradp_bcs ( bx.loVect(), bx.hiVect(),
                           BL_TO_FORTRAN_ANYD(gp_tmp[mfi]),
                           bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                           bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                           bc_klo.dataPtr(), bc_khi.dataPtr(),
                           domain.loVect(), domain.hiVect(),
                           &nghost);
       }

       // Extrapolate velocity Dirichlet bc's to ghost cells
       int extrap_dir_bcs = 1; 
       mfix_set_velocity_bcs(lev, extrap_dir_bcs);

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
       {
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();
   
           calc_drag_particle(
                           BL_TO_FORTRAN_ANYD(       gp_tmp[pti]),
                           BL_TO_FORTRAN_ANYD((  *gp0[lev])[pti]),
                           BL_TO_FORTRAN_ANYD((*vel_g[lev])[pti]),
                           &np, particles.data(), &dx, &dy, &dz);
       }

       // Reset velocity Dirichlet bc's to face values
       extrap_dir_bcs = 0; 
       mfix_set_velocity_bcs(lev, extrap_dir_bcs);
    }
#if 0
    else 
    {

       BoxArray            pba = pc->ParticleBoxArray(lev);
       DistributionMapping pdm = pc->ParticleDistributionMap(lev);

       ng = gp[lev]->nGrow();
       std::unique_ptr<MultiFab> gp_tmp(new MultiFab(pba,pdm,gp[lev]->nComp(),ng));
       gp_tmp->copy(*gp[lev],0,0,gp[lev]->nComp(),ng,ng);
       gp_tmp->FillBoundary(geom[lev].periodicity());

       ng = gp0[lev]->nGrow();
       std::unique_ptr<MultiFab> gp0_tmp(new MultiFab(pba,pdm,gp0[lev]->nComp(),ng));
       gp0_tmp->copy(*gp0[lev],0,0,gp0[lev]->nComp(),ng,ng);
       gp0_tmp->FillBoundary(geom[lev].periodicity());

       ng = vel_g[lev]->nGrow();
       std::unique_ptr<MultiFab> vel_g_tmp(new MultiFab(tmp,pdm,vel_g[lev]->nComp(),ng));
       vel_g_tmp->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),ng,ng);
       vel_g_tmp->FillBoundary(geom[lev].periodicity());

       //
       // NOTE -- it is essential that we call set_gradp_bcs after calling FillBoundary
       //         because the set_gradp_bcs call hopefully sets the ghost cells exterior 
       //         to the domain from ghost cells interior to the domain
       //
#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIter mfi(gp_tmp, true); mfi.isValid(); ++mfi)
       {
           set_gradp_bcs ( bx.loVect(), bx.hiVect(),
                           BL_TO_FORTRAN_ANYD(gp_tmp[mfi]),
                           bc_ilo.dataPtr(), bc_ihi.dataPtr(),
                           bc_jlo.dataPtr(), bc_jhi.dataPtr(),
                           bc_klo.dataPtr(), bc_khi.dataPtr(),
                           domain.loVect(), domain.hiVect(),
                           &nghost);
       }

       int extrap_dir_bcs = 1; 
       mfix_set_velocity_bcs(lev, extrap_dir_bcs);

#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
       {
           auto& particles = pti.GetArrayOfStructs();
           const int np = particles.size();
   
           calc_drag_particle(
                           BL_TO_FORTRAN_ANYD(      gp_tmp[pti]),
                           BL_TO_FORTRAN_ANYD((   *gp0_tmp[pti]),
                           BL_TO_FORTRAN_ANYD((*vel_g_tmp)[pti]),
                           &np, particles.data(), &dx, &dy, &dz);
       }
    }
#endif
}


//
// Subroutine to compute norm0 of EB multifab
//
Real
mfix_level::mfix_norm0 ( const Vector< std::unique_ptr<MultiFab>>& mf, int lev, int comp )
{
   MultiFab mf_tmp( mf[lev]->boxArray(), mf[lev]->DistributionMap(), mf[lev]->nComp(),
                    0,  MFInfo(), *ebfactory[lev]);
  
   MultiFab::Copy( mf_tmp, *mf[lev], comp, comp, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );
   
   return mf_tmp.norm0( comp );
}

Real
mfix_level::mfix_norm0 ( MultiFab& mf, int lev, int comp )
{
   MultiFab mf_tmp( mf.boxArray(), mf.DistributionMap(), mf.nComp(),
                    0,  MFInfo(), *ebfactory[lev]);
  
   MultiFab::Copy( mf_tmp, mf, comp, comp, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );
   
   return mf_tmp.norm0( comp );
}


//
// Subroutine to compute norm1 of EB multifab
//
Real
mfix_level::mfix_norm1 ( const Vector< std::unique_ptr<MultiFab>>& mf, int lev, int comp )
{
   MultiFab mf_tmp( mf[lev]->boxArray(), mf[lev]->DistributionMap(), mf[lev]->nComp(),
                    0,  MFInfo(), *ebfactory[lev]);
  
   MultiFab::Copy( mf_tmp, *mf[lev], comp, comp, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );
   
   return mf_tmp.norm1( comp, geom[lev].periodicity() );
}

Real
mfix_level::mfix_norm1 ( MultiFab& mf, int lev, int comp )
{
   MultiFab mf_tmp( mf.boxArray(), mf.DistributionMap(), mf.nComp(),
                    0,  MFInfo(), *ebfactory[lev]);
  
   MultiFab::Copy( mf_tmp, mf, comp, comp, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );
   
   return mf_tmp.norm1( comp, geom[lev].periodicity() );
}

void
mfix_level::mfix_compute_vort (int lev ) 
{
    BL_PROFILE("mfix_level::mfix_compute_vort");
    Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*vel_g[lev],true); mfi.isValid(); ++mfi) 
    {
       // Tilebox
       Box bx = mfi.tilebox ();

       // This is to check efficiently if this tile contains any eb stuff
       const EBFArrayBox&  vel_fab = dynamic_cast<EBFArrayBox const&>((*vel_g[lev])[mfi]);
       const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();
 
       if (flags.getType(amrex::grow(bx,0)) == FabType::regular )
       {
         compute_vort (
                     BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_ANYD((* vort[lev])[mfi]),
                     BL_TO_FORTRAN_ANYD((*vel_g[lev])[mfi]),
                     geom[lev].CellSize());
       } else {
          vort[lev]->setVal( 0.0, bx, 0, 1);
       }
    }
}
