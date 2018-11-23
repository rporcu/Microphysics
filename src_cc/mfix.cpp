#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>

std::string mfix::particle_init_type   = "AsciiFile";
std::string mfix::load_balance_type    = "FixedSize";
std::string mfix::knapsack_weight_type = "RunTimeCosts";

// Define unit vectors for easily convert indeces
amrex::IntVect mfix::e_x(1,0,0);
amrex::IntVect mfix::e_y(0,1,0);
amrex::IntVect mfix::e_z(0,0,1);

int mfix::nlev = 1;

EBSupport mfix::m_eb_support_level = EBSupport::full;


mfix::~mfix ()
{
};


mfix::mfix ()
{
    // Geometry on all levels has just been defined in the AmrCore constructor

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    nlev = maxLevel() + 1;
    amrex::Print() << "Number of levels: " << nlev << std::endl;

    istep.resize(nlev, 0);

#if 0
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
    if (use_amr_ls) {
        pc = std::unique_ptr<MFIXParticleContainer> (new MFIXParticleContainer(amr_level_set.get()));
    } else {
        pc = std::unique_ptr<MFIXParticleContainer> (new MFIXParticleContainer(this));
    }

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
mfix::usr3()
{
    if (solve_fluid)
    {
       for (int lev = 0; lev < nlev; lev++)
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

    set_bc_type(bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                domain.loVect(),domain.hiVect(),
                &dx, &dy, &dz, &xlen, &ylen, &zlen, &nghost);
}

void
mfix::fill_mf_bc(int lev, MultiFab& mf)
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
		 bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(), 
                 bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
		 bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(), 
                 domain.loVect(), domain.hiVect(), &nghost);
    }
}

void mfix::mfix_calc_volume_fraction(Real& sum_vol)
{
    BL_PROFILE("mfix::mfix_calc_volume_fraction()");

    if (solve_dem)
    {
       // This re-calculates the volume fraction within the domain
       // but does not change the values outside the domain

       // This call simply deposits the particle volume onto the grid in a PIC-like manner
       for (int lev = 0; lev < nlev; lev++)
           pc->CalcVolumeFraction(*ep_g[lev], *particle_ebfactory[lev],
                                  *bc_ilo[lev],*bc_ihi[lev], 
                                  *bc_jlo[lev],*bc_jhi[lev],
                                  *bc_klo[lev],*bc_khi[lev],
                                  nghost);
    }
    else
    {
       for (int lev = 0; lev < nlev; lev++)
          ep_g[lev]->setVal(1.);
    }

    for (int lev = 0; lev < nlev; lev++)
    {

       // Now define rop_g = ro_g * ep_g
       MultiFab::Copy(*rop_g[lev], *ro_g[lev], 0, 0, 1, ro_g[lev]->nGrow());
       MultiFab::Multiply((*rop_g[lev]), (*ep_g[lev]), 0, 0, 1, rop_g[lev]->nGrow());

       // This sets the values outside walls or periodic boundaries
       fill_mf_bc(lev,*ep_g[lev]);
       fill_mf_bc(lev,*rop_g[lev]);
   }

    // Sum up all the values of ep_g[lev] 
    // HACK  -- THIS SHOULD BE a multilevel sum
    for (int lev = 0; lev < nlev; lev++)
       sum_vol = ep_g[lev]->sum();
}
