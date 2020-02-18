#include <mfix_F.H>
#include <mfix.H>
#include <MFIX_BC_Parms.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_NodalProjector.H>

#include <MFIX_MFHelpers.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

using namespace amrex;

void
mfix::mfix_apply_nodal_projection (Vector< MultiFab* >& a_depdt,
                                   Real a_time,
                                   Real a_dt,
                                   bool proj_2 )
{
    BL_PROFILE("mfix::mfix_apply_nodal_projection");

    bool proj_for_small_dt = false;

    // If we have dropped the dt substantially for whatever reason, use a different form of the approximate
    // projection that projects (U^*-U^n + dt Gp) rather than (U^* + dt Gp)

    if (a_time > 0 && a_dt < 0.1 * prev_dt)
       proj_for_small_dt      = true;

    for (int lev(0); lev < nlev; ++lev)
    {
        // Here we add (dt * (1/rho gradp)) to ustar
        if (proj_2)
        {
            // Convert velocities to momenta
            for (int n(0); n < 3; ++n)
                MultiFab::Multiply(*vel_g[lev], *ro_g[lev] , 0, n, 1, vel_g[lev]->nGrow() );

            MultiFab::Saxpy(*vel_g[lev], a_dt, *gp[lev], 0, 0, 3, vel_g[lev]->nGrow());

            // Convert momenta back to velocities
            for (int n(0); n < 3; n++)
                MultiFab::Divide(*vel_g[lev], *ro_g[lev], 0, n, 1, vel_g[lev]->nGrow() );
        }

        // Print level infos
        if (proj_for_small_dt)
           amrex::Print() << "Before projection (with small dt modification):" << std::endl;
        else
           amrex::Print() << "Before projection:" << std::endl;

        mfix_print_max_vel(lev);
        mfix_print_max_gp(lev);
        amrex::Print() << "Min and Max of ep_g "
                       << (m_leveldata[lev]->ep_g)->min(0) << " "
                       << (m_leveldata[lev]->ep_g)->max(0) << std::endl;
  }

    // Set velocities BC before projection
    mfix_set_velocity_bcs(a_time, vel_g, 0);

    // Define "vel" to be U^* - U^n rather than U^*
    if (proj_for_small_dt)
    {
       mfix_set_velocity_bcs(a_time, vel_go, 0);

       for(int lev = 0; lev <= finest_level; lev++)
          MultiFab::Saxpy(*vel_g[lev], -1.0, *vel_go[lev], 0, 0, AMREX_SPACEDIM, vel_g[lev]->nGrow());
    }

    //
    // Compute epu
    //
    Vector< MultiFab* > epu(nlev);

    for (int lev(0); lev < nlev; ++lev)
    {
        // We only need one ghost cell here -- so no need to make it bigger
        int nghost(1);
        epu[lev] = new MultiFab(grids[lev], dmap[lev], 3, 1 , MFInfo(), *ebfactory[lev]);

        epu[lev]->setVal(1.e200);

        MultiFab::Copy(*epu[lev], *vel_g[lev], 0, 0, 3, epu[lev]->nGrow());

        for (int n(0); n < 3; n++)
          MultiFab::Multiply(*epu[lev], *(m_leveldata[lev]->ep_g), 0, n, 1, epu[lev]->nGrow());

        epu[lev]->FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        // Extrapolate Dirichlet values to ghost cells -- but do it differently in that
        //  no-slip walls are treated exactly like slip walls --
        // Note that this routine is essential to impose the correct inflow bc's on
        //  the product ep  * vel
        for (MFIter mfi(*epu[lev], false); mfi.isValid(); ++mfi)
        {
            // Why are we using this instead of simply multiplying vel and ep with their BCs in place
            // already?
            set_vec_bcs(lev, (*epu[lev])[mfi], geom[lev].Domain());
        }

        epu[lev]->FillBoundary(geom[lev].periodicity());

        // We set these to zero because if the values in the covered cells are undefined,
        //   even though they are multiplied by zero in the divu computation, we can still get NaNs
        EB_set_covered(*epu[lev], 0, epu[lev]->nComp(), 1, 0.0);
    }

    //
    // Create sigma
    //
    Vector<MultiFab> sigma(nlev);

    for (int lev = 0; lev < nlev; ++lev )
    {
        sigma[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
        MultiFab::Copy(sigma[lev], *(m_leveldata[lev]->ep_g), 0, 0, 1, 0);
        MultiFab::Divide(sigma[lev], *ro_g[lev], 0, 0, 1, 0);
    }

    //
    // Setup the nodal projector
    //
    Box domain(geom[0].Domain());

    LPInfo info;
    info.setMaxCoarseningLevel(nodal_mg_max_coarsening_level);

    Vector<std::unique_ptr<MultiFab>> ep_g_vec(m_leveldata.size());
    for (int lev = 0; lev < m_leveldata.size(); ++lev)
    {
      MultiFab& ep_g = *(m_leveldata[lev]->ep_g);
      ep_g_vec[lev].reset(new MultiFab(ep_g.boxArray(), ep_g.DistributionMap(),
                                       ep_g.nComp(), ep_g.nGrow(), MFInfo(), ep_g.Factory()));
      MultiFab::Copy(*ep_g_vec[lev], ep_g, 0, 0, ep_g.nComp(), ep_g.nGrow());
    }

    nodal_projector.reset(new NodalProjector(vel_g, GetVecOfConstPtrs(sigma), geom, info));
    nodal_projector->setDomainBC(BC::ppe_lobc, BC::ppe_hibc);
    nodal_projector->setAlpha(GetVecOfConstPtrs(ep_g_vec));

    nodal_projector->computeRHS(diveu, epu, a_depdt);
    nodal_projector->setCustomRHS(GetVecOfConstPtrs(diveu));

    nodal_projector->project();


    // Define "vel" to be U^{n+1} rather than (U^{n+1}-U^n)
    if (proj_for_small_dt)
    {
       for(int lev = 0; lev <= finest_level; lev++)
          MultiFab::Saxpy(*vel_g[lev], 1.0, *vel_go[lev], 0, 0, AMREX_SPACEDIM, vel_g[lev]->nGrow());
    }

    // Get phi and fluxes
    Vector< const amrex::MultiFab* > phi(nlev);
    Vector< const amrex::MultiFab* > gradphi(nlev);

    phi     = nodal_projector->getPhi();
    gradphi = nodal_projector->getGradPhi();

    // Compute diveu to print it out
    nodal_projector->computeRHS(diveu, epu, a_depdt);

    // Since I did not pass dt, I have to normalize here
    Real qdt(1.0/a_dt);
    for (int lev(0); lev < nlev; ++lev)
    {
        if (proj_2)
        {
            // p := phi
            MultiFab::Copy(*p_g[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
            MultiFab::Copy( *gp[lev], *gradphi[lev], 0, 0, 3, gradphi[lev]->nGrow());
            p_g[lev]->mult(qdt);
            gp[lev]->mult(qdt);
        }
        else
        {
            // p := p + phi
            MultiFab::Saxpy(*p_g[lev], qdt, *phi[lev], 0, 0, 1, phi[lev]->nGrow());
            MultiFab::Saxpy(*gp[lev],  qdt, *gradphi[lev], 0, 0, 3, gradphi[lev]->nGrow());
        }
    }


    //
    // This part is just to plot diveu
    //
    mfix_set_velocity_bcs(a_time, vel_g, 0);

    for (int lev(0); lev < nlev; ++lev)
    {
        epu[lev]->setVal(1.e200);

        MultiFab::Copy(*epu[lev], *vel_g[lev], 0, 0, 3, epu[lev]->nGrow());

        for (int n(0); n < 3; n++)
            MultiFab::Multiply(*epu[lev], *(m_leveldata[lev]->ep_g), 0, n, 1, epu[lev]->nGrow());

        epu[lev]->FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        // Extrapolate Dirichlet values to ghost cells -- but do it differently in that
        //  no-slip walls are treated exactly like slip walls --
        // Note that this routine is essential to impose the correct inflow bc's on
        //  the product ep  * vel
        for (MFIter mfi(*epu[lev], false); mfi.isValid(); ++mfi)
        {
            // Why are we using this instead of simply multiplying vel and ep with their BCs in place
            // already?
            set_vec_bcs(lev, (*epu[lev])[mfi], geom[lev].Domain() );
        }

        epu[lev]->FillBoundary(geom[lev].periodicity());

        // We set these to zero because if the values in the covered cells are undefined,
        //   even though they are multiplied by zero in the divu computation, we can still get NaNs
        EB_set_covered(*epu[lev], 0, epu[lev]->nComp(), 1, 0.0);
    }

    nodal_projector->computeRHS(diveu, epu, a_depdt);

    for (int lev = nlev-1; lev > 0; lev--)
    {
        avgDown(lev-1, *vel_g[lev], *vel_g[lev-1]);
        avgDown(lev-1, *gp[lev],    *gp[lev-1]);
    }

    // Swap ghost cells and apply BCs to velocity
    mfix_set_velocity_bcs(a_time, vel_g, 0);

    // Print level info after projection
    for (int lev(0); lev < nlev; lev++)
    {
        if (proj_for_small_dt)
           amrex::Print() << "After  projection (with small dt modification):" << std::endl;
        else
           amrex::Print() << "After  projection:" << std::endl;

        mfix_print_max_vel(lev);
    }

    for (int lev(0); lev < nlev; lev++)
      delete epu[lev];
}
