#include <mfix.H>
#include <mfix_bc_parms.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_NodalProjector.H>

#include <mfix_mf_helpers.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

using namespace amrex;

void
mfix::mfix_apply_nodal_projection (Vector< MultiFab* >& a_S_cc,
                                   Real a_time,
                                   Real a_dt,
                                   Real a_prev_dt,
                                   bool proj_2,
                                   amrex::Vector<amrex::MultiFab const*> const& density)
{
    BL_PROFILE("mfix::mfix_apply_nodal_projection");

    bool proj_for_small_dt = false;

    // If we have dropped the dt substantially for whatever reason, use a different form of the approximate
    // projection that projects (U^*-U^n + dt Gp) rather than (U^* + dt Gp)

    if (a_time > 0 && a_dt < 0.1 * a_prev_dt)
       proj_for_small_dt      = true;

    if(m_use_drag_in_projection)
      amrex::Print() << "Adding drag to vel in nodal projection\n";
    else
      amrex::Print() << "NOT adding drag to vel in nodal projection\n";

    // Create sigma
    Vector<MultiFab> sigma_mf(nlev);

    for (int lev(0); lev < nlev; ++lev)
    {
      sigma_mf[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*m_leveldata[lev]->vel_g,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Tilebox
        Box bx = mfi.tilebox();

        Array4<Real      > const&  vel_g = m_leveldata[lev]->vel_g->array(mfi);
        Array4<Real      > const&  sigma = sigma_mf[lev].array(mfi);

        Array4<Real const> const&  rho_g = density[lev]->const_array(mfi);
        Array4<Real const> const&  ep_g  = m_leveldata[lev]->ep_g->const_array(mfi);
        Array4<Real const> const&  gp    = m_leveldata[lev]->gp->const_array(mfi);
        Array4<Real const> const&  txfr  = m_leveldata[lev]->txfr->const_array(mfi);

        if(m_use_drag_in_projection) {
          amrex::ParallelFor(bx,[a_dt, proj_2, vel_g, sigma, rho_g, ep_g, gp, txfr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {

            Real const beta = txfr(i,j,k,3);

            sigma(i,j,k) = ep_g(i,j,k)/(ep_g(i,j,k)*rho_g(i,j,k) + a_dt*beta);

            if(proj_2) {
              vel_g(i,j,k,0) += sigma(i,j,k)*a_dt*gp(i,j,k,0);
              vel_g(i,j,k,1) += sigma(i,j,k)*a_dt*gp(i,j,k,1);
              vel_g(i,j,k,2) += sigma(i,j,k)*a_dt*gp(i,j,k,2);
            }

            sigma(i,j,k) *= ep_g(i,j,k);

          });

        } else {
          amrex::ParallelFor(bx,[a_dt, proj_2, vel_g, sigma, rho_g, ep_g, gp]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {

            sigma(i,j,k) = 1.0/rho_g(i,j,k);

            if (proj_2) {
              vel_g(i,j,k,0) += sigma(i,j,k)*a_dt*gp(i,j,k,0);
              vel_g(i,j,k,1) += sigma(i,j,k)*a_dt*gp(i,j,k,1);
              vel_g(i,j,k,2) += sigma(i,j,k)*a_dt*gp(i,j,k,2);
            }

            sigma(i,j,k) *= ep_g(i,j,k);

          });
        }

     }

      // Print level infos
      if (proj_for_small_dt)
        amrex::Print() << "Before projection (with small dt modification):" << std::endl;
      else
        amrex::Print() << "Before projection:" << std::endl;

      mfix_print_max_vel(lev);
      mfix_print_max_gp(lev);
      amrex::Print() << "Min and Max of ep_g "
                     << m_leveldata[lev]->ep_g->min(0) << " "
                     << m_leveldata[lev]->ep_g->max(0) << std::endl;
    }

    // Set velocities BC before projection
    mfix_set_velocity_bcs(a_time, get_vel_g(), 0);

    // Define "vel" to be U^* - U^n rather than U^*
    if (proj_for_small_dt)
    {
       mfix_set_velocity_bcs(a_time, get_vel_g_old(), 0);

       for(int lev = 0; lev <= finest_level; lev++)
          MultiFab::Saxpy(*m_leveldata[lev]->vel_g, -1.0,
                          *m_leveldata[lev]->vel_go, 0, 0, 3,
                          m_leveldata[lev]->vel_g->nGrow());
    }

    //
    // Compute epu
    //
    Vector< MultiFab* > epu(nlev);

    for (int lev(0); lev < nlev; ++lev)
    {
        // We only need one ghost cell here -- so no need to make it bigger
        epu[lev] = new MultiFab(grids[lev], dmap[lev], 3, 1, MFInfo(), *ebfactory[lev]);

        epu[lev]->setVal(1.e200);

        MultiFab::Copy(*epu[lev], *m_leveldata[lev]->vel_g, 0, 0, 3, epu[lev]->nGrow());

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

    for (int lev(0); lev <= finest_level; lev++) {
      EB_set_covered(*a_S_cc[lev], 0, a_S_cc[lev]->nComp(), 1, 0.0);
    }

    //
    // Setup the nodal projector
    //

    LPInfo info;
    info.setMaxCoarseningLevel(nodal_mg_max_coarsening_level);

    nodal_projector.reset(new NodalProjector(get_vel_g(),
                                             GetVecOfConstPtrs(sigma_mf),
                                             geom, info,
                                             a_S_cc));

    nodal_projector->setDomainBC(BC::ppe_lobc, BC::ppe_hibc);
    nodal_projector->setAlpha(GetVecOfConstPtrs(get_ep_g()));

    nodal_projector->computeRHS(get_diveu(), epu, a_S_cc);
    nodal_projector->setCustomRHS(GetVecOfConstPtrs(get_diveu()));

    nodal_projector->project(nodal_mg_rtol, nodal_mg_atol);


    // Define "vel" to be U^{n+1} rather than (U^{n+1}-U^n)
    if (proj_for_small_dt)
    {
       for(int lev = 0; lev <= finest_level; lev++)
          MultiFab::Saxpy(*m_leveldata[lev]->vel_g, 1.0,
                          *m_leveldata[lev]->vel_go, 0, 0, 3,
                          m_leveldata[lev]->vel_g->nGrow());
    }

    // Get phi and fluxes
    Vector< const amrex::MultiFab* > phi(nlev);
    Vector< const amrex::MultiFab* > gradphi(nlev);

    phi     = nodal_projector->getPhi();
    gradphi = nodal_projector->getGradPhi();

    // Compute diveu to print it out
    nodal_projector->computeRHS(get_diveu(), epu, a_S_cc);

    // Since I did not pass dt, I have to normalize here
    Real qdt(1.0/a_dt);
    for (int lev(0); lev < nlev; ++lev)
    {
        if (proj_2)
        {
            // p := phi
            MultiFab::Copy(*m_leveldata[lev]->p_g, *phi[lev], 0, 0, 1,
                           phi[lev]->nGrow());
            MultiFab::Copy(*m_leveldata[lev]->gp, *gradphi[lev], 0, 0, 3,
                           gradphi[lev]->nGrow());
            m_leveldata[lev]->p_g->mult(qdt);
            m_leveldata[lev]->gp->mult(qdt);
        }
        else
        {
            // p := p + phi/dt
            MultiFab::Saxpy(*m_leveldata[lev]->p_g, qdt, *phi[lev], 0, 0, 1,
                            phi[lev]->nGrow());
            MultiFab::Saxpy(*m_leveldata[lev]->gp, qdt, *gradphi[lev], 0, 0, 3,
                            gradphi[lev]->nGrow());
        }
    }


    //
    // This part is just to plot diveu
    //
    mfix_set_velocity_bcs(a_time, get_vel_g(), 0);

    for (int lev(0); lev < nlev; ++lev)
    {
        epu[lev]->setVal(1.e200);

        MultiFab::Copy(*epu[lev], *m_leveldata[lev]->vel_g, 0, 0, 3, epu[lev]->nGrow());

        for (int n(0); n < 3; n++)
            MultiFab::Multiply(*epu[lev], *m_leveldata[lev]->ep_g, 0, n, 1, epu[lev]->nGrow());

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

    nodal_projector->computeRHS(get_diveu(), epu, a_S_cc);

    for (int lev = nlev-1; lev > 0; lev--)
    {
        avgDown(lev-1, *m_leveldata[lev]->vel_g, *m_leveldata[lev-1]->vel_g);
        avgDown(lev-1, *m_leveldata[lev]->gp, *m_leveldata[lev-1]->gp);
    }

    // Swap ghost cells and apply BCs to velocity
    mfix_set_velocity_bcs(a_time, get_vel_g(), 0);

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
