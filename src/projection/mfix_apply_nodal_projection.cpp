#include <mfix.H>
#include <mfix_bc.H>
#include <mfix_eb.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <hydro_NodalProjector.H>

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
                                   Vector<MultiFab*      > const& vel_g_old_in,
                                   Vector<MultiFab*      > const& vel_g_in,
                                   Vector<MultiFab*      > const& p_g_in,
                                   Vector<MultiFab*      > const& gp_in,
                                   Vector<MultiFab*      > const& ep_g_in,
                                   Vector<MultiFab*      > const& txfr_in,
                                   Vector<MultiFab const*> const& density,
                                   Vector<MultiFab const*> const& eb_vel)
{
    BL_PROFILE("mfix::mfix_apply_nodal_projection");

    bool proj_for_small_dt = false;

    // If we have dropped the dt substantially for whatever reason, use a different form of the approximate
    // projection that projects (U^*-U^n + dt Gp) rather than (U^* + dt Gp)

    if (a_time > 0 && a_dt < 0.1 * a_prev_dt)
       proj_for_small_dt = true;

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
      for (MFIter mfi(*vel_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Tilebox
        Box bx = mfi.tilebox();

        Array4<Real      > const&  vel_g = vel_g_in[lev]->array(mfi);
        Array4<Real      > const&  sigma = sigma_mf[lev].array(mfi);

        Array4<Real const> const&  rho_g = density[lev]->const_array(mfi);
        Array4<Real const> const&  ep_g  = ep_g_in[lev]->const_array(mfi);
        Array4<Real const> const&  gp    = gp_in[lev]->const_array(mfi);
        Array4<Real const> const&  txfr  = txfr_in[lev]->const_array(mfi);

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

      if (m_run_type != RunType::PIC2DEM) {

        // Print level infos
        if (proj_for_small_dt)
          amrex::Print() << "Before projection (with small dt modification):" << std::endl;
        else
          amrex::Print() << "Before projection:" << std::endl;

        m_rw->mfix_print_max_vel(lev, vel_g_in, p_g_in);
        m_rw->mfix_print_max_gp(lev, gp_in);
        amrex::Print() << "Min and Max of ep_g "
                       << ep_g_in[lev]->min(0) << " "
                       << ep_g_in[lev]->max(0) << std::endl;
      }
    }

    // Set velocities BC before projection
    m_boundary_conditions.set_velocity_bcs(a_time, vel_g_in, 0);

    // Define "vel" to be U^* - U^n rather than U^*
    if (proj_for_small_dt)
    {
       m_boundary_conditions.set_velocity_bcs(a_time, vel_g_old_in, 0);

       for(int lev = 0; lev <= finest_level; lev++)
          MultiFab::Saxpy(*vel_g_in[lev], -1.0, *vel_g_old_in[lev], 0, 0, 3, vel_g_in[lev]->nGrow());
    }

    //
    // Compute epu
    //
    Vector< MultiFab* > epu(nlev);

    for (int lev(0); lev < nlev; ++lev) {
        // We only need one ghost cell here -- so no need to make it bigger
        epu[lev] = new MultiFab(grids[lev], dmap[lev], 3, 1, MFInfo(), *ebfactory[lev]);

        epu[lev]->setVal(1.e200);

        MultiFab::Copy(*epu[lev], *vel_g_in[lev], 0, 0, 3, epu[lev]->nGrow());

        for (int n(0); n < 3; n++)
            MultiFab::Multiply(*epu[lev], *ep_g_in[lev], 0, n, 1, epu[lev]->nGrow());
    }

    // Extrapolate Dirichlet values to ghost cells -- but do it differently in that
    //  no-slip walls are treated exactly like slip walls --
    // Note that this routine is essential to impose the correct inflow bc's on
    //  the product ep  * vel
    //
    //  TODO
    // Why are we using this instead of simply multiplying vel and ep with their BCs in place
    // already?
    m_boundary_conditions.set_vec_bcs(a_time, epu);

    for (int lev(0); lev < nlev; ++lev) {
        // We set these to zero because if the values in the covered cells are undefined,
        //   even though they are multiplied by zero in the divu computation, we can still get NaNs
        EB_set_covered(*epu[lev], 0, epu[lev]->nComp(), 1, 0.0);
    }

    for (int lev(0); lev <= finest_level; lev++) {
      EB_set_covered(*a_S_cc[lev], 0, a_S_cc[lev]->nComp(), 1, 0.0);
    }

    if ( m_redistribute_before_nodal_proj ) {
      PreProjectionRedistribution(a_time);
    }

    //
    // Setup the nodal projector
    //

    LPInfo info;
    info.setMaxCoarseningLevel(nodalproj_options->max_coarsening_level);
    info.setAgglomerationGridSize(agg_grid_size);

    nodal_projector = std::make_unique<Hydro::NodalProjector>(vel_g_in,
                                                       GetVecOfConstPtrs(sigma_mf),
                                                       geom,
                                                       info,
                                                       a_S_cc);

    nodalproj_options->apply(*nodal_projector);
    nodal_projector->setDomainBC(m_boundary_conditions.ppe_lobc(), m_boundary_conditions.ppe_hibc());

    // By setting alpha = ep_g, the nodal projection will correct the velocity by
    // (sigma / alpha) grad(phi) rather than sigma grad phi
    nodal_projector->setAlpha(GetVecOfConstPtrs(ep_g_in));

    if (m_embedded_boundaries.has_flow()) {
       for (int lev(0); lev < nlev; ++lev) {
          nodal_projector->getLinOp().setEBInflowVelocity(lev, *eb_vel[lev]);
       }
    }

    nodal_projector->computeRHS(get_diveu(), epu, a_S_cc);
    nodal_projector->setCustomRHS(GetVecOfConstPtrs(get_diveu()));

    nodal_projector->project(nodalproj_options->mg_rtol, nodalproj_options->mg_atol);

    // Define "vel" to be U^{n+1} rather than (U^{n+1}-U^n)
    if (proj_for_small_dt)
    {
       for(int lev = 0; lev <= finest_level; lev++)
          MultiFab::Saxpy(*vel_g_in[lev], 1.0, *vel_g_old_in[lev], 0, 0, 3, vel_g_in[lev]->nGrow());
    }

    // Get phi and fluxes
    Vector< const MultiFab* > phi(nlev);
    Vector< const MultiFab* > gradphi(nlev);

    phi     = nodal_projector->getPhiConst();
    gradphi = nodal_projector->getGradPhiConst();

    // Since I did not pass dt, I have to normalize here
    Real qdt(1.0/a_dt);
    for (int lev(0); lev < nlev; ++lev)
    {
        if (proj_2)
        {
            // p := phi/dt
            MultiFab::Copy(*p_g_in[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
            MultiFab::Copy(*gp_in[lev], *gradphi[lev], 0, 0, 3, gradphi[lev]->nGrow());
            p_g_in[lev]->mult(qdt);
            gp_in[lev]->mult(qdt);
        }
        else
        {
            // p := p + phi/dt
            MultiFab::Saxpy(*p_g_in[lev], qdt, *phi[lev], 0, 0, 1, phi[lev]->nGrow());
            MultiFab::Saxpy(*gp_in[lev], qdt, *gradphi[lev], 0, 0, 3, gradphi[lev]->nGrow());
        }
    }

    // Perform the redistribution operation on the updated (projected) velocity field -- and
    //     update gp to maintain consistency
    // This has been currently disabled since it seems to cause magnification of differences on
    // different grids. Needs to be revisited.
    if ( m_redistribute_nodal_proj ) {
      PostProjectionRedistribution(a_time, a_dt, GetVecOfPtrs(sigma_mf));
    }

    // Compute diveu for diagnostics only
    PostProjectionDiagnostics(a_time, epu, vel_g_in, p_g_in, gp_in, ep_g_in, a_S_cc, proj_for_small_dt);

    for (int lev(0); lev < nlev; lev++)
      delete epu[lev];
}

void
mfix::PostProjectionDiagnostics(Real a_time,
                                Vector<MultiFab*> const& epu,
                                Vector<MultiFab*> const& vel_g_in,
                                Vector<MultiFab*> const& p_g_in,
                                Vector<MultiFab*> const& gp_in,
                                Vector<MultiFab*> const& ep_g_in,
                                Vector<MultiFab*> const& a_S_cc,
                                bool proj_for_small_dt)
{
    //
    // This part is just to plot diveu
    //
    m_boundary_conditions.set_velocity_bcs(a_time, vel_g_in, 0);

    for (int lev(0); lev < nlev; ++lev) {
        epu[lev]->setVal(1.e200);

        MultiFab::Copy(*epu[lev], *vel_g_in[lev], 0, 0, 3, epu[lev]->nGrow());

        for (int n(0); n < 3; n++)
            MultiFab::Multiply(*epu[lev], *ep_g_in[lev], 0, n, 1, epu[lev]->nGrow());
    }

    // Extrapolate Dirichlet values to ghost cells -- but do it differently in that
    //  no-slip walls are treated exactly like slip walls --
    // Note that this routine is essential to impose the correct inflow bc's on
    //  the product ep  * vel
    //
    // TODO
    // Why are we using this instead of simply multiplying vel and ep with their BCs in place
    // already?
    m_boundary_conditions.set_vec_bcs(a_time, epu);

    for (int lev(0); lev < nlev; ++lev) {
        // We set these to zero because if the values in the covered cells are undefined,
        //   even though they are multiplied by zero in the divu computation, we can still get NaNs
        EB_set_covered(*epu[lev], 0, epu[lev]->nComp(), 1, 0.0);
    }

    nodal_projector->computeRHS(get_diveu(), epu, a_S_cc);

    for (int lev = nlev-1; lev > 0; lev--)
    {
        avgDown(lev-1, *vel_g_in[lev], *vel_g_in[lev-1]);
        avgDown(lev-1, *gp_in[lev], *gp_in[lev-1]);
    }

    // Swap ghost cells and apply BCs to velocity
    m_boundary_conditions.set_velocity_bcs(a_time, vel_g_in, 0);

    if (m_run_type != RunType::PIC2DEM) {

      // Print level info after projection
      for (int lev(0); lev < nlev; lev++)
      {
          if (proj_for_small_dt)
             amrex::Print() << "After  projection (with small dt modification):" << std::endl;
          else
             amrex::Print() << "After  projection:" << std::endl;

          m_rw->mfix_print_max_vel(lev, vel_g_in, p_g_in);
      }
    }
}
