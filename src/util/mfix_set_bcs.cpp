#include <mfix.H>
#include <mfix_fluid.H>


using namespace amrex;


void
MFIXBoundaryConditions::set_epg_bcs (Real time,
                                     const Vector< MultiFab* >& ep_g_in,
                                     const int dir_bc)
{
  BL_PROFILE("MFIXBoundaryConditions::set_ep_g_bcs()");

  const int nlev = ep_g_in.size();

  for (int lev = 0; lev < nlev; lev++) {
    // Set all values outside the domain to covered_val just to avoid use of
    // undefined
    ep_g_in[lev]->setDomainBndry(mfix::covered_val, m_geom[lev]);
  }

  {
    Real* p_bc_ep_g = m_bc_ep_g.data();

    auto bcs_function = [p_bc_ep_g]
      AMREX_GPU_DEVICE (const int bct,
                        const int bcv,
                        const IntVect& ijk,
                        const IntVect& dom_ijk,
                        const int n,
                        const Array4<Real>& mf_arr,
                        const Array4<const EBCellFlag>& /*flags_arr*/,
                        const int /*dir*/)
    {
      if(bct == BCList::pinf || bct == BCList::pout) {
        mf_arr(ijk,n) = mf_arr(dom_ijk,n);
      }
      else if (bct == BCList::minf) {
        mf_arr(ijk,n) = p_bc_ep_g[bcv];
      }
    };

    set_bcs(time, bcs_function, ep_g_in);
  }

  {
    auto bcs_function = [dir_bc]
      AMREX_GPU_DEVICE (const int bct,
                        const IntVect& ijk,
                        const IntVect& dom_ijk,
                        const IntVect& near_ijk,
                        const int n,
                        const Array4<Real>& mf_arr)
    {
      if (bct == BCList::minf) {
        if (dir_bc == 1) {
          mf_arr(ijk,n) = 2*mf_arr(ijk,n) - mf_arr(near_ijk,n);
        } else if (dir_bc == 2) {
          mf_arr(ijk,n) = mf_arr(dom_ijk,n);
        }
      }
    };

    set_bcs_2D(time, bcs_function, ep_g_in);
  }

  for (int lev = 0; lev < nlev; lev++) {
    EB_set_covered(*ep_g_in[lev], 0, ep_g_in[lev]->nComp(), ep_g_in[lev]->nGrow(),
        mfix::covered_val);

    // Do this after as well as before to pick up terms that got updated in the
    // call above
    ep_g_in[lev]->FillBoundary(m_geom[lev].periodicity());
  }
}


void
MFIXBoundaryConditions::set_temperature_bcs (Real time,
                                             MFIXFluidPhase& fluid,
                                             Vector< MultiFab* > const& T_g_in)
{
  BL_PROFILE("MFIXBoundaryConditions::set_temperature_bcs()");

  set_temperature_bc_values(time, fluid);
  Real* p_bc_t_g = m_bc_t_g.data();

  auto bcs_function = [p_bc_t_g]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& /*flags_arr*/,
                      const int /*dir*/)
  {
    if(bct == BCList::pout) {
      mf_arr(ijk,n) = mf_arr(dom_ijk,n);
    }
    else if (bct == BCList::minf || bct == BCList::pinf) {
      mf_arr(ijk,n) = p_bc_t_g[bcv];
    }
  };

  set_bcs(time, bcs_function, T_g_in);
}


void
MFIXBoundaryConditions::set_enthalpy_bcs (Real time,
                                             MFIXFluidPhase& fluid,
                                             Vector< MultiFab* > const& h_g_in)
{
  BL_PROFILE("MFIXBoundaryConditions::set_enthalpy_bcs()");

  const int nspecies_g = fluid.nspecies();
  const int fluid_is_a_mixture = fluid.isMixture();

  set_temperature_bc_values(time, fluid);
  Real* p_bc_t_g = m_bc_t_g.data();

  if(fluid_is_a_mixture)
    set_species_bc_values(time, fluid);

  Real** p_bc_X_gk = fluid_is_a_mixture ? m_bc_X_gk_ptr.data() : nullptr;

  const auto& fluid_parms = fluid.parameters();

  auto bcs_function = [p_bc_t_g,p_bc_X_gk,fluid_parms,fluid_is_a_mixture,nspecies_g]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& flags_arr,
                      const int /*dir*/)
  {
    if(bct == BCList::pout) {
      mf_arr(ijk,n) = mf_arr(dom_ijk,n);
    }
    else if (bct == BCList::minf || bct == BCList::pinf) {
      const int cell_is_covered = static_cast<int>(flags_arr(ijk).isCovered());

      if (!fluid_is_a_mixture) {
        mf_arr(ijk,n) = fluid_parms.calc_h_g<run_on>(p_bc_t_g[bcv], cell_is_covered);
      } else {
        Real h_g_sum(0);
        for (int n_g(0); n_g < nspecies_g; n_g++) {
          const Real h_gk = fluid_parms.calc_h_gk<run_on>(p_bc_t_g[bcv], n_g, cell_is_covered);
          h_g_sum += p_bc_X_gk[n_g][bcv]*h_gk;
        }
        mf_arr(ijk,n) = h_g_sum;
      }
    }
  };

  set_bcs(time, bcs_function, h_g_in);
}


void
MFIXBoundaryConditions::set_velocity_bcs (Real time,
                                          Vector< MultiFab* > const& vel_g_in,
                                          int extrap_dir_bcs)
{
  BL_PROFILE("MFIXBoundaryConditions::set_velocity_bcs()");

  // Set all values outside the domain to covered_val just to avoid use of
  // undefined
  const int nlev = vel_g_in.size();
  for (int lev(0); lev < nlev; ++lev) {
    vel_g_in[lev]->setDomainBndry(mfix::covered_val, m_geom[lev]);
  }

  set_velocity_bc_values(time);

  mfix::mfix_usr1(time);

  const Real* p_bc_u_g = m_bc_u_g.data();
  const Real* p_bc_v_g = m_bc_v_g.data();
  const Real* p_bc_w_g = m_bc_w_g.data();

  {
    auto bcs_function = [p_bc_u_g,p_bc_v_g,p_bc_w_g]
      AMREX_GPU_DEVICE (const int bct,
                        const int bcv,
                        const IntVect& ijk,
                        const IntVect& dom_ijk,
                        const int n,
                        const Array4<Real>& mf_arr,
                        const Array4<const EBCellFlag>& /*flags_arr*/,
                        const int dir)
    {
      GpuArray<const Real*,3> p_bc_vel_g = {p_bc_u_g, p_bc_v_g, p_bc_w_g};

      if(bct == BCList::pinf || bct == BCList::pout) {
        mf_arr(ijk,n) = mf_arr(dom_ijk,n);
      }
      else if (bct == BCList::minf) {
        if (n == dir) {
          mf_arr(ijk,n) = p_bc_vel_g[dir][bcv];
        } else {
          mf_arr(ijk,n) = 0;
        }
      }
    };

    set_bcs(time, bcs_function, vel_g_in);
  }

  if (extrap_dir_bcs > 0) {

    auto bcs_function = []
      AMREX_GPU_DEVICE (const int bct,
                        const IntVect& ijk,
                        const IntVect& /*dom_ijk*/,
                        const IntVect& near_ijk,
                        const int n,
                        const Array4<Real>& mf_arr)
    {
      if(bct == BCList::minf) {
        mf_arr(ijk,n) = 2*mf_arr(ijk,n) - mf_arr(near_ijk,n);
      }
    };

    set_bcs_2D(time, bcs_function, vel_g_in);

  }
}


void
MFIXBoundaryConditions::set_vec_bcs (Real time,
                                     Vector< MultiFab* > const& mf_in)
{
  BL_PROFILE("MFIXBoundaryConditions::set_vec_bcs()");

  const Real* p_bc_u_g = m_bc_u_g.data();
  const Real* p_bc_v_g = m_bc_v_g.data();
  const Real* p_bc_w_g = m_bc_w_g.data();

  const Real* p_bc_ep_g = m_bc_ep_g.data();

  auto bcs_function = [p_bc_u_g,p_bc_v_g,p_bc_w_g,p_bc_ep_g]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& /*flags_arr*/,
                      const int dir)
  {
    GpuArray<const Real*,3> p_bc_vel_g = {p_bc_u_g, p_bc_v_g, p_bc_w_g};

    if(bct == BCList::pinf || bct == BCList::pout) {
      mf_arr(ijk,n) = mf_arr(dom_ijk,n);
    }
    else if (bct == BCList::minf) {
      mf_arr(ijk,n) = p_bc_ep_g[bcv] * p_bc_vel_g[dir][bcv];
    }
  };

  set_bcs(time, bcs_function, mf_in);
}


void
MFIXBoundaryConditions::set_tracer_bcs (Real time,
                                        MFIXFluidPhase& /*fluid*/,
                                        Vector< MultiFab* > const& trac_in)
{
  BL_PROFILE("MFIXBoundaryConditions::set_tracer_bcs()");

//  set_tracer_bc_values(time, fluid);
  Real* p_bc_trac = m_bc_tracer.data();

  auto bcs_function = [p_bc_trac]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& /*flags_arr*/,
                      const int /*dir*/)
  {
    if(bct == BCList::pinf || bct == BCList::pout) {
      mf_arr(ijk,n) = mf_arr(dom_ijk,n);
    }
    else if (bct == BCList::minf) {
      mf_arr(ijk,n) = p_bc_trac[bcv];
    }
  };

  set_bcs(time, bcs_function, trac_in);
}


void
MFIXBoundaryConditions::set_species_bcs (Real time,
                                         MFIXFluidPhase& fluid,
                                         Vector< MultiFab* > const& X_gk_in)
{
  BL_PROFILE("MFIXBoundaryConditions::set_species_bcs()");

  set_species_bc_values(time, fluid);
  Real** p_bc_X_gk = m_bc_X_gk_ptr.data();

  auto bcs_function = [p_bc_X_gk]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& /*flags_arr*/,
                      const int /*dir*/)
  {
    if(bct == BCList::pout) {
      mf_arr(ijk,n) = mf_arr(dom_ijk,n);
    }
    else if (bct == BCList::minf || bct == BCList::pinf) {
      mf_arr(ijk,n) = p_bc_X_gk[n][bcv];
    }
  };

  set_bcs(time, bcs_function, X_gk_in);
}


void
MFIXBoundaryConditions::set_density_bcs (Real time,
                                         Vector< MultiFab* > const& ro_g_in)
{
  BL_PROFILE("MFIXBoundaryConditions::set_density_bcs()");

  set_density_bc_values(time);
  Real* p_bc_ro_g = m_bc_ro_g.data();

  auto bcs_function = [p_bc_ro_g]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& /*flags_arr*/,
                      const int /*dir*/)
  {
    if(bct == BCList::pinf || bct == BCList::pout) {
      mf_arr(ijk,n) = mf_arr(dom_ijk,n);
    }
    else if (bct == BCList::minf) {
      mf_arr(ijk,n) = p_bc_ro_g[bcv];
    }
  };

  set_bcs(time, bcs_function, ro_g_in);
}
