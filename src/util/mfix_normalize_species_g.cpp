#include <mfix.H>

#include <AMReX_VisMF.H>
#include <mfix_mf_helpers.H>
#include <mfix_fluid_parms.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif


void
mfix::mfix_normalize_fluid_species(const Vector< MultiFab* >& X_gk)
{
  const int nspecies_g = fluid.nspecies;

  for (int lev = 0; lev <= finest_level; lev++) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!X_gk[lev]->contains_nan(),
        "The species computation contains NaN");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!X_gk[lev]->contains_inf(),
        "The species computation contains Inf");
  }

#if defined(AMREX_DEBUG)
  {
    Vector< Real > X_gk_max(finest_level+1, std::numeric_limits<Real>::min());
    Vector< Real > X_gk_min(finest_level+1, std::numeric_limits<Real>::max());

    for (int lev = 0; lev <= finest_level; lev++)
    {
      for (int n = 0; n < nspecies_g; ++n)
      {
        X_gk_max[lev] = amrex::max(X_gk_max[lev], X_gk[lev]->max(n));
        X_gk_min[lev] = amrex::min(X_gk_min[lev], X_gk[lev]->min(n));
      }
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
      amrex::Print() << "On level " << lev << " :" << std::endl;
      amrex::Print() << "Max fluid species mass fraction = " << X_gk_max[lev] << std::endl;
      amrex::Print() << "Min fluid species mass fraction = " << X_gk_min[lev] << std::endl;
    }
  }
#endif

  for (int lev(0); lev <= finest_level; lev++)
  {
    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(X_gk[lev]->Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*X_gk[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      Array4<Real> const& X_gk_arr = X_gk[lev]->array(mfi);
      auto const& flags_arr = flags.const_array(mfi);

      amrex::ParallelFor(bx, [flags_arr,X_gk_arr,nspecies_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        for (int n(0); n < nspecies_g; ++n) {
          Real value = X_gk_arr(i,j,k,n);
          X_gk_arr(i,j,k,n) = value > 1 ? 1 : (value < 0 ? 0 : value);
        }

        if (!flags_arr(i,j,k).isCovered())
        {
          Real sum(0);

          for (int n(0); n < nspecies_g; ++n)
            sum += X_gk_arr(i,j,k,n);

          for (int n(0); n < nspecies_g; ++n)
            X_gk_arr(i,j,k,n) /= sum;
        }
      });
    }
  }

  // Update ghost cells
  for (int lev = 0; lev <= finest_level; lev++) {
    X_gk[lev]->FillBoundary(geom[lev].periodicity());
  }

#if defined(AMREX_DEBUG)
  {
    Vector< Real > X_gk_max(finest_level+1, std::numeric_limits<Real>::min());
    Vector< Real > X_gk_min(finest_level+1, std::numeric_limits<Real>::max());

    for (int lev = 0; lev <= finest_level; lev++)
    {
      for (int n = 0; n < nspecies_g; ++n)
      {
        X_gk_max[lev] = amrex::max(X_gk_max[lev], X_gk[lev]->max(n));
        X_gk_min[lev] = amrex::min(X_gk_min[lev], X_gk[lev]->min(n));
      }
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
      amrex::Print() << "On level " << lev << " :" << std::endl;
      amrex::Print() << "Max fluid species mass fraction AFTER normalization = " << X_gk_max[lev] << std::endl;
      amrex::Print() << "Min fluid species mass fraction AFTER normalization = " << X_gk_min[lev] << std::endl;

      if (X_gk_max[lev] > 1 || X_gk_min[lev] < 0)
        amrex::Abort("Fluid mass fractions rescaling FAILED");
    }
  }
#endif

}


void
mfix::mfix_update_fluid_and_species(const Vector< MultiFab* >& h_g,
                                    const Vector< MultiFab* >& T_g,
                                    const Vector< MultiFab* >& X_gk)
{
  const int nspecies_g = fluid.nspecies;

  auto& fluid_parms = *fluid.parameters;

  if (fluid.is_a_mixture)
  {
    // Set covered values
    for (int lev(0); lev <= finest_level; lev++) {
      EB_set_covered(*X_gk[lev], 0, nspecies_g, 0, covered_val);
    }

    for (int lev(0); lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for(MFIter mfi(*T_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();

        if (advect_enthalpy)
        {
          Array4< Real       > const& h_g_arr  = h_g[lev]->array(mfi);
          Array4< Real const > const& X_gk_arr = X_gk[lev]->array(mfi);
          Array4< Real const > const& T_g_arr  = T_g[lev]->array(mfi);

          ParallelFor(bx, [h_g_arr,nspecies_g,X_gk_arr,T_g_arr,fluid_parms]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            Real h_g_sum(0);

            const Real Tg_loc = T_g_arr(i,j,k);

            for (int n(0); n < nspecies_g; n++) {
              h_g_sum += X_gk_arr(i,j,k,n) * fluid_parms.calc_h_gk(Tg_loc,n);
            }

            h_g_arr(i,j,k) = h_g_sum;
          });
        }
      }
    }
  }

  // Set covered values
  for (int lev(0); lev <= finest_level; lev++) {
    if (advect_enthalpy) {
      EB_set_covered(*h_g[lev], 0, 1, 0, covered_val);
    }
  }

  // Fill boundary values
  for (int lev(0); lev <= finest_level; lev++) {
    if (advect_enthalpy) {
      h_g[lev]->FillBoundary(geom[lev].periodicity());
    }
  }
}
