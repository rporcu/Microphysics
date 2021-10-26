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

      Array4<Real> dummy_array4;

      Array4<Real> const& X_gk_arr = X_gk[lev]->array(mfi);

      auto const& flags_arr = flags.const_array(mfi);

      auto& fluid_parms = *fluid.parameters;
      const int adv_enthalpy = advect_enthalpy;

      amrex::ParallelFor(bx, [flags_arr,X_gk_arr,nspecies_g,
          adv_enthalpy,fluid_parms]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        for (int n(0); n < nspecies_g; ++n) {
          Real value = X_gk_arr(i,j,k,n);
          X_gk_arr(i,j,k,n) = value > 1 ? 1 : (value < 0 ? 0 : value);
        }

        if (!flags_arr(i,j,k).isCovered()) {
          Real X_gk_sum(0);

          for (int n(0); n < nspecies_g; ++n)
            X_gk_sum += X_gk_arr(i,j,k,n);

          for (int n(0); n < nspecies_g; ++n)
            X_gk_arr(i,j,k,n) /= X_gk_sum;
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
