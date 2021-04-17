#include <mfix.H>

#include <AMReX_VisMF.H>
#include <mfix_mf_helpers.H>
#include <mfix_fluid_parms.H>
#include <mfix_calc_fluid_coeffs.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif


void
mfix::mfix_normalize_fluid_species(const Vector< MultiFab* >& X_gk)
{
  const int nspecies_g = FLUID::nspecies;

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
mfix::mfix_update_fluid_and_species(const Vector< MultiFab* >& cp_gk,
                                     const Vector< MultiFab* >& h_gk,
                                     const Vector< MultiFab* >& MW_g,
                                     const Vector< MultiFab* >& cp_g,
                                     const Vector< MultiFab* >& h_g,
                                     const Vector< MultiFab* >& T_g,
                                     const Vector< MultiFab* >& X_gk)
{
  if (advect_enthalpy)
  {
    for (int lev(0); lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for(MFIter mfi(*MW_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();

        FArrayBox& cp_gk_fab = (*cp_gk[lev])[mfi];
        FArrayBox& h_gk_fab  = (*h_gk[lev])[mfi];
        FArrayBox& T_g_fab   = (*T_g[lev])[mfi];

        // Update species specific heat
        calc_cp_gk(bx, cp_gk_fab, T_g_fab);

        // Update species enthalpy
        calc_h_gk(bx, h_gk_fab, cp_gk_fab, T_g_fab);
      }
    }
  }

  const int nspecies_g = FLUID::nspecies;

  if (FLUID::is_a_mixture)
  {
    Gpu::DeviceVector< Real > MW_gk_d(nspecies_g);
    Gpu::copyAsync(Gpu::hostToDevice, FLUID::MW_gk0.begin(), FLUID::MW_gk0.end(), MW_gk_d.begin());
    Gpu::synchronize();

    Real* p_MW_gk = MW_gk_d.data();

    // Set covered values
    for (int lev(0); lev <= finest_level; lev++) {
      EB_set_covered(*X_gk[lev], 0, nspecies_g, 0, covered_val);
    }

    for (int lev(0); lev <= finest_level; lev++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for(MFIter mfi(*MW_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();

        Array4< Real > const& MW_g_arr = MW_g[lev]->array(mfi);
        Array4< Real > const& X_gk_arr = X_gk[lev]->array(mfi);

        ParallelFor(bx, [nspecies_g,X_gk_arr,p_MW_gk,MW_g_arr]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real MW_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_g_sum += X_gk_arr(i,j,k,n)*p_MW_gk[n];
          }

          MW_g_arr(i,j,k) = 1./MW_g_sum;
        });

        if (advect_enthalpy)
        {
          Array4< Real > const& cp_gk_arr = cp_gk[lev]->array(mfi);
          Array4< Real > const& h_gk_arr  = h_gk[lev]->array(mfi);
          Array4< Real > const& cp_g_arr  = cp_g[lev]->array(mfi);
          Array4< Real > const& h_g_arr   = h_g[lev]->array(mfi);

          ParallelFor(bx, [cp_gk_arr,h_gk_arr,cp_g_arr,h_g_arr,nspecies_g,
              X_gk_arr,MW_g_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            Real cp_g_sum(0);
            Real h_g_sum(0);

            for (int n(0); n < nspecies_g; n++) {
              cp_g_sum += X_gk_arr(i,j,k,n)*cp_gk_arr(i,j,k,n);
              h_g_sum += X_gk_arr(i,j,k,n)*h_gk_arr(i,j,k,n);
            }

            cp_g_arr(i,j,k) = cp_g_sum;
            h_g_arr(i,j,k) = h_g_sum;
          });
        }
      }
    }
  }

  // Set covered values
  for (int lev(0); lev <= finest_level; lev++) {
    if (FLUID::is_a_mixture)
      EB_set_covered(*MW_g[lev], 0, 1, 0, covered_val);

    if (advect_enthalpy) {
      EB_set_covered(*cp_g[lev], 0, 1, 0, covered_val);
      EB_set_covered(*h_g[lev], 0, 1, 0, covered_val);
      EB_set_covered(*cp_gk[lev], 0, nspecies_g, 0, covered_val);
      EB_set_covered(*h_gk[lev], 0, nspecies_g, 0, covered_val);
    }
  }

  // Fill boundary values
  for (int lev(0); lev <= finest_level; lev++) {
    if (FLUID::is_a_mixture)
      MW_g[lev]->FillBoundary(geom[lev].periodicity());

    if (advect_enthalpy) {
      cp_g[lev]->FillBoundary(geom[lev].periodicity());
      h_g[lev]->FillBoundary(geom[lev].periodicity());
      cp_gk[lev]->FillBoundary(geom[lev].periodicity());
      h_gk[lev]->FillBoundary(geom[lev].periodicity());
    }
  }
}
