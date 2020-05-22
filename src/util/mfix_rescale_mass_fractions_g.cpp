#include <mfix.H>
#include <mfix_rescale_mass_fractions_g.H>

#include <AMReX_VisMF.H>
#include <MFIX_MFHelpers.H>
#include <MFIX_FLUID_Parms.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

void rescale_species(const Vector< MultiFab* >& X_g)
{
  const int nspecies_g = X_g[0]->nComp();
  const int finest_level = X_g.size() - 1;

  for (int lev = 0; lev <= finest_level; lev++) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(not X_g[lev]->contains_nan(),
        "The species computation contains NaN");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(not X_g[lev]->contains_inf(),
        "The species computation contains Inf");
  }

  Vector< Real > X_g_max(finest_level+1, std::numeric_limits<Real>::min());
  Vector< Real > X_g_min(finest_level+1, std::numeric_limits<Real>::max());

// FOR DEBUG PURPOSES
////////////////////////////////////////////////////////
  for (int lev = 0; lev <= finest_level; lev++)
  {
    for (int n = 0; n < nspecies_g; ++n)
    {
      X_g_max[lev] = amrex::max(X_g_max[lev], X_g[lev]->max(n));
      X_g_min[lev] = amrex::min(X_g_min[lev], X_g[lev]->min(n));
    }
  }

  //const Real epsilon = std::numeric_limits<Real>::epsilon();

  for (int lev = 0; lev <= finest_level; lev++)
  {
    amrex::Print() << "On level " << lev << " :" << std::endl;
    amrex::Print() << "Max fluid species mass fraction = " << X_g_max[lev] << std::endl;
    amrex::Print() << "Min fluid species mass fraction = " << X_g_min[lev] << std::endl;

    // For now we take out the following sanity check
    //if (((X_g_max[lev]-1) > epsilon) or (X_g_min[lev] < -epsilon)) {
    //  amrex::Abort("Fluid mass fractions exceed admissibile domain [0,1] of definition");
    //}
  }
////////////////////////////////////////////////////////
// END DEBUG

// This is a cheat: we try to make this rescaling more numerically robust
// if X_g(i,j,k,n) > 1 --> X_g(i,j,k,n) = 1;
// if X_g(i,j,k,n) < 0 --> X_g(i,j,k,n) = 0
  for (int lev = 0; lev <= finest_level; lev++)
  {
    for (MFIter mfi(*X_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      // Tilebox
      Box bx = mfi.tilebox();

      const auto& X_g_array = X_g[lev]->array(mfi);

      amrex::ParallelFor(bx, nspecies_g,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        if (X_g_array(i,j,k,n) > 1)
          X_g_array(i,j,k,n) = 1;
        
        if (X_g_array(i,j,k,n) < 0)
          X_g_array(i,j,k,n) = 0;
      });
    }
  }

  // The following MultiFabs store the cellwise sum of species mass fractions
  Vector< MultiFab* > X_g_sum(finest_level+1, nullptr);

  // We allocate memory
  for (int lev = 0; lev <= finest_level; lev++)
  {
    X_g_sum[lev] = new MultiFab (X_g[lev]->boxArray(),
        X_g[lev]->DistributionMap(), 1, 0, MFInfo(), X_g[lev]->Factory());

    X_g_sum[lev]->setVal(0.);
  }

  // Cellwise sum fluid species mass fractions on every level
  for (int lev = 0; lev <= finest_level; lev++)
  {
    for (int n = 0; n < nspecies_g; n++) {
      MultiFab::Add(*X_g_sum[lev], *X_g[lev], n, 0, 1, 0);
    }
  }

// FOR DEBUG PURPOSES
///////////////////////////////////////////////////////////
  {
    Vector< Real > X_g_sum_max(finest_level+1, std::numeric_limits<Real>::min());
    Vector< Real > X_g_sum_min(finest_level+1, std::numeric_limits<Real>::max());
    
    for (int lev = 0; lev <= finest_level; lev++)
    {
      X_g_sum_max[lev] = amrex::max(X_g_sum_max[lev], X_g_sum[lev]->max(0));
      X_g_sum_min[lev] = amrex::min(X_g_sum_min[lev], X_g_sum[lev]->min(0));
    }

    // For now we take out the following sanity check
    //const Real min_X_g_sum = 0.5; // TODO: check which value we want in here
    //const Real max_X_g_sum = 1.5; // TODO: check which value we want in here
    
    for (int lev = 0; lev <= finest_level; lev++)
    {
      amrex::Print() << "On level " << lev << " :" << std::endl;
      amrex::Print() << "Max SUMMED fluid species mass fractions BEFORE rescale = " << X_g_sum_max[lev] << std::endl;
      amrex::Print() << "Min SUMMED fluid species mass fractions BEFORE rescale = " << X_g_sum_min[lev] << std::endl;

      // For now we take out the following sanity check
      //if ((X_g_sum_max[lev] > max_X_g_sum) or (X_g_sum_min[lev] < min_X_g_sum)) {
      //  amrex::Abort("Fluid summed mass fractions exceed admissibile values");
      //}
    }
  }
///////////////////////////////////////////////////////////
// END DEBUG

  // Rescale species mass fractions in order to have their cellwise sum = 1
  for (int lev = 0; lev <= finest_level; lev++)
  {
    for (MFIter mfi(*X_g_sum[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      const auto& X_g_arr     = X_g[lev]->array(mfi);
      const auto& X_g_sum_arr = X_g_sum[lev]->array(mfi);

      const EBFArrayBox& X_g_fab =
        static_cast<EBFArrayBox const&>((*X_g[lev])[mfi]);

      const EBCellFlagFab& flags_fab = X_g_fab.getEBCellFlagFab();

      const auto& flags = flags_fab.array();

      amrex::ParallelFor(bx, nspecies_g,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        if (not flags(i,j,k).isCovered())
          X_g_arr(i,j,k,n) /= X_g_sum_arr(i,j,k);
      });
    }
  }

// FOR DEBUG PURPOSES
///////////////////////////////////////////////////////////
  // We recompute the sum
  for (int lev = 0; lev <= finest_level; lev++)
  {
    X_g_sum[lev]->setVal(0.);

    for (int n = 0; n < nspecies_g; n++) {
      MultiFab::Add(*X_g_sum[lev], *X_g[lev], n, 0, 1, 0);
    }
  }

  {
    Vector< Real > X_g_sum_max(finest_level+1, std::numeric_limits<Real>::min());
    Vector< Real > X_g_sum_min(finest_level+1, std::numeric_limits<Real>::max());
    
    for (int lev = 0; lev <= finest_level; lev++)
    {
      X_g_sum_max[lev] = amrex::max(X_g_sum_max[lev], X_g_sum[lev]->max(0));
      X_g_sum_min[lev] = amrex::min(X_g_sum_min[lev], X_g_sum[lev]->min(0));
    }

    //const Real epsilon = std::numeric_limits<Real>::epsilon();
    const Real epsilon = 1.e-12;

    for (int lev = 0; lev <= finest_level; lev++)
    {
      amrex::Print() << "On level " << lev << " :" << std::endl;
      amrex::Print() << "Max SUMMED fluid species mass fractions AFTER rescale = " << X_g_sum_max[lev] << std::endl;
      amrex::Print() << "Min SUMMED fluid species mass fractions AFTER rescale = " << X_g_sum_min[lev] << std::endl;

      if (std::fabs(X_g_sum_max[lev]-1) > epsilon or
          std::fabs(X_g_sum_min[lev]-1) > epsilon)
      {
        amrex::Abort("Fluid mass fractions rescaling FAILED");
      }
    }
  }
///////////////////////////////////////////////////////////
// END DEBUG

  // Deallocate memory
  for (int lev = 0; lev <= finest_level; lev++)
  {
    delete X_g_sum[lev];
  }

}
