#include <mfix_calc_mu_g.H>
#include <MFIX_FLUID_Parms.H>

void calc_mu_g (const Box& bx,
                FArrayBox& mu_g_fab)
{

  const Real mu_g0 = FLUID::mu_g0;

  Array4<Real> const& mu_g = mu_g_fab.array();

  amrex::ParallelFor(bx, [mu_g,mu_g0]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {mu_g(i,j,k) = mu_g0;});
}
