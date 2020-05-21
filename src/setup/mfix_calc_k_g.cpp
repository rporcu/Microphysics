#include <mfix_calc_k_g.H>
#include <MFIX_FLUID_Parms.H>

void calc_k_g (const Box& bx,
               FArrayBox& temp_fab,
               FArrayBox& k_g_fab)
{

  const Real k_g0 = FLUID::k_g0;

  Array4<Real> const& k_g = k_g_fab.array();

  amrex::ParallelFor(bx, [k_g,k_g0]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {k_g(i,j,k) = k_g0;});
}
