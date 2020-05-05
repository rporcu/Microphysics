#include <mfix_calc_cp_g.H>
#include <eos_mod.H>
#include <param_mod_F.H>
#include <MFIX_FLUID_Parms.H>

void calc_cp_g (const Box& bx,
                FArrayBox& temp_fab,
                FArrayBox& cp_g_fab)
{

  const Real cp_g0 = FLUID::Cp_g0;

  Array4<Real> const& cp_g = cp_g_fab.array();

  amrex::ParallelFor(bx, [cp_g,cp_g0]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {cp_g(i,j,k) = cp_g0;});
}
