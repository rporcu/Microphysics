#include <mfix_calc_cp_g.H>
#include <eos_mod.H>
#include <param_mod_F.H>
#include <MFIX_FLUID_Parms.H>

void calc_cp_g (const Box& bx,
                FArrayBox& temp_fab,
                FArrayBox& cp_g_fab)
{
  Real cp_val;

  const Real cp_g0 = FLUID::Cp_g0;

  Array4<Real> const& cp_g = cp_g_fab.array();

  // Set the initial Cp_g
  if (is_undefined_db_cpp(cp_g0))
    amrex::Abort("Right now we need to define Cp_g0");
  else
    cp_val = cp_g0;

  amrex::ParallelFor(bx, [cp_g,cp_val]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {cp_g(i,j,k) = cp_val;});
}
