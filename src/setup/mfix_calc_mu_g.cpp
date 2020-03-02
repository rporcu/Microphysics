#include <mfix_calc_mu_g.H>
#include <eos_mod.H>
#include <param_mod_F.H>
#include <MFIX_FLUID_Parms.H>

void calc_mu_g (const Box& bx,
                FArrayBox& mu_g_fab)
{
  Real mu_val;

  const Real mu_g0 = FLUID::mu_g0;

  Array4<Real> const& mu_g = mu_g_fab.array();

  // Set the initial viscosity
  if(is_undefined_db_cpp(mu_g0))
    mu_val = sutherland(293.15);
  else
    mu_val = mu_g0;

  amrex::ParallelFor(bx, [mu_g,mu_val]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {mu_g(i,j,k) = mu_val;});
}
