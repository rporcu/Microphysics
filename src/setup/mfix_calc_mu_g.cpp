#include <mfix_calc_mu_g.hpp>
#include <fld_constants_mod_F.H>
#include <eos_mod.hpp>
#include <param_mod_F.H>

void calc_mu_g(const Box& bx,
               FArrayBox& mu_g_fab)
{
  Real mu_val;

  const Real mu_g0 = get_mu_g0();

  Array4<Real> const& mu_g = mu_g_fab.array();

  // Set the initial viscosity
  if(is_undefined_db(mu_g0))
    mu_val = sutherland(293.15);
  else
    mu_val = mu_g0;

  AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k, {mu_g(i,j,k) = mu_val;});

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif
}
