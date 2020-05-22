#include <mfix_calc_fluid_coeffs.H>
#include <eos_mod.H>
#include <MFIX_FLUID_Parms.H>

void calc_mu_g (const Box& bx,
                FArrayBox& /*T_g_fab*/,
                FArrayBox& mu_g_fab)
{
  const Real mu_g0 = FLUID::mu_g0;

  Array4<Real> const& mu_g = mu_g_fab.array();

  amrex::ParallelFor(bx, [mu_g,mu_g0]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {mu_g(i,j,k) = mu_g0;});
}

void calc_cp_g (const Box& bx,
                FArrayBox& /*T_g_fab*/,
                FArrayBox& cp_g_fab)
{

  const Real cp_g0 = FLUID::cp_g0;

  Array4<Real> const& cp_g = cp_g_fab.array();

  amrex::ParallelFor(bx, [cp_g,cp_g0]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {cp_g(i,j,k) = cp_g0;});
}

void calc_k_g (const Box& bx,
               FArrayBox& /*T_g_fab*/,
               FArrayBox& k_g_fab)
{

  const Real k_g0 = FLUID::k_g0;

  Array4<Real> const& k_g = k_g_fab.array();

  amrex::ParallelFor(bx, [k_g,k_g0]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {k_g(i,j,k) = k_g0;});
}

void calc_D_g (const Box& bx,
               FArrayBox& /*T_g_fab*/,
               FArrayBox& D_g_fab)
{
  Real* D_g0;
  D_g0 = new Real [FLUID::nspecies_g];
  for (int n(0); n < FLUID::nspecies_g; n++)
    D_g0[n] = FLUID::D_g0[n];

  Array4<Real> const& D_g = D_g_fab.array();

  amrex::ParallelFor(bx, FLUID::nspecies_g, [D_g, D_g0]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {D_g(i,j,k,n) = D_g0[n];});

  delete D_g0;
}
