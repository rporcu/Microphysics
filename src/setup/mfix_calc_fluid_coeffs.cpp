#include <mfix_calc_fluid_coeffs.H>
#include <mfix_eos_mod.H>
#include <mfix_fluid_parms.H>


void calc_mu_g (const Box& bx,
                FArrayBox& mu_g_fab)
{
  const Real mu_g0 = FLUID::mu_g0;

  Array4<Real> const& mu_g = mu_g_fab.array();

  amrex::ParallelFor(bx, [mu_g,mu_g0]
  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  { mu_g(i,j,k) = mu_g0; });
}


void calc_cp_g (const Box& bx,
                FArrayBox& cp_g_fab,
                FArrayBox& /*T_g_fab*/)
{
  const Real cp_g0 = FLUID::cp_g0;

  Array4<Real> const& cp_g = cp_g_fab.array();

  amrex::ParallelFor(bx, [cp_g,cp_g0]
  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  { cp_g(i,j,k) = cp_g0; });
}


void calc_h_g (const Box& bx,
               FArrayBox& h_g_fab,
               FArrayBox& cp_g_fab,
               FArrayBox& T_g_fab)
{
  Array4<Real> const& h_g  = h_g_fab.array();
  Array4<Real> const& cp_g = cp_g_fab.array();
  Array4<Real> const& T_g  = T_g_fab.array();

  amrex::ParallelFor(bx, [h_g,cp_g,T_g]
  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  { h_g(i,j,k) = cp_g(i,j,k)*T_g(i,j,k); });
}


void calc_k_g (const Box& bx,
               FArrayBox& k_g_fab,
               FArrayBox& /*T_g_fab*/)
{
  const Real k_g0 = FLUID::k_g0;

  Array4<Real> const& k_g = k_g_fab.array();

  amrex::ParallelFor(bx, [k_g,k_g0]
  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  { k_g(i,j,k) = k_g0; });
}


void calc_D_gk (const Box& bx,
                FArrayBox& D_gk_fab)
{
  const int nspecies_g = FLUID::nspecies;
  Gpu::DeviceVector< Real> D_gk0(nspecies_g);
  Gpu::copyAsync(Gpu::hostToDevice, FLUID::D_gk0.begin(), FLUID::D_gk0.end(), D_gk0.begin());

  Real* p_D_gk0 = D_gk0.data();

  Array4<Real> const& D_gk = D_gk_fab.array();

  amrex::ParallelFor(bx, nspecies_g, [D_gk, p_D_gk0]
  AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  { D_gk(i,j,k,n) = p_D_gk0[n]; });

  Gpu::synchronize();
}


void calc_cp_gk (const Box& bx,
                 FArrayBox& cp_gk_fab,
                 FArrayBox& /*T_g_fab*/)
{
  const int nspecies_g = FLUID::nspecies;
  Gpu::DeviceVector< Real> cp_gk0(nspecies_g);
  Gpu::copyAsync(Gpu::hostToDevice, FLUID::cp_gk0.begin(), FLUID::cp_gk0.end(), cp_gk0.begin());

  Real* p_cp_gk0 = cp_gk0.data();

  Array4<Real> const& cp_gk = cp_gk_fab.array();

  amrex::ParallelFor(bx, nspecies_g, [cp_gk, p_cp_gk0]
  AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  { cp_gk(i,j,k,n) = p_cp_gk0[n]; });

  Gpu::synchronize();
}


void calc_h_gk (const Box& bx,
                FArrayBox& h_gk_fab,
                FArrayBox& cp_gk_fab,
                FArrayBox& T_g_fab)
{
  const int nspecies_g = FLUID::nspecies;

  Array4<Real> const& h_gk = h_gk_fab.array();
  Array4<Real> const& cp_gk = cp_gk_fab.array();
  Array4<Real> const& T_g = T_g_fab.array();

  amrex::ParallelFor(bx, nspecies_g, [h_gk,cp_gk,T_g]
  AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  { h_gk(i,j,k,n) = cp_gk(i,j,k,n)*T_g(i,j,k); });
}
