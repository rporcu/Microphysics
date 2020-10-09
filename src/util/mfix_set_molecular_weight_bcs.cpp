#include <mfix.H>
#include <mfix_fluid_parms.H>

using namespace amrex;

void
mfix::set_molecular_weight_bcs (Real time,
                                const int lev,
                                FArrayBox& scal_fab,
                                const Box& domain)
{
  BL_PROFILE("mfix::set_molecular_weight_bcs()");

  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  Array4<Real> const& scal_arr = scal_fab.array();

  Real bc0 = FLUID::MW_g0;

  const int nspecies_g = FLUID::nspecies;

  Gpu::DeviceVector< Real > MW_gk0_d(nspecies_g);
  Gpu::copyAsync(Gpu::hostToDevice, FLUID::MW_gk0.begin(), FLUID::MW_gk0.end(), MW_gk0_d.begin());

  // Flag to understand if fluid is a mixture
  const int fluid_is_a_mixture = FLUID::is_a_mixture;

  Real** p_bc_X_gk = fluid_is_a_mixture ? m_bc_X_gk_ptr.data() : nullptr;
  Real* p_MW_gk0 = fluid_is_a_mixture ? MW_gk0_d.data() : nullptr;

  IntVect scal_lo(scal_fab.loVect());
  IntVect scal_hi(scal_fab.hiVect());

  const int nlft = amrex::max(0, dom_lo[0]-scal_lo[0]);
  const int nbot = amrex::max(0, dom_lo[1]-scal_lo[1]);
  const int ndwn = amrex::max(0, dom_lo[2]-scal_lo[2]);

  const int nrgt = amrex::max(0, scal_hi[0]-dom_hi[0]);
  const int ntop = amrex::max(0, scal_hi[1]-dom_hi[1]);
  const int nup  = amrex::max(0, scal_hi[2]-dom_hi[2]);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(scal_hi), bx_xz_lo_hi_3D(scal_hi), bx_xy_lo_hi_3D(scal_hi);
  IntVect bx_yz_hi_lo_3D(scal_lo), bx_xz_hi_lo_3D(scal_lo), bx_xy_hi_lo_3D(scal_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for GPU loops
  const Box bx_yz_lo_3D(scal_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, scal_hi);

  const Box bx_xz_lo_3D(scal_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, scal_hi);

  const Box bx_xy_lo_3D(scal_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, scal_hi);

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D,
      [bct_ilo,dom_lo,bc0,pinf,pout,minf,scal_arr,
       fluid_is_a_mixture,nspecies_g,p_MW_gk0,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if (bct == pout)
      {
        scal_arr(i,j,k) = scal_arr(dom_lo[0],j,k);
      }
      else if (bct == minf or bct == pinf)
      {
        if (not fluid_is_a_mixture) {
          scal_arr(i,j,k) = bc0;
        }
        else {
          Real MW_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_g_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          scal_arr(i,j,k) = 1./MW_g_sum;
        }
      }
    });
  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D,
      [bct_ihi,dom_hi,bc0,pinf,pout,minf,scal_arr,
       fluid_is_a_mixture,nspecies_g,p_MW_gk0,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if (bct == pout)
      {
        scal_arr(i,j,k) = scal_arr(dom_hi[0],j,k);
      }
      else if (bct == minf or bct == pinf)
      {
        if (not fluid_is_a_mixture) {
          scal_arr(i,j,k) = bc0;
        }
        else {
          Real MW_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_g_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          scal_arr(i,j,k) = 1./MW_g_sum;
        }
      }
    });
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D,
      [bct_jlo,dom_lo,bc0,pinf,pout,minf,scal_arr,
       fluid_is_a_mixture,nspecies_g,p_MW_gk0,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if (bct == pout)
      {
        scal_arr(i,j,k) = scal_arr(i,dom_lo[1],k);
      }
      else if (bct == minf or bct == pinf)
      {
        if (not fluid_is_a_mixture) {
          scal_arr(i,j,k) = bc0;
        }
        else {
          Real MW_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_g_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          scal_arr(i,j,k) = 1./MW_g_sum;
        }
      }
    });
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D,
      [bct_jhi,dom_hi,bc0,pinf,pout,minf,scal_arr,
       fluid_is_a_mixture,nspecies_g,p_MW_gk0,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if (bct == pout)
      {
        scal_arr(i,j,k) = scal_arr(i,dom_hi[1],k);
      }
      else if (bct == minf or bct == pinf)
      {
        if (not fluid_is_a_mixture) {
          scal_arr(i,j,k) = bc0;
        }
        else {
          Real MW_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_g_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          scal_arr(i,j,k) = 1./MW_g_sum;
        }
      }
    });
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D,
      [bct_klo,dom_lo,bc0,pinf,pout,minf,scal_arr,
       fluid_is_a_mixture,nspecies_g,p_MW_gk0,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if (bct == pout)
      {
        scal_arr(i,j,k) = scal_arr(i,j,dom_lo[2]);
      }
      else if (bct == minf or bct == pinf)
      {
        if (not fluid_is_a_mixture) {
          scal_arr(i,j,k) = bc0;
        }
        else {
          Real MW_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_g_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          scal_arr(i,j,k) = 1./MW_g_sum;
        }
      }
    });
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D,
      [bct_khi,dom_hi,bc0,pinf,pout,minf,scal_arr,
       fluid_is_a_mixture,nspecies_g,p_MW_gk0,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if (bct == pout)
      {
        scal_arr(i,j,k) = scal_arr(i,j,dom_hi[2]);
      }
      else if (bct == minf or bct == pinf)
      {
        if (not fluid_is_a_mixture) {
          scal_arr(i,j,k) = bc0;
        }
        else {
          Real MW_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            MW_g_sum += p_bc_X_gk[n][bcv]/p_MW_gk0[n];
          }

          scal_arr(i,j,k) = 1./MW_g_sum;
        }
      }
    });
  }

  Gpu::synchronize();
}
