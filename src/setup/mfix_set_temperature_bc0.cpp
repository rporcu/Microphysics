#include <mfix.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_mf_helpers.H>

void
mfix::set_temperature_bc0 (const Box& sbx,
                           MFIter* mfi,
                           const int lev,
                           const Box& domain)
{
  const Real cp_g0  = FLUID::cp_g0;
  const Real k_g0   = FLUID::k_g0;

  amrex::Real* p_bc_t_g = m_bc_t_g.data();
 
  const int nspecies_g = FLUID::nspecies;

  Gpu::DeviceVector< Real > cp_gk0_d(nspecies_g);
  Gpu::copyAsync(Gpu::hostToDevice, FLUID::cp_gk0.begin(), FLUID::cp_gk0.end(), cp_gk0_d.begin());

  // Flag to understand if fluid is a mixture
  const int fluid_is_a_mixture = FLUID::is_a_mixture;

  Real** p_bc_X_gk = fluid_is_a_mixture ? m_bc_X_gk_ptr.data() : nullptr;
  Real* p_cp_gk0 = fluid_is_a_mixture ? cp_gk0_d.data() : nullptr;

  Array4<Real> const& a_T_g  = m_leveldata[lev]->T_g->array(*mfi);
  Array4<Real> const& a_cp_g = m_leveldata[lev]->cp_g->array(*mfi);
  Array4<Real> const& a_k_g  = m_leveldata[lev]->k_g->array(*mfi);
  Array4<Real> const& a_h_g  = m_leveldata[lev]->h_g->array(*mfi);

  const IntVect sbx_lo(sbx.loVect());
  const IntVect sbx_hi(sbx.hiVect());

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  Array4<const int> const& a_bc_ilo = bc_ilo[lev]->array();
  Array4<const int> const& a_bc_ihi = bc_ihi[lev]->array();
  Array4<const int> const& a_bc_jlo = bc_jlo[lev]->array();
  Array4<const int> const& a_bc_jhi = bc_jhi[lev]->array();
  Array4<const int> const& a_bc_klo = bc_klo[lev]->array();
  Array4<const int> const& a_bc_khi = bc_khi[lev]->array();

  const int nlft = amrex::max(0,dom_lo[0]-sbx_lo[0]);
  const int nbot = amrex::max(0,dom_lo[1]-sbx_lo[1]);
  const int ndwn = amrex::max(0,dom_lo[2]-sbx_lo[2]);

  const int nrgt = amrex::max(0,sbx_hi[0]-dom_hi[0]);
  const int ntop = amrex::max(0,sbx_hi[1]-dom_hi[1]);
  const int nup  = amrex::max(0,sbx_hi[2]-dom_hi[2]);

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  if (nlft > 0)
  {
    IntVect bx_yz_lo_hi_3D(sbx_hi);

    // Fix lo and hi limits
    bx_yz_lo_hi_3D[0] = dom_lo[0]-1;

    const Box bx_yz_lo_3D(sbx_lo, bx_yz_lo_hi_3D);

    ParallelFor(bx_yz_lo_3D, [a_bc_ilo,dom_lo,pinf,pout,minf,a_T_g,a_k_g,p_bc_t_g,
        k_g0,a_cp_g,a_h_g,cp_g0,fluid_is_a_mixture,nspecies_g,p_bc_X_gk,
        p_cp_gk0]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_ilo(dom_lo[0]-1,j,k,1);
      const int bct = a_bc_ilo(dom_lo[0]-1,j,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_k_g(i,j,k)  = k_g0;

        if (not fluid_is_a_mixture) {
          a_cp_g(i,j,k) = cp_g0;
          a_h_g(i,j,k)  = cp_g0*p_bc_t_g[bcv];
        }
        else {
          Real cp_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            cp_g_sum += p_bc_X_gk[n][bcv]*p_cp_gk0[n];
          }

          a_cp_g(i,j,k) = cp_g_sum;
          a_h_g(i,j,k)  = cp_g_sum * p_bc_t_g[bcv];
        }
      }
    });
  }

  if (nrgt > 0)
  {
    IntVect bx_yz_hi_lo_3D(sbx_lo);

    // Fix lo and hi limits
    bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

    const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, sbx_hi);

    ParallelFor(bx_yz_hi_3D, [a_bc_ihi,dom_hi,pinf,pout,minf,a_T_g,a_k_g,p_bc_t_g,
        k_g0,a_cp_g,a_h_g,cp_g0,fluid_is_a_mixture,nspecies_g,p_bc_X_gk,
        p_cp_gk0]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_ihi(dom_hi[0]+1,j,k,1);
      const int bct = a_bc_ihi(dom_hi[0]+1,j,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_k_g(i,j,k)  = k_g0;

        if (not fluid_is_a_mixture) {
          a_cp_g(i,j,k) = cp_g0;
          a_h_g(i,j,k)  = cp_g0*p_bc_t_g[bcv];
        }
        else {
          Real cp_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            cp_g_sum += p_bc_X_gk[n][bcv]*p_cp_gk0[n];
          }

          a_cp_g(i,j,k) = cp_g_sum;
          a_h_g(i,j,k)  = cp_g_sum * p_bc_t_g[bcv];
        }
      }
    });
  }

  if (nbot > 0)
  {
    IntVect bx_xz_lo_hi_3D(sbx_hi);

    // Fix lo and hi limits
    bx_xz_lo_hi_3D[1] = dom_lo[1]-1;

    const Box bx_xz_lo_3D(sbx_lo, bx_xz_lo_hi_3D);

    ParallelFor(bx_xz_lo_3D, [a_bc_jlo,dom_lo,pinf,pout,minf,a_T_g,a_k_g,p_bc_t_g,
        k_g0,a_cp_g,a_h_g,cp_g0,fluid_is_a_mixture,nspecies_g,p_bc_X_gk,
        p_cp_gk0]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_jlo(i,dom_lo[1]-1,k,1);
      const int bct = a_bc_jlo(i,dom_lo[1]-1,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_k_g(i,j,k)  = k_g0;

        if (not fluid_is_a_mixture) {
          a_cp_g(i,j,k) = cp_g0;
          a_h_g(i,j,k)  = cp_g0*p_bc_t_g[bcv];
        }
        else {
          Real cp_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            cp_g_sum += p_bc_X_gk[n][bcv]*p_cp_gk0[n];
          }

          a_cp_g(i,j,k) = cp_g_sum;
          a_h_g(i,j,k)  = cp_g_sum * p_bc_t_g[bcv];
        }
      }
    });
  }

  if (ntop > 0)
  {
    IntVect bx_xz_hi_lo_3D(sbx_lo);

    // Fix lo and hi limits
    bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

    const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, sbx_hi);

    ParallelFor(bx_xz_hi_3D, [a_bc_jhi,dom_hi,pinf,pout,minf,a_T_g,a_k_g,p_bc_t_g,
        k_g0,a_cp_g,a_h_g,cp_g0,fluid_is_a_mixture,nspecies_g,p_bc_X_gk,
        p_cp_gk0]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_jhi(i,dom_hi[1]+1,k,1);
      const int bct = a_bc_jhi(i,dom_hi[1]+1,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_k_g(i,j,k)  = k_g0;

        if (not fluid_is_a_mixture) {
          a_cp_g(i,j,k) = cp_g0;
          a_h_g(i,j,k)  = cp_g0*p_bc_t_g[bcv];
        }
        else {
          Real cp_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            cp_g_sum += p_bc_X_gk[n][bcv]*p_cp_gk0[n];
          }

          a_cp_g(i,j,k) = cp_g_sum;
          a_h_g(i,j,k)  = cp_g_sum * p_bc_t_g[bcv];
        }
      }
    });
  }

  if (ndwn > 0)
  {
    IntVect bx_xy_lo_hi_3D(sbx_hi);

    // Fix lo and hi limits
    bx_xy_lo_hi_3D[2] = dom_lo[2]-1;

    const Box bx_xy_lo_3D(sbx_lo, bx_xy_lo_hi_3D);

    ParallelFor(bx_xy_lo_3D, [a_bc_klo,dom_lo,pinf,pout,minf,a_T_g,a_k_g,p_bc_t_g,
        k_g0,a_cp_g,a_h_g,cp_g0,fluid_is_a_mixture,nspecies_g,p_bc_X_gk,
        p_cp_gk0]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_klo(i,j,dom_lo[2]-1,1);
      const int bct = a_bc_klo(i,j,dom_lo[2]-1,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_k_g(i,j,k)  = k_g0;

        if (not fluid_is_a_mixture) {
          a_cp_g(i,j,k) = cp_g0;
          a_h_g(i,j,k)  = cp_g0*p_bc_t_g[bcv];
        }
        else {
          Real cp_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            cp_g_sum += p_bc_X_gk[n][bcv]*p_cp_gk0[n];
          }

          a_cp_g(i,j,k) = cp_g_sum;
          a_h_g(i,j,k)  = cp_g_sum * p_bc_t_g[bcv];
        }
      }
    });
  }

  if (nup > 0)
  {
    IntVect bx_xy_hi_lo_3D(sbx_lo);

    // Fix lo and hi limits
    bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

    const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, sbx_hi);

    ParallelFor(bx_xy_hi_3D, [a_bc_khi,dom_hi,pinf,pout,minf,a_T_g,a_k_g,p_bc_t_g,
        k_g0,a_cp_g,a_h_g,cp_g0,fluid_is_a_mixture,nspecies_g,p_bc_X_gk,
        p_cp_gk0]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_khi(i,j,dom_hi[2]+1,1);
      const int bct = a_bc_khi(i,j,dom_hi[2]+1,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        a_T_g(i,j,k)  = p_bc_t_g[bcv];
        a_k_g(i,j,k)  = k_g0;

        if (not fluid_is_a_mixture) {
          a_cp_g(i,j,k) = cp_g0;
          a_h_g(i,j,k)  = cp_g0*p_bc_t_g[bcv];
        }
        else {
          Real cp_g_sum(0);

          for (int n(0); n < nspecies_g; n++) {
            cp_g_sum += p_bc_X_gk[n][bcv]*p_cp_gk0[n];
          }

          a_cp_g(i,j,k) = cp_g_sum;
          a_h_g(i,j,k)  = cp_g_sum * p_bc_t_g[bcv];
        }
      }
    });
  }

  Gpu::synchronize();

}
