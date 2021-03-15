#include <mfix.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_mf_helpers.H>

void
mfix::set_species_bc0 (const Box& sbx,
                       MFIter* mfi,
                       const int lev,
                       const Box& domain)
{
  const int nspecies_g = FLUID::nspecies;

  Gpu::DeviceVector< Real > D_gk0_d(nspecies_g);
  Gpu::copyAsync(Gpu::hostToDevice, FLUID::D_gk0.begin(), FLUID::D_gk0.end(), D_gk0_d.begin());

  Real** p_bc_X_gk = m_bc_X_gk_ptr.data();
  Real* p_D_gk0 = D_gk0_d.data();

  // In case of advect_enthalpy
  Gpu::DeviceVector< Real > cp_gk0_d;
  if (advect_enthalpy) {
      cp_gk0_d.resize(nspecies_g);
      Gpu::copyAsync(Gpu::hostToDevice, FLUID::cp_gk0.begin(), FLUID::cp_gk0.end(), cp_gk0_d.begin());
  }

  const int loc_advect_enthalpy = advect_enthalpy;

  Real* p_cp_gk0 = loc_advect_enthalpy ? cp_gk0_d.data() : nullptr;
  Real* p_bc_t_g = loc_advect_enthalpy ? m_bc_t_g.data() : nullptr;

  // Get data
  Array4<Real> const& a_X_gk = m_leveldata[lev]->X_gk->array(*mfi);
  Array4<Real> const& a_D_gk = m_leveldata[lev]->D_gk->array(*mfi);

  Array4<Real> const& a_cp_gk = loc_advect_enthalpy ?
    m_leveldata[lev]->cp_gk->array(*mfi) : Array4<Real>();

  Array4<Real> const& a_h_gk = loc_advect_enthalpy ?
    m_leveldata[lev]->h_gk->array(*mfi) : Array4<Real>();

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

    ParallelFor(bx_yz_lo_3D, [a_bc_ilo,dom_lo,pinf,pout,minf,nspecies_g,a_X_gk,
        a_D_gk,a_cp_gk,a_h_gk,p_bc_X_gk,p_D_gk0,p_cp_gk0,p_bc_t_g,loc_advect_enthalpy]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_ilo(dom_lo[0]-1,j,k,1);
      const int bct = a_bc_ilo(dom_lo[0]-1,j,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        for (int n(0); n < nspecies_g; n++) {
          a_X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
          a_D_gk(i,j,k,n) = p_D_gk0[n];

          if (loc_advect_enthalpy) {
            a_cp_gk(i,j,k,n) = p_cp_gk0[n];
            a_h_gk(i,j,k,n) = p_cp_gk0[n] * p_bc_t_g[bcv];
          }
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

    ParallelFor(bx_yz_hi_3D, [a_bc_ihi,dom_hi,pinf,pout,minf,nspecies_g,a_X_gk,
        a_D_gk,a_cp_gk,a_h_gk,p_bc_X_gk,p_D_gk0,p_cp_gk0,p_bc_t_g,loc_advect_enthalpy]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_ihi(dom_hi[0]+1,j,k,1);
      const int bct = a_bc_ihi(dom_hi[0]+1,j,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        for (int n(0); n < nspecies_g; n++) {
          a_X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
          a_D_gk(i,j,k,n) = p_D_gk0[n];

          if (loc_advect_enthalpy) {
            a_cp_gk(i,j,k,n) = p_cp_gk0[n];
            a_h_gk(i,j,k,n) = p_cp_gk0[n] * p_bc_t_g[bcv];
          }
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

    ParallelFor(bx_xz_lo_3D, [a_bc_jlo,dom_lo,pinf,pout,minf,nspecies_g,a_X_gk,
        a_D_gk,a_cp_gk,a_h_gk,p_bc_X_gk,p_D_gk0,p_cp_gk0,p_bc_t_g,loc_advect_enthalpy]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_jlo(i,dom_lo[1]-1,k,1);
      const int bct = a_bc_jlo(i,dom_lo[1]-1,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        for (int n(0); n < nspecies_g; n++) {
          a_X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
          a_D_gk(i,j,k,n) = p_D_gk0[n];

          if (loc_advect_enthalpy) {
            a_cp_gk(i,j,k,n) = p_cp_gk0[n];
            a_h_gk(i,j,k,n) = p_cp_gk0[n] * p_bc_t_g[bcv];
          }
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

    ParallelFor(bx_xz_hi_3D, [a_bc_jhi,dom_hi,pinf,pout,minf,nspecies_g,a_X_gk,
        a_D_gk,a_cp_gk,a_h_gk,p_bc_X_gk,p_D_gk0,p_cp_gk0,p_bc_t_g,loc_advect_enthalpy]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_jhi(i,dom_hi[1]+1,k,1);
      const int bct = a_bc_jhi(i,dom_hi[1]+1,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        for (int n(0); n < nspecies_g; n++) {
          a_X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
          a_D_gk(i,j,k,n) = p_D_gk0[n];

          if (loc_advect_enthalpy) {
            a_cp_gk(i,j,k,n) = p_cp_gk0[n];
            a_h_gk(i,j,k,n) = p_cp_gk0[n] * p_bc_t_g[bcv];
          }
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

    ParallelFor(bx_xy_lo_3D, [a_bc_klo,dom_lo,pinf,pout,minf,nspecies_g,a_X_gk,
        a_D_gk,a_cp_gk,a_h_gk,p_bc_X_gk,p_D_gk0,p_cp_gk0,p_bc_t_g,loc_advect_enthalpy]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_klo(i,j,dom_lo[2]-1,1);
      const int bct = a_bc_klo(i,j,dom_lo[2]-1,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        for (int n(0); n < nspecies_g; n++) {
          a_X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
          a_D_gk(i,j,k,n) = p_D_gk0[n];

          if (loc_advect_enthalpy) {
            a_cp_gk(i,j,k,n) = p_cp_gk0[n];
            a_h_gk(i,j,k,n) = p_cp_gk0[n] * p_bc_t_g[bcv];
          }
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

    ParallelFor(bx_xy_hi_3D, [a_bc_khi,dom_hi,pinf,pout,minf,nspecies_g,a_X_gk,
        a_D_gk,a_cp_gk,a_h_gk,p_bc_X_gk,p_D_gk0,p_cp_gk0,p_bc_t_g,loc_advect_enthalpy]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_khi(i,j,dom_hi[2]+1,1);
      const int bct = a_bc_khi(i,j,dom_hi[2]+1,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
      {
        for (int n(0); n < nspecies_g; n++) {
          a_X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
          a_D_gk(i,j,k,n) = p_D_gk0[n];

          if (loc_advect_enthalpy) {
            a_cp_gk(i,j,k,n) = p_cp_gk0[n];
            a_h_gk(i,j,k,n) = p_cp_gk0[n] * p_bc_t_g[bcv];
          }
        }
      }
    });
  }

  Gpu::synchronize();

}
