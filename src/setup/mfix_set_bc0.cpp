#include <mfix.H>
#include <eos_mod.H>
#include <param_mod_F.H>

#include <MFIX_FLUID_Parms.H>

void 
mfix::set_bc0 (const Box& sbx,
               MFIter* mfi,
               const int lev,
               const Box& domain)
{
  const Real ro_g0  = FLUID::ro_g0;
  const Real mu_g0  = FLUID::mu_g0;
  const Real trac_0 = FLUID::trac_0;

  Array4<Real> const& a_ep_g = m_leveldata[lev]->ep_g->array(*mfi);
  Array4<Real> const& a_ro_g = ro_g[lev]->array(*mfi);
  Array4<Real> const& a_trac = trac[lev]->array(*mfi);
  Array4<Real> const& a_mu_g = mu_g[lev]->array(*mfi);

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

  const int nlft = std::max(0,dom_lo[0]-sbx_lo[0]);
  const int nbot = std::max(0,dom_lo[1]-sbx_lo[1]);
  const int ndwn = std::max(0,dom_lo[2]-sbx_lo[2]);

  const int nrgt = std::max(0,sbx_hi[0]-dom_hi[0]);
  const int ntop = std::max(0,sbx_hi[1]-dom_hi[1]);
  const int nup  = std::max(0,sbx_hi[2]-dom_hi[2]);

  const Real undefined = get_undefined();

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  amrex::Real* p_bc_ep_g = m_bc_ep_g.data();
  amrex::Real* p_bc_t_g = m_bc_t_g.data();

  if (nlft > 0)
  {
    IntVect bx_yz_lo_hi_3D(sbx_hi);
  
    // Fix lo and hi limits
    bx_yz_lo_hi_3D[0] = dom_lo[0]-1;

    const Box bx_yz_lo_3D(sbx_lo, bx_yz_lo_hi_3D);

    amrex::ParallelFor(bx_yz_lo_3D,
        [a_bc_ilo,dom_lo,pinf,pout,minf,ro_g0,trac_0,mu_g0,undefined,
         p_bc_t_g,p_bc_ep_g,a_ep_g,a_ro_g,a_trac,a_mu_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = a_bc_ilo(dom_lo[0]-1,j,k,1);
          const int bct = a_bc_ilo(dom_lo[0]-1,j,k,0);

          if((bct == pinf) or (bct == pout) or (bct == minf))
          {
            Real bc_ro_g(ro_g0);
            Real bc_trac(trac_0);
            Real bc_mu_g(0);

            if (is_equal(mu_g0, undefined))
              bc_mu_g = sutherland(p_bc_t_g[bcv]);
            else
              bc_mu_g = mu_g0;

            a_ep_g(i,j,k) = p_bc_ep_g[bcv];
            a_ro_g(i,j,k) = bc_ro_g;
            a_trac(i,j,k) = bc_trac;
            a_mu_g(i,j,k) = bc_mu_g;
          }
        });
  }
  
  if (nrgt > 0)
  {
    IntVect bx_yz_hi_lo_3D(sbx_lo);

    // Fix lo and hi limits
    bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

    const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, sbx_hi);
    
    amrex::ParallelFor(bx_yz_hi_3D,
        [a_bc_ihi,dom_hi,pinf,pout,minf,ro_g0,trac_0,mu_g0,undefined,
         p_bc_t_g,p_bc_ep_g,a_ep_g,a_ro_g,a_trac,a_mu_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = a_bc_ihi(dom_hi[0]+1,j,k,1);
          const int bct = a_bc_ihi(dom_hi[0]+1,j,k,0);

          if((bct == pinf) or (bct == pout) or (bct == minf))
          {
            Real bc_ro_g(ro_g0);
            Real bc_trac(trac_0);
            Real bc_mu_g(0);

            if (is_equal(mu_g0, undefined))
              bc_mu_g = sutherland(p_bc_t_g[bcv]);
            else
              bc_mu_g = mu_g0;

            a_ep_g(i,j,k) = p_bc_ep_g[bcv];
            a_ro_g(i,j,k) = bc_ro_g;
            a_trac(i,j,k) = bc_trac;
              a_mu_g(i,j,k) = bc_mu_g;
          }
        });
  }

  if (nbot > 0)
  {
    IntVect bx_xz_lo_hi_3D(sbx_hi);
  
    // Fix lo and hi limits
    bx_xz_lo_hi_3D[1] = dom_lo[1]-1;

    const Box bx_xz_lo_3D(sbx_lo, bx_xz_lo_hi_3D);

    amrex::ParallelFor(bx_xz_lo_3D,
        [a_bc_jlo,dom_lo,pinf,pout,minf,ro_g0,trac_0,mu_g0,undefined,
         p_bc_t_g,p_bc_ep_g,a_ep_g,a_ro_g,a_trac,a_mu_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = a_bc_jlo(i,dom_lo[1]-1,k,1);
          const int bct = a_bc_jlo(i,dom_lo[1]-1,k,0);

          if((bct == pinf) or (bct == pout) or (bct == minf))
          {
            Real bc_ro_g(ro_g0);
            Real bc_trac(trac_0);
            Real bc_mu_g(0);

            if (is_equal(mu_g0, undefined))
               bc_mu_g = sutherland(p_bc_t_g[bcv]);
            else
               bc_mu_g = mu_g0;
      
            a_ep_g(i,j,k) = p_bc_ep_g[bcv];
            a_ro_g(i,j,k) = bc_ro_g;
            a_trac(i,j,k) = bc_trac;
              a_mu_g(i,j,k) = bc_mu_g;
          }
        });
  }

  if (ntop > 0)
  {
    IntVect bx_xz_hi_lo_3D(sbx_lo);
  
    // Fix lo and hi limits
    bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

    const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, sbx_hi);

    amrex::ParallelFor(bx_xz_hi_3D,
        [a_bc_jhi,dom_hi,pinf,pout,minf,ro_g0,trac_0,mu_g0,undefined,
         p_bc_t_g,p_bc_ep_g,a_ep_g,a_ro_g,a_trac,a_mu_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = a_bc_jhi(i,dom_hi[1]+1,k,1);
          const int bct = a_bc_jhi(i,dom_hi[1]+1,k,0);

          if((bct == pinf) or (bct == pout) or (bct == minf))
          {
            Real bc_ro_g(ro_g0);
            Real bc_trac(trac_0);
            Real bc_mu_g(0);

            if (is_equal(mu_g0, undefined))
               bc_mu_g = sutherland(p_bc_t_g[bcv]);
            else
               bc_mu_g = mu_g0;

             a_ep_g(i,j,k) = p_bc_ep_g[bcv];
            a_ro_g(i,j,k) = bc_ro_g;
            a_trac(i,j,k) = bc_trac;
              a_mu_g(i,j,k) = bc_mu_g;
          }
        });
  }

  if (ndwn > 0)
  {
    IntVect bx_xy_lo_hi_3D(sbx_hi);
  
    // Fix lo and hi limits
    bx_xy_lo_hi_3D[2] = dom_lo[2]-1;

    const Box bx_xy_lo_3D(sbx_lo, bx_xy_lo_hi_3D);
    
    amrex::ParallelFor(bx_xy_lo_3D,
        [a_bc_klo,dom_lo,pinf,pout,minf,ro_g0,trac_0,mu_g0,undefined,
         p_bc_t_g,p_bc_ep_g,a_ep_g,a_ro_g,a_trac,a_mu_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = a_bc_klo(i,j,dom_lo[2]-1,1);
          const int bct = a_bc_klo(i,j,dom_lo[2]-1,0);

          if((bct == pinf) or (bct == pout) or (bct == minf))
          {
            Real bc_ro_g(ro_g0);
            Real bc_trac(trac_0);
            Real bc_mu_g(0);

            if (is_equal(mu_g0, undefined))
               bc_mu_g = sutherland(p_bc_t_g[bcv]);
            else
               bc_mu_g = mu_g0;

            a_ep_g(i,j,k) = p_bc_ep_g[bcv];
            a_ro_g(i,j,k) = bc_ro_g;
            a_trac(i,j,k) = bc_trac;
              a_mu_g(i,j,k) = bc_mu_g;
          }
        });
  }

  if (nup > 0)
  {
    IntVect bx_xy_hi_lo_3D(sbx_lo);

    // Fix lo and hi limits
    bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

    const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, sbx_hi);
    
    amrex::ParallelFor(bx_xy_hi_3D,
        [a_bc_khi,dom_hi,pinf,pout,minf,ro_g0,trac_0,mu_g0,undefined,
         p_bc_t_g,p_bc_ep_g,a_ep_g,a_ro_g,a_trac,a_mu_g]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = a_bc_khi(i,j,dom_hi[2]+1,1);
          const int bct = a_bc_khi(i,j,dom_hi[2]+1,0);

          if((bct == pinf) or (bct == pout) or (bct == minf))
          {
            Real bc_ro_g(ro_g0);
            Real bc_trac(trac_0);
            Real bc_mu_g(0);

            if (is_equal(mu_g0, undefined))
               bc_mu_g = sutherland(p_bc_t_g[bcv]);
            else
               bc_mu_g = mu_g0;

            a_ep_g(i,j,k) = p_bc_ep_g[bcv];
            a_ro_g(i,j,k) = bc_ro_g;
            a_trac(i,j,k) = bc_trac;
              a_mu_g(i,j,k) = bc_mu_g;
          }
        });
  }
}
