#include <mfix.H>
#include <bc_mod_F.H>
#include <eos_mod.hpp>
#include <fld_constants_mod_F.H>
#include <param_mod_F.H>

using namespace amrex;

//
// Set the BCs for all the variables EXCEPT pressure or velocity.
//
void
mfix::mfix_set_scalar_bcs (Real time,
                           Vector< std::unique_ptr<MultiFab> > & ro_g_in,
                           Vector< std::unique_ptr<MultiFab> > & trac_in,
                           Vector< std::unique_ptr<MultiFab> > & ep_g_in,
                           Vector< std::unique_ptr<MultiFab> > & mu_g_in)
{
  BL_PROFILE("mfix::mfix_set_scalar_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*ep_g[lev], true); mfi.isValid(); ++mfi)
     {
        set_scalar_bcs(time, lev, (*ro_g_in[lev])[mfi], 0, domain);
        set_scalar_bcs(time, lev, (*trac_in[lev])[mfi], 1, domain);
        set_scalar_bcs(time, lev, (*ep_g_in[lev])[mfi], 2, domain);
        set_scalar_bcs(time, lev, (*mu_g_in[lev])[mfi], 3, domain);
     }

     ro_g_in[lev] -> FillBoundary (geom[lev].periodicity());
     ep_g_in[lev] -> FillBoundary (geom[lev].periodicity());
     trac_in[lev] -> FillBoundary (geom[lev].periodicity());
     mu_g_in[lev] -> FillBoundary (geom[lev].periodicity());

     EB_set_covered(*ro_g_in[lev], 0, ro_g_in[lev]->nComp(), ro_g_in[lev]->nGrow(), covered_val);
     EB_set_covered(*trac_in[lev], 0, trac_in[lev]->nComp(), trac_in[lev]->nGrow(), covered_val);
     EB_set_covered(*ep_g_in[lev], 0, ep_g_in[lev]->nComp(), ep_g_in[lev]->nGrow(), covered_val);
     EB_set_covered(*mu_g_in[lev], 0, mu_g_in[lev]->nComp(), mu_g_in[lev]->nGrow(), covered_val);
  }
}

void 
mfix::set_scalar_bcs(Real time,
                     const int lev,
                     FArrayBox& scal_fab,
                     int comp,
                     const Box& domain)

{
  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  Array4<Real> const& scal_arr = scal_fab.array();

  Real bc0;

  if (comp == 0) {
     bc0 = get_ro_g0();
  } else if (comp == 1) {
     bc0 = get_trac0();
  } else if (comp == 2) {
    // We don't have a default value for ep_g
  } else if (comp == 3) {
     bc0 = get_mu_g0();
  }

  IntVect scal_lo(scal_fab.loVect());
  IntVect scal_hi(scal_fab.hiVect());

  const int nlft = std::max(0, dom_lo[0]-scal_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-scal_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-scal_lo[2]);

  const int nrgt = std::max(0, scal_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, scal_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, scal_hi[2]-dom_hi[2]);

  // Create InVects for following 2D Boxes
  IntVect bx_yz_lo_lo_2D(scal_lo), bx_yz_lo_hi_2D(scal_hi);
  IntVect bx_yz_hi_lo_2D(scal_lo), bx_yz_hi_hi_2D(scal_hi);
  IntVect bx_xz_lo_lo_2D(scal_lo), bx_xz_lo_hi_2D(scal_hi);
  IntVect bx_xz_hi_lo_2D(scal_lo), bx_xz_hi_hi_2D(scal_hi);
  IntVect bx_xy_lo_lo_2D(scal_lo), bx_xy_lo_hi_2D(scal_hi);
  IntVect bx_xy_hi_lo_2D(scal_lo), bx_xy_hi_hi_2D(scal_hi);

  // Fix lo and hi limits
  bx_yz_lo_lo_2D[0] = dom_lo[0]-1;
  bx_yz_lo_hi_2D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_2D[0] = dom_hi[0]+1;
  bx_yz_hi_hi_2D[0] = dom_hi[0]+1;

  bx_xz_lo_lo_2D[1] = dom_lo[1]-1;
  bx_xz_lo_hi_2D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_2D[1] = dom_hi[1]+1;
  bx_xz_hi_hi_2D[1] = dom_hi[1]+1;

  bx_xy_lo_lo_2D[2] = dom_lo[2]-1;
  bx_xy_lo_hi_2D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_2D[2] = dom_hi[2]+1;
  bx_xy_hi_hi_2D[2] = dom_hi[2]+1;

  // Create 2D boxes for CUDA loops
  const Box bx_yz_lo_2D(bx_yz_lo_lo_2D, bx_yz_lo_hi_2D);
  const Box bx_yz_hi_2D(bx_yz_hi_lo_2D, bx_yz_hi_hi_2D);

  const Box bx_xz_lo_2D(bx_xz_lo_lo_2D, bx_xz_lo_hi_2D);
  const Box bx_xz_hi_2D(bx_xz_hi_lo_2D, bx_xz_hi_hi_2D);

  const Box bx_xy_lo_2D(bx_xy_lo_lo_2D, bx_xy_lo_hi_2D);
  const Box bx_xy_hi_2D(bx_xy_hi_lo_2D, bx_xy_hi_hi_2D);

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

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(scal_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, scal_hi);

  const Box bx_xz_lo_3D(scal_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, scal_hi);

  const Box bx_xy_lo_3D(scal_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, scal_hi);

  const Real undefined = get_undefined();

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  amrex::Real* p_bc_ep_g = m_bc_ep_g.data();
  amrex::Real* p_bc_t_g  = m_bc_t_g.data();

  if (nlft > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_lo_3D, i, j, k,
    {
      Real bc_scal(bc0);
      int my_comp(comp);

      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if ((bct == pinf) or (bct == pout))
      {
        scal_arr(i,j,k) = scal_arr(dom_lo[0],j,k);
      }
      else if (bct == minf)
      {
        if (my_comp == 3 && is_equal(bc0, undefined))
          bc_scal = sutherland(p_bc_t_g[bcv]);
        else
          bc_scal = bc0;

        if (my_comp != 2)
           scal_arr(i,j,k) = bc_scal;
      }
    });

    if (comp == 2)
    {
       AMREX_HOST_DEVICE_FOR_3D(bx_yz_lo_2D, i, j, k,
       {
         const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
         const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
   
         if(bct == minf)
           scal_arr(i,j,k) = 2*p_bc_ep_g[bcv] - scal_arr(i+1,j,k);
       });
     }
  }

  if (nrgt > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_hi_3D, i, j, k,
    {
      Real bc_scal(bc0);
      int my_comp(comp);

      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == pinf) or (bct == pout))
      {
        scal_arr(i,j,k) = scal_arr(dom_hi[0],j,k);
      }
      else if(bct == minf)
      {
        if (my_comp == 3 && is_equal(bc0, undefined))
          bc_scal = sutherland(p_bc_t_g[bcv]);
        else
          bc_scal = bc0;

        if (my_comp != 2)
           scal_arr(i,j,k) = bc_scal;
      }
    });

    if (comp == 2)
    {
       AMREX_HOST_DEVICE_FOR_3D(bx_yz_hi_2D, i, j, k,
       {
         const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
         const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
   
         if(bct == minf)
           scal_arr(i,j,k) = 2*p_bc_ep_g[bcv] - scal_arr(i-1,j,k);
       });
     }
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (nbot > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xz_lo_3D, i, j, k,
    {
      Real bc_scal(bc0);
      int my_comp(comp);

      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == pinf) or (bct == pout))
      {
        scal_arr(i,j,k) = scal_arr(i,dom_lo[1],k);
      }
      else if(bct == minf)
      {
        if (my_comp == 3 && is_equal(bc0, undefined))
          bc_scal = sutherland(p_bc_t_g[bcv]);
        else
          bc_scal = bc0;

        if (my_comp != 2)
           scal_arr(i,j,k) = bc_scal;
      }
    });

    if (comp == 2)
    {
       AMREX_HOST_DEVICE_FOR_3D(bx_xz_lo_2D, i, j, k,
       {
         const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
         const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
   
         if(bct == minf)
           scal_arr(i,j,k) = 2*p_bc_ep_g[bcv] - scal_arr(i,j+1,k);
       });
     }
  }

  if (ntop > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xz_hi_3D, i, j, k,
    {
      Real bc_scal(bc0);
      int my_comp(comp);

      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == pinf) or (bct == pout))
      {
        scal_arr(i,j,k) = scal_arr(i,dom_hi[1],k);
      }
      else if(bct == minf)
      {
        if (my_comp == 3 && is_equal(bc0, undefined))
          bc_scal = sutherland(p_bc_t_g[bcv]);
        else
          bc_scal = bc0;

        if (my_comp != 2)
           scal_arr(i,j,k) = bc_scal;
      }
    });

    if (comp == 2)
    {
       AMREX_HOST_DEVICE_FOR_3D(bx_xz_hi_2D, i, j, k,
       {
         const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
         const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

         if(bct == minf)
           scal_arr(i,j,k) = 2*p_bc_ep_g[bcv] - scal_arr(i,j-1,k);
       });
     }
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (ndwn > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xy_lo_3D, i, j, k,
    {
      Real bc_scal(bc0);
      int my_comp(comp);

      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == pinf) or (bct == pout))
      {
        scal_arr(i,j,k) = scal_arr(i,j,dom_lo[2]);
      }
      else if(bct == minf)
      {
        if (my_comp == 3 && is_equal(bc0, undefined))
          bc_scal = sutherland(p_bc_t_g[bcv]);
        else
          bc_scal = bc0;

        if (my_comp != 2)
           scal_arr(i,j,k) = bc_scal;
      }
    });

    if (comp == 2)
    {
       AMREX_HOST_DEVICE_FOR_3D(bx_xy_lo_2D, i, j, k,
       {
         const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
         const int bct = bct_klo(i,j,dom_lo[2]-1,0);
   
         if(bct == minf)
           scal_arr(i,j,k) = 2*p_bc_ep_g[bcv] - scal_arr(i,j,k+1);
       });
     }
  }

  if (nup > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xy_hi_3D, i, j, k,
    {
      Real bc_scal(bc0);
      int my_comp(comp);

      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if ((bct == pinf) or (bct == pout))
      {
        scal_arr(i,j,k) = scal_arr(i,j,dom_hi[2]);
      }
      else if (bct == minf)
      {
        if (my_comp == 3 && is_equal(bc0, undefined))
          bc_scal = sutherland(p_bc_t_g[bcv]);
        else
          bc_scal = bc0;

        if (my_comp != 2)
           scal_arr(i,j,k) = bc_scal;
      }
    });

    if (comp == 2)
    {
       AMREX_HOST_DEVICE_FOR_3D(bx_xy_hi_2D, i, j, k,
       {
         const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
         const int bct = bct_khi(i,j,dom_hi[2]+1,0);
   
         if(bct == minf)
           scal_arr(i,j,k) = 2*p_bc_ep_g[bcv] - scal_arr(i,j,k-1);
       });
    }
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif
}
