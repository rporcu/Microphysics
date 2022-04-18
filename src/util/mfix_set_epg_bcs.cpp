#include <mfix.H>

//
//  These subroutines set the BCs for ep_g only
//

using namespace amrex;

void
mfix::mfix_set_epg_bcs (const Vector< MultiFab* >& epg_in, const int dir_bc) const
{
  BL_PROFILE("mfix::mfix_set_epg_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     // Set all values outside the domain to covered_val just to avoid use of undefined
     epg_in[lev]->setDomainBndry(covered_val,geom[lev]);

     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*epg_in[lev], false); mfi.isValid(); ++mfi)
       set_epg_bcs(lev, (*epg_in[lev])[mfi], domain, &dir_bc);

     EB_set_covered(*epg_in[lev], 0, epg_in[lev]->nComp(), epg_in[lev]->nGrow(), covered_val);

     // Do this after as well as before to pick up terms that got updated in the call above
     epg_in[lev]->FillBoundary(geom[lev].periodicity());
  }
}

void
mfix::set_epg_bcs (const int lev,
                   FArrayBox& epg_fab,
                   const Box& domain,
                   const int* dir_bc) const
{
  BL_PROFILE("mfix::set_epg_bcs()");

  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<Real> const& epg = epg_fab.array();

  IntVect epg_lo(epg_fab.loVect());
  IntVect epg_hi(epg_fab.hiVect());

  Array4<const int> const& bct_ilo = bc_list.bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_list.bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_list.bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_list.bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_list.bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_list.bc_khi[lev]->array();

  const int nlft = amrex::max(0, dom_lo[0]-epg_lo[0]);
  const int nbot = amrex::max(0, dom_lo[1]-epg_lo[1]);
  const int ndwn = amrex::max(0, dom_lo[2]-epg_lo[2]);

  const int nrgt = amrex::max(0, epg_hi[0]-dom_hi[0]);
  const int ntop = amrex::max(0, epg_hi[1]-dom_hi[1]);
  const int nup  = amrex::max(0, epg_hi[2]-dom_hi[2]);

  // Create InVects for following 2D Boxes
  IntVect bx_yz_lo_lo_2D(epg_lo), bx_yz_lo_hi_2D(epg_hi);
  IntVect bx_yz_hi_lo_2D(epg_lo), bx_yz_hi_hi_2D(epg_hi);
  IntVect bx_xz_lo_lo_2D(epg_lo), bx_xz_lo_hi_2D(epg_hi);
  IntVect bx_xz_hi_lo_2D(epg_lo), bx_xz_hi_hi_2D(epg_hi);
  IntVect bx_xy_lo_lo_2D(epg_lo), bx_xy_lo_hi_2D(epg_hi);
  IntVect bx_xy_hi_lo_2D(epg_lo), bx_xy_hi_hi_2D(epg_hi);

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

  // Create 2D boxes for GPU loops
  const Box bx_yz_lo_2D(bx_yz_lo_lo_2D, bx_yz_lo_hi_2D);
  const Box bx_yz_hi_2D(bx_yz_hi_lo_2D, bx_yz_hi_hi_2D);

  const Box bx_xz_lo_2D(bx_xz_lo_lo_2D, bx_xz_lo_hi_2D);
  const Box bx_xz_hi_2D(bx_xz_hi_lo_2D, bx_xz_hi_hi_2D);

  const Box bx_xy_lo_2D(bx_xy_lo_lo_2D, bx_xy_lo_hi_2D);
  const Box bx_xy_hi_2D(bx_xy_hi_lo_2D, bx_xy_hi_hi_2D);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(epg_hi), bx_xz_lo_hi_3D(epg_hi), bx_xy_lo_hi_3D(epg_hi);
  IntVect bx_yz_hi_lo_3D(epg_lo), bx_xz_hi_lo_3D(epg_lo), bx_xy_hi_lo_3D(epg_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for GPU loops
  const Box bx_yz_lo_3D(epg_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, epg_hi);

  const Box bx_xz_lo_3D(epg_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, epg_hi);

  const Box bx_xy_lo_3D(epg_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, epg_hi);

  const Real* p_bc_ep_g = m_bc_ep_g.data();

  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D,
      [bct_ilo,dom_lo,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);

      if((bct == BCList::pinf) || (bct == BCList::pout))
        epg(i,j,k) = epg(dom_lo[0],j,k);
      else if (bct == BCList::minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*dir_bc == 1 )
    {
      amrex::ParallelFor(bx_yz_lo_2D,
        [bct_ilo,dom_lo,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

        if(bct == BCList::minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i+1,j,k);
      });
    }
    else if(*dir_bc == 2)
    {
      amrex::ParallelFor(bx_yz_lo_2D,
        [bct_ilo,dom_lo,epg]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

        if(bct == BCList::minf)
          epg(i,j,k) = epg(dom_lo[0],j,k);
      });
    }
  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D,
      [bct_ihi,dom_hi,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);

      if((bct == BCList::pinf) || (bct == BCList::pout))
        epg(i,j,k) = epg(dom_hi[0],j,k);
      else if (bct == BCList::minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*dir_bc == 1)
    {
      amrex::ParallelFor(bx_yz_hi_2D,
        [bct_ihi,dom_hi,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

        if(bct == BCList::minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i-1,j,k);
      });
    }
    else if(*dir_bc == 2)
    {
      amrex::ParallelFor(bx_yz_hi_2D,
        [bct_ihi,dom_hi,epg]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

        if(bct == BCList::minf)
          epg(i,j,k) = epg(dom_hi[0],j,k);
      });
    }
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D,
      [bct_jlo,dom_lo,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);

      if((bct == BCList::pinf) || (bct == BCList::pout))
        epg(i,j,k) = epg(i,dom_lo[1],k);
      else if (bct == BCList::minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*dir_bc == 1)
    {

      amrex::ParallelFor(bx_xz_lo_2D,
        [bct_jlo,dom_lo,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

        if(bct == BCList::minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i,j+1,k);
      });
    }
    else if (*dir_bc == 2)
    {
      amrex::ParallelFor(bx_xz_lo_2D,
        [bct_jlo,dom_lo,epg]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

        if(bct == BCList::minf)
          epg(i,j,k) = epg(i,dom_lo[1],k);
      });
    }
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D,
      [bct_jhi,dom_hi,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);

      if((bct == BCList::pinf) || (bct == BCList::pout))
        epg(i,j,k) = epg(i,dom_hi[1],k);
      else if (bct == BCList::minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*dir_bc == 1)
    {

      amrex::ParallelFor(bx_xz_hi_2D,
        [bct_jhi,dom_hi,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

        if(bct == BCList::minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i,j-1,k);
      });
    }
    else if (*dir_bc == 2)
    {
      amrex::ParallelFor(bx_xz_hi_2D,
        [bct_jhi,dom_hi,epg]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

        if(bct == BCList::minf)
          epg(i,j,k) = epg(i,dom_hi[1],k);
      });
    }
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D,
      [bct_klo,dom_lo,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);

      if((bct == BCList::pinf) || (bct == BCList::pout))
        epg(i,j,k) = epg(i,j,dom_lo[2]);
      else if (bct == BCList::minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*dir_bc == 1)
    {

      amrex::ParallelFor(bx_xy_lo_2D,
        [bct_klo,dom_lo,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_klo(i,j,dom_lo[2]-1,0);

        if(bct == BCList::minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i,j,k+1);
      });
    }
    else if (*dir_bc == 2)
    {
      amrex::ParallelFor(bx_xy_lo_2D,
        [bct_klo,dom_lo,epg]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_klo(i,j,dom_lo[2]-1,0);

        if(bct == BCList::minf)
          epg(i,j,k) = epg(i,j,dom_lo[2]);
      });
    }
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D,
      [bct_khi,dom_hi,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);

      if((bct == BCList::pinf) || (bct == BCList::pout))
        epg(i,j,k) = epg(i,j,dom_hi[2]);
      else if (bct == BCList::minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*dir_bc == 1)
    {
      amrex::ParallelFor(bx_xy_hi_2D,
        [bct_khi,dom_hi,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_khi(i,j,dom_hi[2]+1,0);

        if(bct == BCList::minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i,j,k-1);
      });
    }
    else if (*dir_bc == 2)
    {
      amrex::ParallelFor(bx_xy_hi_2D,
        [bct_khi,dom_hi,epg]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_khi(i,j,dom_hi[2]+1,0);

        if(bct == BCList::minf)
          epg(i,j,k) = epg(i,j,dom_hi[2]);
      });
    }

  }
}
