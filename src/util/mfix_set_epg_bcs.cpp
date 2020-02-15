#include <mfix.H>
#include <MFIX_LevelData.H>

//
//  These subroutines set the BCs for ep_g only
//

using namespace amrex;
using namespace std;

void
mfix::mfix_set_epg_bcs (const Vector< MultiFab* >& ep_g) const
{
  BL_PROFILE("mfix::mfix_set_epg_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
    // Set all values outside the domain to covered_val just to avoid use of
    // undefined
    ep_g[lev]->setDomainBndry(covered_val,geom[lev]);

    ep_g[lev]->FillBoundary(geom[lev].periodicity());
    Box domain(geom[lev].Domain());

    const int extrap_dir_bcs(1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ep_g[lev], false); mfi.isValid(); ++mfi)
      set_epg_bcs(lev, (*ep_g[lev])[mfi], domain, &extrap_dir_bcs);

    EB_set_covered(*(ep_g[lev]), 0, ep_g[lev]->nComp(),
      ep_g[lev]->nGrow(), covered_val);

    // Do this after as well as before to pick up terms that got updated in the
    // call above
    ep_g[lev]->FillBoundary(geom[lev].periodicity());
  }
}

void
mfix::set_epg_bcs (const int lev,
                   FArrayBox& epg_fab,
                   const Box& domain,
                   const int* extrap_dir_bcs) const
{
  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<Real> const& epg = epg_fab.array();

  IntVect epg_lo(epg_fab.loVect());
  IntVect epg_hi(epg_fab.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  const int nlft = std::max(0, dom_lo[0]-epg_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-epg_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-epg_lo[2]);

  const int nrgt = std::max(0, epg_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, epg_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, epg_hi[2]-dom_hi[2]);

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

  // Create 2D boxes for CUDA loops
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

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(epg_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, epg_hi);

  const Box bx_xz_lo_3D(epg_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, epg_hi);

  const Box bx_xy_lo_3D(epg_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, epg_hi);

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  amrex::Real* p_bc_ep_g = m_bc_ep_g.data();

  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D, 
      [bct_ilo,dom_lo,minf,pinf,pout,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);

      if((bct == pinf) or (bct == pout))
        epg(i,j,k) = epg(dom_lo[0],j,k);
      else if (bct == minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*extrap_dir_bcs > 0)
    {
      amrex::ParallelFor(bx_yz_lo_2D, 
        [bct_ilo,dom_lo,minf,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

        if(bct == minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i+1,j,k);
      });
    }
  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D, 
      [bct_ihi,dom_hi,minf,pinf,pout,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);

      if((bct == pinf) or (bct == pout))
        epg(i,j,k) = epg(dom_hi[0],j,k);
      else if (bct == minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*extrap_dir_bcs > 0)
    {
      amrex::ParallelFor(bx_yz_hi_2D, 
        [bct_ihi,dom_hi,minf,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

        if(bct == minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i-1,j,k);
      });
    }
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D, 
      [bct_jlo,dom_lo,minf,pinf,pout,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);

      if((bct == pinf) or (bct == pout))
        epg(i,j,k) = epg(i,dom_lo[1],k);
      else if (bct == minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*extrap_dir_bcs > 0)
    {

      amrex::ParallelFor(bx_xz_lo_2D, 
        [bct_jlo,dom_lo,minf,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

        if(bct == minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i,j+1,k);
      });
    }
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D, 
      [bct_jhi,dom_hi,minf,pinf,pout,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);

      if((bct == pinf) or (bct == pout))
        epg(i,j,k) = epg(i,dom_hi[1],k);
      else if (bct == minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*extrap_dir_bcs > 0)
    {

      amrex::ParallelFor(bx_xz_hi_2D, 
        [bct_jhi,dom_hi,minf,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

        if(bct == minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i,j-1,k);
      });
    }
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D, 
      [bct_klo,dom_lo,minf,pinf,pout,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);

      if((bct == pinf) or (bct == pout))
        epg(i,j,k) = epg(i,j,dom_lo[2]);
      else if (bct == minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*extrap_dir_bcs > 0)
    {

      amrex::ParallelFor(bx_xy_lo_2D, 
        [bct_klo,dom_lo,minf,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_klo(i,j,dom_lo[2]-1,0);

        if(bct == minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i,j,k+1);
      });
    }
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D, 
      [bct_khi,dom_hi,minf,pinf,pout,p_bc_ep_g,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);

      if((bct == pinf) or (bct == pout))
        epg(i,j,k) = epg(i,j,dom_hi[2]);
      else if (bct == minf)
        epg(i,j,k) = p_bc_ep_g[bcv];
    });

    if(*extrap_dir_bcs > 0)
    {
      amrex::ParallelFor(bx_xy_hi_2D, 
        [bct_khi,dom_hi,minf,epg] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_khi(i,j,dom_hi[2]+1,0);

        if(bct == minf)
          epg(i,j,k) = 2*epg(i,j,k) - epg(i,j,k-1);
      });
    }
  }
}
