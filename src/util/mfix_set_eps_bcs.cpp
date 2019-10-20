#include <mfix.H>
#include <bc_mod_F.H>

//
//  These subroutines set the BCs for the epsocity components only.
//

void
mfix::mfix_set_eps_bcs (const amrex::Vector< std::unique_ptr<MultiFab> > & eps_in) const
{
  BL_PROFILE("mfix::mfix_set_eps_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     // Set all values outside the domain to covered_val just to avoid use of undefined
     eps_in[lev]->setDomainBndry(covered_val,geom[lev]);

     eps_in[lev] -> FillBoundary (geom[lev].periodicity());
     Box domain(geom[lev].Domain());

     const int extrap_dir_bcs(1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*eps_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       set_eps_bcs(lev, (*eps_in[lev])[mfi], domain, &extrap_dir_bcs);

     EB_set_covered(*eps_in[lev], 0, eps_in[lev]->nComp(), eps_in[lev]->nGrow(), covered_val);

     // Do this after as well as before to pick up terms that got updated in the call above
     eps_in[lev] -> FillBoundary (geom[lev].periodicity());
  }
}

void
mfix::set_eps_bcs(const int lev,
                  FArrayBox& eps_fab,
                  const Box& domain,
                  const int* extrap_dir_bcs) const
{
  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<Real> const& eps = eps_fab.array();

  IntVect eps_lo(eps_fab.loVect());
  IntVect eps_hi(eps_fab.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  const int nlft = std::max(0, dom_lo[0]-eps_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-eps_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-eps_lo[2]);

  const int nrgt = std::max(0, eps_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, eps_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, eps_hi[2]-dom_hi[2]);

  // Create InVects for following 2D Boxes
  IntVect bx_yz_lo_lo_2D(eps_lo), bx_yz_lo_hi_2D(eps_hi);
  IntVect bx_yz_hi_lo_2D(eps_lo), bx_yz_hi_hi_2D(eps_hi);
  IntVect bx_xz_lo_lo_2D(eps_lo), bx_xz_lo_hi_2D(eps_hi);
  IntVect bx_xz_hi_lo_2D(eps_lo), bx_xz_hi_hi_2D(eps_hi);
  IntVect bx_xy_lo_lo_2D(eps_lo), bx_xy_lo_hi_2D(eps_hi);
  IntVect bx_xy_hi_lo_2D(eps_lo), bx_xy_hi_hi_2D(eps_hi);

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
  IntVect bx_yz_lo_hi_3D(eps_hi), bx_xz_lo_hi_3D(eps_hi), bx_xy_lo_hi_3D(eps_hi);
  IntVect bx_yz_hi_lo_3D(eps_lo), bx_xz_hi_lo_3D(eps_lo), bx_xy_hi_lo_3D(eps_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(eps_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, eps_hi);

  const Box bx_xz_lo_3D(eps_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, eps_hi);

  const Box bx_xy_lo_3D(eps_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, eps_hi);


  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  const amrex::Real* p_bc_u_g = m_bc_u_g.data();
  const amrex::Real* p_bc_v_g = m_bc_v_g.data();
  const amrex::Real* p_bc_w_g = m_bc_w_g.data();

  if (nlft > 0)
  {
    AMREX_FOR_3D(bx_yz_lo_3D, i, j, k,
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
        eps(i,j,k) = eps(dom_lo[0],j,k);
    });

    if(*extrap_dir_bcs > 0)
    {
      AMREX_FOR_3D(bx_yz_lo_2D, i, j, k,
      {
        const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

        if(bct == minf)
          eps(i,j,k) = 2*eps(i,j,k) - eps(i+1,j,k);
      });
    }
  }

  if (nrgt > 0)
  {
    AMREX_FOR_3D(bx_yz_hi_3D, i, j, k,
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
        eps(i,j,k) = eps(dom_hi[0],j,k);
    });

    if(*extrap_dir_bcs > 0)
    {
      AMREX_FOR_3D(bx_yz_hi_2D, i, j, k,
      {
        const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

        if(bct == minf)
          eps(i,j,k) = 2*eps(i,j,k) - eps(i-1,j,k);
      });
    }
  }

  if (nbot > 0)
  {
    AMREX_FOR_3D(bx_xz_lo_3D, i, j, k,
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
        eps(i,j,k) = eps(i,dom_lo[1],k);
    });

    if(*extrap_dir_bcs > 0)
    {

      AMREX_FOR_3D(bx_xz_lo_2D, i, j, k,
      {
        const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

        if(bct == minf)
          eps(i,j,k) = 2*eps(i,j,k) - eps(i,j+1,k);
      });
    }
  }

  if (ntop > 0)
  {
    AMREX_FOR_3D(bx_xz_hi_3D, i, j, k,
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
        eps(i,j,k) = eps(i,dom_hi[1],k);
    });

    if(*extrap_dir_bcs > 0)
    {

      AMREX_FOR_3D(bx_xz_hi_2D, i, j, k,
      {
        const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

        if(bct == minf)
          eps(i,j,k) = 2*eps(i,j,k) - eps(i,j-1,k);
      });
    }
  }

  if (ndwn > 0)
  {
    AMREX_FOR_3D(bx_xy_lo_3D, i, j, k,
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
        eps(i,j,k) = eps(i,j,dom_lo[2]);
    });

    if(*extrap_dir_bcs > 0)
    {

      AMREX_FOR_3D(bx_xy_lo_2D, i, j, k,
      {
        const int bct = bct_klo(i,j,dom_lo[2]-1,0);

        if(bct == minf)
          eps(i,j,k) = 2*eps(i,j,k) - eps(i,j,k+1);
      });
    }
  }

  if (nup > 0)
  {
    AMREX_FOR_3D(bx_xy_hi_3D, i, j, k,
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == pinf) or (bct == pout) or (bct == minf))
        eps(i,j,k) = eps(i,j,dom_hi[2]);
    });

    if(*extrap_dir_bcs > 0)
    {
      AMREX_FOR_3D(bx_xy_hi_2D, i, j, k,
      {
        const int bct = bct_khi(i,j,dom_hi[2]+1,0);

        if(bct == minf)
          eps(i,j,k) = 2*eps(i,j,k) - eps(i,j,k-1);
      });
    }
  }

  Gpu::synchronize();
}
