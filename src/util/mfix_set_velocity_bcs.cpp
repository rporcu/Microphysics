#include <mfix.H>
#include <bc_mod_F.H>

//
//  These subroutines set the BCs for the velocity components only.
//

void
mfix::mfix_set_velocity_bcs (Real time,
                             Vector< std::unique_ptr<MultiFab> > & vel_in,
                             int extrap_dir_bcs) const
{
  BL_PROFILE("mfix::mfix_set_velocity_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     // Set all values outside the domain to covered_val just to avoid use of undefined
     vel_in[lev]->setDomainBndry(covered_val,geom[lev]);

     vel_in[lev] -> FillBoundary (geom[lev].periodicity());
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        set_velocity_bcs(time, lev, (*vel_in[lev])[mfi], domain, &extrap_dir_bcs);

     EB_set_covered(*vel_in[lev], 0, vel_in[lev]->nComp(), vel_in[lev]->nGrow(), covered_val);

     // Do this after as well as before to pick up terms that got updated in the call above
     vel_in[lev] -> FillBoundary (geom[lev].periodicity());
  }
}

void
mfix::set_velocity_bcs(Real time,
                       const int lev,
                       FArrayBox& vel_fab,
                       const Box& domain,
                       const int* extrap_dir_bcs) const
{
  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<Real> const& vel = vel_fab.array();

  IntVect vel_lo(vel_fab.loVect());
  IntVect vel_hi(vel_fab.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  const int nlft = std::max(0, dom_lo[0]-vel_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-vel_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-vel_lo[2]);

  const int nrgt = std::max(0, vel_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, vel_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, vel_hi[2]-dom_hi[2]);

  // Create InVects for following 2D Boxes
  IntVect bx_yz_lo_lo_2D(vel_lo), bx_yz_lo_hi_2D(vel_hi);
  IntVect bx_yz_hi_lo_2D(vel_lo), bx_yz_hi_hi_2D(vel_hi);
  IntVect bx_xz_lo_lo_2D(vel_lo), bx_xz_lo_hi_2D(vel_hi);
  IntVect bx_xz_hi_lo_2D(vel_lo), bx_xz_hi_hi_2D(vel_hi);
  IntVect bx_xy_lo_lo_2D(vel_lo), bx_xy_lo_hi_2D(vel_hi);
  IntVect bx_xy_hi_lo_2D(vel_lo), bx_xy_hi_hi_2D(vel_hi);

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
  IntVect bx_yz_lo_hi_3D(vel_hi), bx_xz_lo_hi_3D(vel_hi), bx_xy_lo_hi_3D(vel_hi);
  IntVect bx_yz_hi_lo_3D(vel_lo), bx_xz_hi_lo_3D(vel_lo), bx_xy_hi_lo_3D(vel_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(vel_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, vel_hi);

  const Box bx_xz_lo_3D(vel_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, vel_hi);

  const Box bx_xy_lo_3D(vel_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, vel_hi);

  mfix_usr1_cpp(time);

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  const amrex::Real* p_bc_u_g = m_bc_u_g.data();
  const amrex::Real* p_bc_v_g = m_bc_v_g.data();
  const amrex::Real* p_bc_w_g = m_bc_w_g.data();

  if (nlft > 0)
  {
    AMREX_FOR_4D(bx_yz_lo_3D, 3, i, j, k, n,
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == pinf) or (bct == pout))
        vel(i,j,k,n) = vel(dom_lo[0],j,k,n);
      else if(bct == minf)
      {
        if(n == 0)
          vel(i,j,k,n) = p_bc_u_g[bcv];
        else
          vel(i,j,k,n) = 0;
      }
    });

    if(*extrap_dir_bcs > 0)
    {
      AMREX_FOR_4D(bx_yz_lo_2D, 3, i, j, k, n,
      {
        const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

        if(bct == minf)
          vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i+1,j,k,n);
      });
    }
  }

  if (nrgt > 0)
  {
    AMREX_FOR_4D(bx_yz_hi_3D, 3, i, j, k, n,
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == pinf) or (bct == pout))
        vel(i,j,k,n) = vel(dom_hi[0],j,k,n);
      else if(bct == minf)
      {
        if(n == 0)
          vel(i,j,k,n) = p_bc_u_g[bcv];
        else
          vel(i,j,k,n) = 0;
      }
    });

    if(*extrap_dir_bcs > 0)
    {
      AMREX_FOR_4D(bx_yz_hi_2D, 3, i, j, k, n,
      {
        const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

        if(bct == minf)
          vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i-1,j,k,n);
      });
    }
  }

  if (nbot > 0)
  {
    AMREX_FOR_4D(bx_xz_lo_3D, 3, i, j, k, n,
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == pinf) or (bct == pout))
        vel(i,j,k,n) = vel(i,dom_lo[1],k,n);
      else if(bct == minf)
      {
        if(n == 1)
          vel(i,j,k,n) = p_bc_v_g[bcv];
        else
          vel(i,j,k,n) = 0;
      }
    });

    if(*extrap_dir_bcs > 0)
    {

      AMREX_FOR_4D(bx_xz_lo_2D, 3, i, j, k, n,
      {
        const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

        if(bct == minf)
          vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i,j+1,k,n);
      });
    }
  }

  if (ntop > 0)
  {
    AMREX_FOR_4D(bx_xz_hi_3D, 3, i, j, k, n,
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == pinf) or (bct == pout))
        vel(i,j,k,n) = vel(i,dom_hi[1],k,n);
      else if(bct == minf)
      {
        if(n == 1)
          vel(i,j,k,n) = p_bc_v_g[bcv];
        else
          vel(i,j,k,n) = 0;
      }
    });

    if(*extrap_dir_bcs > 0)
    {

      AMREX_FOR_4D(bx_xz_hi_2D, 3, i, j, k, n,
      {
        const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

        if(bct == minf)
          vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i,j-1,k,n);
      });
    }
  }

  if (ndwn > 0)
  {
    AMREX_FOR_4D(bx_xy_lo_3D, 3, i, j, k, n,
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == pinf) or (bct == pout))
        vel(i,j,k,n) = vel(i,j,dom_lo[2],n);
      else if(bct == minf)
      {
        if(n == 2)
          vel(i,j,k,n) = p_bc_w_g[bcv];
        else
          vel(i,j,k,n) = 0;
      }
    });

    if(*extrap_dir_bcs > 0)
    {

      AMREX_FOR_4D(bx_xy_lo_2D, 3, i, j, k, n,
      {
        const int bct = bct_klo(i,j,dom_lo[2]-1,0);

        if(bct == minf)
          vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i,j,k+1,n);
      });
    }
  }

  if (nup > 0)
  {
    AMREX_FOR_4D(bx_xy_hi_3D, 3, i, j, k, n,
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == pinf) or (bct == pout))
        vel(i,j,k,n) = vel(i,j,dom_hi[2],n);
      else if(bct == minf)
      {
        if(n == 2)
          vel(i,j,k,n) = p_bc_w_g[bcv];
        else
          vel(i,j,k,n) = 0;
      }
    });

    if(*extrap_dir_bcs > 0)
    {
      AMREX_FOR_4D(bx_xy_hi_2D, 3, i, j, k, n,
      {
        const int bct = bct_khi(i,j,dom_hi[2]+1,0);

        if(bct == minf)
          vel(i,j,k,n) = 2*vel(i,j,k,n) - vel(i,j,k-1,n);
      });
    }
  }
}


void
mfix::set_vec_bcs(const int lev,
                  FArrayBox& vec_fab,
                  const Box& domain) const
{
  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<Real> const& vec = vec_fab.array();

  IntVect vec_lo(vec_fab.loVect());
  IntVect vec_hi(vec_fab.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  const int nlft = std::max(0, dom_lo[0]-vec_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-vec_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-vec_lo[2]);

  const int nrgt = std::max(0, vec_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, vec_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, vec_hi[2]-dom_hi[2]);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(vec_hi), bx_xz_lo_hi_3D(vec_hi), bx_xy_lo_hi_3D(vec_hi);
  IntVect bx_yz_hi_lo_3D(vec_lo), bx_xz_hi_lo_3D(vec_lo), bx_xy_hi_lo_3D(vec_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(vec_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, vec_hi);

  const Box bx_xz_lo_3D(vec_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, vec_hi);

  const Box bx_xy_lo_3D(vec_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, vec_hi);

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  const amrex::Real* p_bc_u_g = m_bc_u_g.data();
  const amrex::Real* p_bc_v_g = m_bc_v_g.data();
  const amrex::Real* p_bc_w_g = m_bc_w_g.data();

  amrex::Real* p_bc_ep_g = m_bc_ep_g.data();

  if (nlft > 0)
  {
    AMREX_FOR_4D(bx_yz_lo_3D, 3, i, j, k, n,
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == pinf) or (bct == pout))
        vec(i,j,k,n) = vec(dom_lo[0],j,k,n);
      else if(bct == minf)
      {
        if(n == 0)
          vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
        if(n == 1)
          vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
        if(n == 2)
          vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }

  if (nrgt > 0)
  {
    AMREX_FOR_4D(bx_yz_hi_3D, 3, i, j, k, n,
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == pinf) or (bct == pout))
        vec(i,j,k,n) = vec(dom_hi[0],j,k,n);
      else if(bct == minf)
      {
        if(n == 0)
          vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
        if(n == 1)
          vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
        if(n == 2)
          vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }

  if (nbot > 0)
  {
    AMREX_FOR_4D(bx_xz_lo_3D, 3, i, j, k, n,
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == pinf) or (bct == pout))
        vec(i,j,k,n) = vec(i,dom_lo[1],k,n);
      else if(bct == minf)
      {
        if(n == 0)
          vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
        if(n == 1)
          vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
        if(n == 2)
          vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }

  if (ntop > 0)
  {
    AMREX_FOR_4D(bx_xz_hi_3D, 3, i, j, k, n,
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == pinf) or (bct == pout))
        vec(i,j,k,n) = vec(i,dom_hi[1],k,n);
      else if(bct == minf)
      {
        if(n == 0)
          vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
        if(n == 1)
          vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
        if(n == 2)
          vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }

  if (ndwn > 0)
  {
    AMREX_FOR_4D(bx_xy_lo_3D, 3, i, j, k, n,
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == pinf) or (bct == pout))
        vec(i,j,k,n) = vec(i,j,dom_lo[2],n);
      else if(bct == minf)
      {
        if(n == 0)
          vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
        if(n == 1)
          vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
        if(n == 2)
          vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }

  if (nup > 0)
  {
    AMREX_FOR_4D(bx_xy_hi_3D, 3, i, j, k, n,
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == pinf) or (bct == pout))
        vec(i,j,k,n) = vec(i,j,dom_hi[2],n);
      else if(bct == minf)
      {
        if(n == 0)
          vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
        if(n == 1)
          vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
        if(n == 2)
          vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }
}
