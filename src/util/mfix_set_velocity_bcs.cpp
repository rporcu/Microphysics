#include <mfix.H>

//
//  These subroutines set the BCs for the velocity components only.
//

void
MFIXBoundaryConditions::set_velocity_bcs (Real time,
                                          Vector< MultiFab* > const& vel_g_in,
                                          int extrap_dir_bcs)
{
  BL_PROFILE("MFIXBoundaryConditions::set_velocity_bcs()");

  const int nlev = vel_g_in.size();

  for (int lev = 0; lev < nlev; lev++) {

    // Set all values outside the domain to covered_val just to avoid use of undefined
    vel_g_in[lev]->setDomainBndry(mfix::covered_val, m_geom[lev]);

    Box domain(m_geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vel_g_in[lev], false); mfi.isValid(); ++mfi) {

      FArrayBox& vel_fab = (*vel_g_in[lev])[mfi];

      IntVect dom_lo(domain.loVect());
      IntVect dom_hi(domain.hiVect());

      Array4<Real> const& vel = vel_fab.array();

      IntVect vel_lo(vel_fab.loVect());
      IntVect vel_hi(vel_fab.hiVect());

      Array4<const int> const& bct_ilo = m_bc_list.bc_ilo[lev]->array();
      Array4<const int> const& bct_ihi = m_bc_list.bc_ihi[lev]->array();
      Array4<const int> const& bct_jlo = m_bc_list.bc_jlo[lev]->array();
      Array4<const int> const& bct_jhi = m_bc_list.bc_jhi[lev]->array();
      Array4<const int> const& bct_klo = m_bc_list.bc_klo[lev]->array();
      Array4<const int> const& bct_khi = m_bc_list.bc_khi[lev]->array();

      const int nlft = amrex::max(0, dom_lo[0]-vel_lo[0]);
      const int nbot = amrex::max(0, dom_lo[1]-vel_lo[1]);
      const int ndwn = amrex::max(0, dom_lo[2]-vel_lo[2]);

      const int nrgt = amrex::max(0, vel_hi[0]-dom_hi[0]);
      const int ntop = amrex::max(0, vel_hi[1]-dom_hi[1]);
      const int nup  = amrex::max(0, vel_hi[2]-dom_hi[2]);

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

      // Create 2D boxes for GPU loops
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

      // Create 3D boxes for GPU loops
      const Box bx_yz_lo_3D(vel_lo, bx_yz_lo_hi_3D);
      const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, vel_hi);

      const Box bx_xz_lo_3D(vel_lo, bx_xz_lo_hi_3D);
      const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, vel_hi);

      const Box bx_xy_lo_3D(vel_lo, bx_xy_lo_hi_3D);
      const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, vel_hi);

      set_velocity_bc_values(time);

      mfix::mfix_usr1(time);

      const Real* p_bc_u_g = m_bc_u_g.data();
      const Real* p_bc_v_g = m_bc_v_g.data();
      const Real* p_bc_w_g = m_bc_w_g.data();

      if (nlft > 0)
      {
        amrex::ParallelFor(bx_yz_lo_3D,
          [bct_ilo,dom_lo,p_bc_u_g,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
          const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

          if((bct == BCList::pinf) || (bct == BCList::pout))
          {
             vel(i,j,k,0) = vel(dom_lo[0],j,k,0);
             vel(i,j,k,1) = vel(dom_lo[0],j,k,1);
             vel(i,j,k,2) = vel(dom_lo[0],j,k,2);
          }
          else if(bct == BCList::minf)
          {
             vel(i,j,k,0) = p_bc_u_g[bcv];
             vel(i,j,k,1) = 0;
             vel(i,j,k,2) = 0;
          }
        });

        if(extrap_dir_bcs > 0)
        {
          amrex::ParallelFor(bx_yz_lo_2D,
            [bct_ilo,dom_lo,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

            if(bct == BCList::minf)
            {
               vel(i,j,k,0) = 2*vel(i,j,k,0) - vel(i+1,j,k,0);
               vel(i,j,k,1) = 2*vel(i,j,k,1) - vel(i+1,j,k,1);
               vel(i,j,k,2) = 2*vel(i,j,k,2) - vel(i+1,j,k,2);
            }
          });
        }
      }

      if (nrgt > 0)
      {
        amrex::ParallelFor(bx_yz_hi_3D,
          [bct_ihi,dom_hi,p_bc_u_g,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
          const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

          if((bct == BCList::pinf) || (bct == BCList::pout))
          {
             vel(i,j,k,0) = vel(dom_hi[0],j,k,0);
             vel(i,j,k,1) = vel(dom_hi[0],j,k,1);
             vel(i,j,k,2) = vel(dom_hi[0],j,k,2);
          }
          else if(bct == BCList::minf)
          {
             vel(i,j,k,0) = p_bc_u_g[bcv];
             vel(i,j,k,1) = 0;
             vel(i,j,k,2) = 0;
          }
        });

        if(extrap_dir_bcs > 0)
        {
          amrex::ParallelFor(bx_yz_hi_2D,
            [bct_ihi,dom_hi,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

            if(bct == BCList::minf)
            {
               vel(i,j,k,0) = 2*vel(i,j,k,0) - vel(i-1,j,k,0);
               vel(i,j,k,1) = 2*vel(i,j,k,1) - vel(i-1,j,k,1);
               vel(i,j,k,2) = 2*vel(i,j,k,2) - vel(i-1,j,k,2);
            }
          });
        }
      }

      if (nbot > 0)
      {
        amrex::ParallelFor(bx_xz_lo_3D,
          [bct_jlo,dom_lo,p_bc_v_g,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
          const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

          if((bct == BCList::pinf) || (bct == BCList::pout))
          {
             vel(i,j,k,0) = vel(i,dom_lo[1],k,0);
             vel(i,j,k,1) = vel(i,dom_lo[1],k,1);
             vel(i,j,k,2) = vel(i,dom_lo[1],k,2);
          }
          else if(bct == BCList::minf)
          {
             vel(i,j,k,0) = 0;
             vel(i,j,k,1) = p_bc_v_g[bcv];
             vel(i,j,k,2) = 0;
          }
        });

        if(extrap_dir_bcs > 0)
        {

          amrex::ParallelFor(bx_xz_lo_2D,
            [bct_jlo,dom_lo,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

            if(bct == BCList::minf)
            {
               vel(i,j,k,0) = 2*vel(i,j,k,0) - vel(i,j+1,k,0);
               vel(i,j,k,1) = 2*vel(i,j,k,1) - vel(i,j+1,k,1);
               vel(i,j,k,2) = 2*vel(i,j,k,2) - vel(i,j+1,k,2);
            }
          });
        }
      }

      if (ntop > 0)
      {
        amrex::ParallelFor(bx_xz_hi_3D,
          [bct_jhi,dom_hi,p_bc_v_g,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
          const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

          if((bct == BCList::pinf) || (bct == BCList::pout))
          {
             vel(i,j,k,0) = vel(i,dom_hi[1],k,0);
             vel(i,j,k,1) = vel(i,dom_hi[1],k,1);
             vel(i,j,k,2) = vel(i,dom_hi[1],k,2);
          }
          else if(bct == BCList::minf)
          {
             vel(i,j,k,0) = 0;
             vel(i,j,k,1) = p_bc_v_g[bcv];
             vel(i,j,k,2) = 0;
          }
        });

        if(extrap_dir_bcs > 0)
        {

          amrex::ParallelFor(bx_xz_hi_2D,
            [bct_jhi,dom_hi,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

            if(bct == BCList::minf)
            {
               vel(i,j,k,0) = 2*vel(i,j,k,0) - vel(i,j-1,k,0);
               vel(i,j,k,1) = 2*vel(i,j,k,1) - vel(i,j-1,k,1);
               vel(i,j,k,2) = 2*vel(i,j,k,2) - vel(i,j-1,k,2);
            }
          });
        }
      }

      if (ndwn > 0)
      {
        amrex::ParallelFor(bx_xy_lo_3D,
          [bct_klo,dom_lo,p_bc_w_g,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
          const int bct = bct_klo(i,j,dom_lo[2]-1,0);

          if((bct == BCList::pinf) || (bct == BCList::pout))
          {
             vel(i,j,k,0) = vel(i,j,dom_lo[2],0);
             vel(i,j,k,1) = vel(i,j,dom_lo[2],1);
             vel(i,j,k,2) = vel(i,j,dom_lo[2],2);
          }
          else if(bct == BCList::minf)
          {
              vel(i,j,k,0) = 0;
              vel(i,j,k,1) = 0;
              vel(i,j,k,2) = p_bc_w_g[bcv];
          }
        });

        if(extrap_dir_bcs > 0)
        {

          amrex::ParallelFor(bx_xy_lo_2D,
            [bct_klo,dom_lo,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int bct = bct_klo(i,j,dom_lo[2]-1,0);

            if(bct == BCList::minf)
            {
               vel(i,j,k,0) = 2*vel(i,j,k,0) - vel(i,j,k+1,0);
               vel(i,j,k,1) = 2*vel(i,j,k,1) - vel(i,j,k+1,1);
               vel(i,j,k,2) = 2*vel(i,j,k,2) - vel(i,j,k+1,2);
            }
          });
        }
      }

      if (nup > 0)
      {
        amrex::ParallelFor(bx_xy_hi_3D,
          [bct_khi,dom_hi,p_bc_w_g,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
          const int bct = bct_khi(i,j,dom_hi[2]+1,0);

          if((bct == BCList::pinf) || (bct == BCList::pout))
          {
             vel(i,j,k,0) = vel(i,j,dom_hi[2],0);
             vel(i,j,k,1) = vel(i,j,dom_hi[2],1);
             vel(i,j,k,2) = vel(i,j,dom_hi[2],2);
          }
          else if(bct == BCList::minf)
          {
              vel(i,j,k,0) = 0;
              vel(i,j,k,1) = 0;
              vel(i,j,k,2) = p_bc_w_g[bcv];
          }
        });

        if(extrap_dir_bcs > 0)
        {
          amrex::ParallelFor(bx_xy_hi_2D,
            [bct_khi,dom_hi,vel] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int bct = bct_khi(i,j,dom_hi[2]+1,0);

            if(bct == BCList::minf)
            {
               vel(i,j,k,0) = 2*vel(i,j,k,0) - vel(i,j,k-1,0);
               vel(i,j,k,1) = 2*vel(i,j,k,1) - vel(i,j,k-1,1);
               vel(i,j,k,2) = 2*vel(i,j,k,2) - vel(i,j,k-1,2);
            }
          });
        }
      }

    } // end MFIter loop

    EB_set_covered(*vel_g_in[lev], 0, vel_g_in[lev]->nComp(),
        vel_g_in[lev]->nGrow(), mfix::covered_val);

    // Do this after as well as before to pick up terms that got updated in the call above
    vel_g_in[lev]->FillBoundary(m_geom[lev].periodicity());
  }
}


void
MFIXBoundaryConditions::set_vec_bcs (const int lev,
                                     FArrayBox& vec_fab,
                                     const Box& domain)
{
  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<Real> const& vec = vec_fab.array();

  IntVect vec_lo(vec_fab.loVect());
  IntVect vec_hi(vec_fab.hiVect());

  Array4<const int> const& bct_ilo = m_bc_list.bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = m_bc_list.bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = m_bc_list.bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = m_bc_list.bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = m_bc_list.bc_klo[lev]->array();
  Array4<const int> const& bct_khi = m_bc_list.bc_khi[lev]->array();

  const int nlft = amrex::max(0, dom_lo[0]-vec_lo[0]);
  const int nbot = amrex::max(0, dom_lo[1]-vec_lo[1]);
  const int ndwn = amrex::max(0, dom_lo[2]-vec_lo[2]);

  const int nrgt = amrex::max(0, vec_hi[0]-dom_hi[0]);
  const int ntop = amrex::max(0, vec_hi[1]-dom_hi[1]);
  const int nup  = amrex::max(0, vec_hi[2]-dom_hi[2]);

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

  // Create 3D boxes for GPU loops
  const Box bx_yz_lo_3D(vec_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, vec_hi);

  const Box bx_xz_lo_3D(vec_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, vec_hi);

  const Box bx_xy_lo_3D(vec_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, vec_hi);

  const Real* p_bc_u_g = m_bc_u_g.data();
  const Real* p_bc_v_g = m_bc_v_g.data();
  const Real* p_bc_w_g = m_bc_w_g.data();

  const Real* p_bc_ep_g = m_bc_ep_g.data();

  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D,
      [bct_ilo,dom_lo,p_bc_u_g,p_bc_v_g,p_bc_w_g,p_bc_ep_g,vec]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
         vec(i,j,k,0) = vec(dom_lo[0],j,k,0);
         vec(i,j,k,1) = vec(dom_lo[0],j,k,1);
         vec(i,j,k,2) = vec(dom_lo[0],j,k,2);
      }
      else if(bct == BCList::minf)
      {
         vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
         vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
         vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D,
      [bct_ihi,dom_hi,p_bc_u_g,p_bc_v_g,p_bc_w_g,p_bc_ep_g,vec]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
         vec(i,j,k,0) = vec(dom_hi[0],j,k,0);
         vec(i,j,k,1) = vec(dom_hi[0],j,k,1);
         vec(i,j,k,2) = vec(dom_hi[0],j,k,2);
      }
      else if(bct == BCList::minf)
      {
         vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
         vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
         vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D,
      [bct_jlo,dom_lo,p_bc_u_g,p_bc_v_g,p_bc_w_g,p_bc_ep_g,vec]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
        vec(i,j,k,0) = vec(i,dom_lo[1],k,0);
        vec(i,j,k,1) = vec(i,dom_lo[1],k,1);
        vec(i,j,k,2) = vec(i,dom_lo[1],k,2);
      }
      else if(bct == BCList::minf)
      {
          vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
          vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
          vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D,
      [bct_jhi,dom_hi,p_bc_u_g,p_bc_v_g,p_bc_w_g,p_bc_ep_g,vec]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
         vec(i,j,k,0) = vec(i,dom_hi[1],k,0);
         vec(i,j,k,1) = vec(i,dom_hi[1],k,1);
         vec(i,j,k,2) = vec(i,dom_hi[1],k,2);
      }
      else if(bct == BCList::minf)
      {
         vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
         vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
         vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D,
      [bct_klo,dom_lo,p_bc_u_g,p_bc_v_g,p_bc_w_g,p_bc_ep_g,vec]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
         vec(i,j,k,0) = vec(i,j,dom_lo[2],0);
         vec(i,j,k,1) = vec(i,j,dom_lo[2],1);
         vec(i,j,k,2) = vec(i,j,dom_lo[2],2);
      }
      else if(bct == BCList::minf)
      {
         vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
         vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
         vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D,
      [bct_khi,dom_hi,p_bc_u_g,p_bc_v_g,p_bc_w_g,p_bc_ep_g,vec]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
         vec(i,j,k,0) = vec(i,j,dom_hi[2],0);
         vec(i,j,k,1) = vec(i,j,dom_hi[2],1);
         vec(i,j,k,2) = vec(i,j,dom_hi[2],2);
      }
      else if(bct == BCList::minf)
      {
         vec(i,j,k,0) = p_bc_ep_g[bcv] * p_bc_u_g[bcv];
         vec(i,j,k,1) = p_bc_ep_g[bcv] * p_bc_v_g[bcv];
         vec(i,j,k,2) = p_bc_ep_g[bcv] * p_bc_w_g[bcv];
      }
    });
  }
}
