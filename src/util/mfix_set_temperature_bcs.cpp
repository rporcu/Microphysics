#include <mfix.H>

#include <param_mod_F.H>

#include <MFIX_FLUID_Parms.H>

using namespace amrex;

//
// Set the BCs for temperature only
//
void
mfix::mfix_set_temperature_bcs (Real time,
                                Vector< MultiFab* > const& T_g_in)
{
  BL_PROFILE("mfix::mfix_set_temperature_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     // Set all values outside the domain to covered_val just to avoid use of undefined
     T_g_in[lev]->setDomainBndry(covered_val,geom[lev]);

     T_g_in[lev]->FillBoundary(geom[lev].periodicity());

     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*T_g_in[lev], false); mfi.isValid(); ++mfi)
       set_temperature_bcs(time, lev, (*T_g_in[lev])[mfi], domain);

     EB_set_covered(*T_g_in[lev], 0, T_g_in[lev]->nComp(), T_g_in[lev]->nGrow(), covered_val);

     // Do this after as well as before to pick up terms that got updated in the
     // call above
     T_g_in[lev]->FillBoundary(geom[lev].periodicity());
  }
}

void
mfix::set_temperature_bcs (Real time,
                           const int lev,
                           FArrayBox& T_g_fab,
                           const Box& domain)
{
  BL_PROFILE("mfix::set_temperature_bcs()");

  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  Array4<Real> const& T_g = T_g_fab.array();

  IntVect T_g_lo(T_g_fab.loVect());
  IntVect T_g_hi(T_g_fab.hiVect());

  const int nlft = std::max(0, dom_lo[0]-T_g_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-T_g_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-T_g_lo[2]);

  const int nrgt = std::max(0, T_g_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, T_g_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, T_g_hi[2]-dom_hi[2]);

  // Create InVects for following 2D Boxes
  IntVect bx_yz_lo_lo_2D(T_g_lo), bx_yz_lo_hi_2D(T_g_hi);
  IntVect bx_yz_hi_lo_2D(T_g_lo), bx_yz_hi_hi_2D(T_g_hi);
  IntVect bx_xz_lo_lo_2D(T_g_lo), bx_xz_lo_hi_2D(T_g_hi);
  IntVect bx_xz_hi_lo_2D(T_g_lo), bx_xz_hi_hi_2D(T_g_hi);
  IntVect bx_xy_lo_lo_2D(T_g_lo), bx_xy_lo_hi_2D(T_g_hi);
  IntVect bx_xy_hi_lo_2D(T_g_lo), bx_xy_hi_hi_2D(T_g_hi);

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
  IntVect bx_yz_lo_hi_3D(T_g_hi), bx_xz_lo_hi_3D(T_g_hi), bx_xy_lo_hi_3D(T_g_hi);
  IntVect bx_yz_hi_lo_3D(T_g_lo), bx_xz_hi_lo_3D(T_g_lo), bx_xy_hi_lo_3D(T_g_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(T_g_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, T_g_hi);

  const Box bx_xz_lo_3D(T_g_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, T_g_hi);

  const Box bx_xy_lo_3D(T_g_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, T_g_hi);

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  amrex::Real* p_bc_t_g = m_bc_t_g.data();

  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);

      if((bct == pout)) {
        T_g(i,j,k) = T_g(dom_lo[0],j,k);
      }
      else if ((bct == minf) or (bct == pinf)) {
        T_g(i,j,k) = p_bc_t_g[bcv];
      }
    });

    // BCs extrapolation
    amrex::ParallelFor(bx_yz_lo_2D,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == minf) or (bct == pinf)) {
        T_g(i,j,k) = 2*T_g(i,j,k) - T_g(i+1,j,k);
      }
    });
  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);

      if(bct == pout) {
        T_g(i,j,k) = T_g(dom_hi[0],j,k);
      }
      else if ((bct == minf) or (bct == pinf)) {
        T_g(i,j,k) = p_bc_t_g[bcv];
      }
    });

    // Extrapolate BCs
    amrex::ParallelFor(bx_yz_hi_2D, 
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == minf) or (bct == pinf)) {
        T_g(i,j,k) = 2*T_g(i,j,k) - T_g(i-1,j,k);
      }
    });
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D, 
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);

      if(bct == pout) {
        T_g(i,j,k) = T_g(i,dom_lo[1],k);
        }
      else if ((bct == minf) or (bct == pinf))
        T_g(i,j,k) = p_bc_t_g[bcv];
    });

    amrex::ParallelFor(bx_xz_lo_2D,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if ((bct == minf) or (bct == pinf))
        T_g(i,j,k) = 2*T_g(i,j,k) - T_g(i,j+1,k);
    });
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D, 
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);

      if(bct == pout)
        T_g(i,j,k) = T_g(i,dom_hi[1],k);
      else if ((bct == minf) or (bct == pinf))
        T_g(i,j,k) = p_bc_t_g[bcv];
    });

    amrex::ParallelFor(bx_xz_hi_2D, 
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if ((bct == minf) or (bct == pinf))
        T_g(i,j,k) = 2*T_g(i,j,k) - T_g(i,j-1,k);
    });
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D, 
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);

      if(bct == pout)
        T_g(i,j,k) = T_g(i,j,dom_lo[2]);
      else if ((bct == minf) or (bct == pinf))
        T_g(i,j,k) = p_bc_t_g[bcv];
    });

    amrex::ParallelFor(bx_xy_lo_2D,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if ((bct == minf) or (bct == pinf))
        T_g(i,j,k) = 2*T_g(i,j,k) - T_g(i,j,k+1);
    });
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);

      if(bct == pout)
        T_g(i,j,k) = T_g(i,j,dom_hi[2]);
      else if ((bct == minf) or (bct == pinf))
        T_g(i,j,k) = p_bc_t_g[bcv];
    });

    amrex::ParallelFor(bx_xy_hi_2D,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == minf) or (bct == pinf))
        T_g(i,j,k) = 2*T_g(i,j,k) - T_g(i,j,k-1);
    });

  }
}
