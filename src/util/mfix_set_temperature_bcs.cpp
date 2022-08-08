#include <mfix.H>
#include <mfix_fluid.H>

using namespace amrex;

//
// Set the BCs for temperature only
//
void
MFIXBoundaryConditions::set_temperature_bcs (Real time,
                                             MFIXFluidPhase& fluid,
                                             Vector< MultiFab* > const& T_g_in)
{
  BL_PROFILE("MFIXBoundaryConditions::set_temperature_bcs()");

  const int nlev = T_g_in.size();

  for (int lev = 0; lev < nlev; lev++) {

    Box domain(m_geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*T_g_in[lev], false); mfi.isValid(); ++mfi) {

      FArrayBox& T_g_fab = (*T_g_in[lev])[mfi];

      IntVect dom_lo(domain.loVect());
      IntVect dom_hi(domain.hiVect());

      Array4<const int> const& bct_ilo = m_bc_list.bc_ilo[lev]->array();
      Array4<const int> const& bct_ihi = m_bc_list.bc_ihi[lev]->array();
      Array4<const int> const& bct_jlo = m_bc_list.bc_jlo[lev]->array();
      Array4<const int> const& bct_jhi = m_bc_list.bc_jhi[lev]->array();
      Array4<const int> const& bct_klo = m_bc_list.bc_klo[lev]->array();
      Array4<const int> const& bct_khi = m_bc_list.bc_khi[lev]->array();

      Array4<Real> const& T_g = T_g_fab.array();

      IntVect T_g_lo(T_g_fab.loVect());
      IntVect T_g_hi(T_g_fab.hiVect());

      const int nlft = amrex::max(0, dom_lo[0]-T_g_lo[0]);
      const int nbot = amrex::max(0, dom_lo[1]-T_g_lo[1]);
      const int ndwn = amrex::max(0, dom_lo[2]-T_g_lo[2]);

      const int nrgt = amrex::max(0, T_g_hi[0]-dom_hi[0]);
      const int ntop = amrex::max(0, T_g_hi[1]-dom_hi[1]);
      const int nup  = amrex::max(0, T_g_hi[2]-dom_hi[2]);

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

      // Create 3D boxes for GPU loops
      const Box bx_yz_lo_3D(T_g_lo, bx_yz_lo_hi_3D);
      const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, T_g_hi);

      const Box bx_xz_lo_3D(T_g_lo, bx_xz_lo_hi_3D);
      const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, T_g_hi);

      const Box bx_xy_lo_3D(T_g_lo, bx_xy_lo_hi_3D);
      const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, T_g_hi);

      set_temperature_bc_values(time, fluid);
      Real* p_bc_t_g = m_bc_t_g.data();

      if (nlft > 0)
      {
        amrex::ParallelFor(bx_yz_lo_3D, [bct_ilo,dom_lo,T_g,p_bc_t_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
          const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);

          if(bct == BCList::pout) {
            T_g(i,j,k) = T_g(dom_lo[0],j,k);
          }
          else if (bct == BCList::minf || bct == BCList::pinf) {
            T_g(i,j,k) = p_bc_t_g[bcv];
          }
        });
      }

      if (nrgt > 0)
      {
        amrex::ParallelFor(bx_yz_hi_3D, [bct_ihi,dom_hi,T_g,p_bc_t_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
          const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);

          if(bct == BCList::pout) {
            T_g(i,j,k) = T_g(dom_hi[0],j,k);
          }
          else if (bct == BCList::minf || bct == BCList::pinf) {
            T_g(i,j,k) = p_bc_t_g[bcv];
          }
        });
      }

      if (nbot > 0)
      {
        amrex::ParallelFor(bx_xz_lo_3D, [bct_jlo,dom_lo,T_g,p_bc_t_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
          const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);

          if(bct == BCList::pout) {
            T_g(i,j,k) = T_g(i,dom_lo[1],k);
          }
          else if (bct == BCList::minf || bct == BCList::pinf) {
            T_g(i,j,k) = p_bc_t_g[bcv];
          }
        });
      }

      if (ntop > 0)
      {
        amrex::ParallelFor(bx_xz_hi_3D, [bct_jhi,dom_hi,T_g,p_bc_t_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
          const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);

          if(bct == BCList::pout) {
            T_g(i,j,k) = T_g(i,dom_hi[1],k);
          }
          else if (bct == BCList::minf || bct == BCList::pinf) {
            T_g(i,j,k) = p_bc_t_g[bcv];
          }
        });
      }

      if (ndwn > 0)
      {
        amrex::ParallelFor(bx_xy_lo_3D, [bct_klo,dom_lo,T_g,p_bc_t_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct = bct_klo(i,j,dom_lo[2]-1,0);
          const int bcv = bct_klo(i,j,dom_lo[2]-1,1);

          if(bct == BCList::pout) {
            T_g(i,j,k) = T_g(i,j,dom_lo[2]);
          }
          else if (bct == BCList::minf || bct == BCList::pinf) {
            T_g(i,j,k) = p_bc_t_g[bcv];
          }
        });
      }

      if (nup > 0)
      {
        amrex::ParallelFor(bx_xy_hi_3D, [bct_khi,dom_hi,T_g,p_bc_t_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct = bct_khi(i,j,dom_hi[2]+1,0);
          const int bcv = bct_khi(i,j,dom_hi[2]+1,1);

          if(bct == BCList::pout) {
            T_g(i,j,k) = T_g(i,j,dom_hi[2]);
          }
          else if (bct == BCList::minf || bct == BCList::pinf) {
            T_g(i,j,k) = p_bc_t_g[bcv];
          }
        });
      }

    } // end MFIter loop

    T_g_in[lev]->FillBoundary(m_geom[lev].periodicity());

    EB_set_covered(*T_g_in[lev], 0, T_g_in[lev]->nComp(), T_g_in[lev]->nGrow(),
        mfix::covered_val);
  }
}
