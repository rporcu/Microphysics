#include <mfix.H>


//
// Set the BCs for face-centroid-based velocity components only
//
void
mfix::set_MAC_velocity_bcs (int lev,
                            Vector< MultiFab const*> const& mac_rhs,
                            Vector< MultiFab*      > const& ep_u_mac,
                            Vector< MultiFab*      > const& ep_v_mac,
                            Vector< MultiFab*      > const& ep_w_mac,
                            Real time)
{
  BL_PROFILE("MacProjection::set_MAC_velocity_bcs()");

  ep_u_mac[lev]->FillBoundary(geom[lev].periodicity());
  ep_v_mac[lev]->FillBoundary(geom[lev].periodicity());
  ep_w_mac[lev]->FillBoundary(geom[lev].periodicity());

  Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(*mac_rhs[lev], false); mfi.isValid(); ++mfi)
  {
    const Box& ubx = (*ep_u_mac[lev])[mfi].box();
    IntVect ubx_lo(ubx.loVect());
    IntVect ubx_hi(ubx.hiVect());

    const Box& vbx = (*ep_v_mac[lev])[mfi].box();
    IntVect vbx_lo(vbx.loVect());
    IntVect vbx_hi(vbx.hiVect());

    const Box& wbx = (*ep_w_mac[lev])[mfi].box();
    IntVect wbx_lo(wbx.loVect());
    IntVect wbx_hi(wbx.hiVect());

    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());

    Array4<Real> const& ep_u = ep_u_mac[lev]->array(mfi);
    Array4<Real> const& ep_v = ep_v_mac[lev]->array(mfi);
    Array4<Real> const& ep_w = ep_w_mac[lev]->array(mfi);

    Array4<int> const& bct_ilo = bc_list.bc_ilo[lev]->array();
    Array4<int> const& bct_ihi = bc_list.bc_ihi[lev]->array();
    Array4<int> const& bct_jlo = bc_list.bc_jlo[lev]->array();
    Array4<int> const& bct_jhi = bc_list.bc_jhi[lev]->array();
    Array4<int> const& bct_klo = bc_list.bc_klo[lev]->array();
    Array4<int> const& bct_khi = bc_list.bc_khi[lev]->array();

    const int nlft = amrex::max(0, dom_lo[0]-ubx_lo[0]);
    const int nbot = amrex::max(0, dom_lo[1]-vbx_lo[1]);
    const int ndwn = amrex::max(0, dom_lo[2]-wbx_lo[2]);

    const int nrgt = amrex::max(0, ubx_hi[0]-dom_hi[0]);
    const int ntop = amrex::max(0, vbx_hi[1]-dom_hi[1]);
    const int nup  = amrex::max(0, wbx_hi[2]-dom_hi[2]);

    mfix_usr1(time);

    const int minf = bc_list.get_minf();

    amrex::Real* p_bc_u_g = m_bc_u_g.data();
    amrex::Real* p_bc_v_g = m_bc_v_g.data();
    amrex::Real* p_bc_w_g = m_bc_w_g.data();
    amrex::Real* p_bc_e_g = m_bc_ep_g.data();

    // NOTE - we only call this for MAC velocities which are only defined on normal faces

    if (nlft > 0)
    {
      // Create InVects for following Box
      IntVect ulo_bx_yz_lo(ubx_lo);
      IntVect ulo_bx_yz_hi(ubx_hi);

      // Fix lo and hi limits
      ulo_bx_yz_lo[0] = dom_lo[0];
      ulo_bx_yz_hi[0] = dom_lo[0];

      const Box ulo_bx_yz(ulo_bx_yz_lo, ulo_bx_yz_hi);

      amrex::ParallelFor(ulo_bx_yz,
        [bct_ilo,dom_lo,minf,p_bc_u_g,p_bc_e_g,ep_u]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
        const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
        if(bct == minf) ep_u(i,j,k) = p_bc_u_g[bcv] * p_bc_e_g[bcv];
      });
    }

    if (nrgt > 0)
    {
      // Create InVects for following Box
      IntVect uhi_bx_yz_lo(ubx_lo);
      IntVect uhi_bx_yz_hi(ubx_hi);

      // Fix lo and hi limits
      uhi_bx_yz_lo[0] = dom_hi[0]+1;
      uhi_bx_yz_hi[0] = dom_hi[0]+1;

      const Box uhi_bx_yz(uhi_bx_yz_lo, uhi_bx_yz_hi);

      amrex::ParallelFor(uhi_bx_yz,
        [bct_ihi,dom_hi,minf,p_bc_u_g,p_bc_e_g,ep_u]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
        const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
        if (bct == minf) ep_u(i,j,k) = p_bc_u_g[bcv] * p_bc_e_g[bcv];
      });
    }

    if (nbot > 0)
    {
      // Create InVects for following Box
      IntVect vlo_bx_xz_lo(vbx_lo);
      IntVect vlo_bx_xz_hi(vbx_hi);

      // Fix lo and hi limits
      vlo_bx_xz_lo[1] = dom_lo[1];
      vlo_bx_xz_hi[1] = dom_lo[1];

      const Box vlo_bx_xz(vlo_bx_xz_lo, vlo_bx_xz_hi);

      amrex::ParallelFor(vlo_bx_xz,
        [bct_jlo,dom_lo,minf,p_bc_v_g,p_bc_e_g,ep_v]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
        const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
        if (bct == minf) ep_v(i,j,k) = p_bc_v_g[bcv] * p_bc_e_g[bcv];
      });
    }

    if (ntop > 0)
    {
      // Create InVects for following Box
      IntVect vhi_bx_xz_lo(vbx_lo);
      IntVect vhi_bx_xz_hi(vbx_hi);

      // Fix lo and hi limits
      vhi_bx_xz_lo[1] = dom_hi[1]+1;
      vhi_bx_xz_hi[1] = dom_hi[1]+1;

      const Box vhi_bx_xz(vhi_bx_xz_lo, vhi_bx_xz_hi);

      amrex::ParallelFor(vhi_bx_xz,
        [bct_jhi,dom_hi,minf,p_bc_v_g,p_bc_e_g,ep_v]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
        const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
        if (bct == minf) ep_v(i,j,k) = p_bc_v_g[bcv] * p_bc_e_g[bcv];
      });
    }

    if (ndwn > 0)
    {
      // Create InVects for following Boxes
      IntVect wlo_bx_xy_lo(wbx_lo);
      IntVect wlo_bx_xy_hi(wbx_hi);

      // Fix lo and hi limits
      wlo_bx_xy_lo[2] = dom_lo[2];
      wlo_bx_xy_hi[2] = dom_lo[2];

      const Box wlo_bx_xy(wlo_bx_xy_lo, wlo_bx_xy_hi);

      amrex::ParallelFor(wlo_bx_xy,
        [bct_klo,dom_lo,minf,p_bc_w_g,p_bc_e_g,ep_w]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
        const int bct = bct_klo(i,j,dom_lo[2]-1,0);
        if (bct == minf) ep_w(i,j,k) = p_bc_w_g[bcv] * p_bc_e_g[bcv];
      });
    }

    if (nup > 0)
    {
      // Create InVects for following Boxes
      IntVect whi_bx_xy_lo(wbx_lo);
      IntVect whi_bx_xy_hi(wbx_hi);

      // Fix lo and hi limits
      whi_bx_xy_lo[2] = dom_hi[2]+1;
      whi_bx_xy_hi[2] = dom_hi[2]+1;

      const Box whi_bx_xy(whi_bx_xy_lo, whi_bx_xy_hi);

      amrex::ParallelFor(whi_bx_xy,
        [bct_khi,dom_hi,minf,p_bc_w_g,p_bc_e_g,ep_w]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
        const int bct = bct_khi(i,j,dom_hi[2]+1,0);
        if (bct == minf) ep_w(i,j,k) = p_bc_w_g[bcv] * p_bc_e_g[bcv];
      });
    }
  }
}
