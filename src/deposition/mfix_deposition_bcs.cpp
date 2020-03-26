#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix.H>
#include <param_mod_F.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>

void mfix::mfix_deposition_bcs_scalar (int lev, amrex::MultiFab& filled_mf)
{
  BL_PROFILE("mfix::mfix_deposition_bcs_scalar()");

  Box domain(Geom(lev).Domain());

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();

  for (MFIter mfi(filled_mf, false); mfi.isValid(); ++mfi) {

    const Box& sbx = filled_mf[mfi].box();

    // Extract the lower and upper boundaries of Box sbx and Domain
    const IntVect sbx_lo(sbx.loVect()), sbx_hi(sbx.hiVect());
    const amrex::Dim3& dom_lo = amrex::lbound(domain);
    const amrex::Dim3& dom_hi = amrex::ubound(domain);

    // Create a 2D Box collapsing sbx on x-direction
    IntVect sbx_yz_hi(sbx.hiVect());
    sbx_yz_hi[0] = sbx_lo[0];
    const Box sbx_yz(sbx_lo, sbx_yz_hi);

    // Create a 2D Box collapsing sbx on y-direction
    IntVect sbx_xz_hi(sbx.hiVect());
    sbx_xz_hi[1] = sbx_lo[1];
    const Box sbx_xz(sbx_lo, sbx_xz_hi);

    // Create a 2D Box collapsing sbx on z-direction
    IntVect sbx_xy_hi(sbx.hiVect());
    sbx_xy_hi[2] = sbx_lo[2];
    const Box sbx_xy(sbx_lo, sbx_xy_hi);

    Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
    Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();
    Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
    Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();
    Array4<int> const& bc_klo_type = bc_klo[lev]->array();
    Array4<int> const& bc_khi_type = bc_khi[lev]->array();

    Array4<Real> const& vol = filled_mf.array(mfi);

    if(sbx_lo[0] < dom_lo.x) {
      const int ilo = dom_lo.x;

      amrex::ParallelFor(sbx_yz, [bc_ilo_type,pinf,minf,vol,ilo]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct_ilo = bc_ilo_type(ilo-1,j,k,0);

          if(bct_ilo == pinf or bct_ilo == minf) 
          {
            vol(ilo,j,k) += vol(ilo-1,j,k);
            vol(ilo-1,j,k) = 0;
          }
        });
    }

    if(sbx_hi[0] > dom_hi.x) {
      const int ihi = dom_hi.x;

      amrex::ParallelFor(sbx_yz, [bc_ihi_type,pinf,minf,vol,ihi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct_ihi = bc_ihi_type(ihi+1,j,k,0);

          if(bct_ihi == pinf or bct_ihi == minf)
          {
            vol(ihi,j,k) += vol(ihi+1,j,k);
            vol(ihi+1,j,k) = 0;
          }
        });
    }

    if(sbx_lo[1] < dom_lo.y) {
      const int jlo = dom_lo.y;

      amrex::ParallelFor(sbx_xz, [bc_jlo_type,pinf,minf,vol,jlo]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct_jlo = bc_jlo_type(i,jlo-1,k,0);

          if(bct_jlo == pinf or bct_jlo == minf)
          {
            vol(i,jlo,k) += vol(i,jlo-1,k);
            vol(i,jlo-1,k) = 0;
          }
        });
    }

    if(sbx_hi[1] > dom_hi.y)
      {
        const int jhi = dom_hi.y;

      amrex::ParallelFor(sbx_xz, [bc_jhi_type,pinf,minf,vol,jhi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct_jhi = bc_jhi_type(i,jhi+1,k,0);

          if(bct_jhi == pinf or bct_jhi == minf)
          {
            vol(i,jhi,k) += vol(i,jhi+1,k);
            vol(i,jhi+1,k) = 0;
          }
        });
      }

    if(sbx_lo[2] < dom_lo.z) {
      const int klo = dom_lo.z;
      
      amrex::ParallelFor(sbx_xy, [bc_klo_type,pinf,minf,vol,klo]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct_klo = bc_klo_type(i,j,klo-1,0);

          if(bct_klo == pinf or bct_klo == minf)
          {
            vol(i,j,klo) += vol(i,j,klo-1);
            vol(i,j,klo-1) = 0;
          }
        });
    }

    if(sbx_hi[2] > dom_hi.z) {
      const int khi = dom_hi.z;
      
      amrex::ParallelFor(sbx_xy, [bc_khi_type,pinf,minf,vol,khi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          const int bct_khi = bc_khi_type(i,j,khi+1,0);

          if(bct_khi == pinf or bct_khi == minf)
          {
            vol(i,j,khi) += vol(i,j,khi+1);
            vol(i,j,khi+1) = 0;
          }
        });
    }

  }
}
