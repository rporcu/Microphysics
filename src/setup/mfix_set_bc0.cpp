#include <mfix_set_bc0.hpp>
#include <bc_mod_F.H>
#include <eos_mod_F.H>
#include <fld_constants_mod_F.H>
#include <param_mod_F.H>

void set_bc0(const Box& sbx,
             const BcList& bc_list,
             FArrayBox& ep_g_fab,
             FArrayBox& ro_g_fab,
             FArrayBox& mu_g_fab,
             const IArrayBox& bct_ilo_fab,
             const IArrayBox& bct_ihi_fab,
             const IArrayBox& bct_jlo_fab,
             const IArrayBox& bct_jhi_fab,
             const IArrayBox& bct_klo_fab,
             const IArrayBox& bct_khi_fab,
             const Box& domain,
             const int* ng)
{
  Real bc_ro_g, bc_mu_g;
  const Real ro_g0(get_ro_g0());
  const Real mu_g0(get_mu_g0());

  Array4<Real> const& ep_g = ep_g_fab.array();
  Array4<Real> const& ro_g = ro_g_fab.array();
  Array4<Real> const& mu_g = mu_g_fab.array();

  const IntVect sbx_lo(sbx.loVect());
  const IntVect sbx_hi(sbx.hiVect());

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  Array4<const int> const& bct_ilo = bct_ilo_fab.array();
  Array4<const int> const& bct_ihi = bct_ihi_fab.array();
  Array4<const int> const& bct_jlo = bct_jlo_fab.array();
  Array4<const int> const& bct_jhi = bct_jhi_fab.array();
  Array4<const int> const& bct_klo = bct_klo_fab.array();
  Array4<const int> const& bct_khi = bct_khi_fab.array();

  const int nlft = std::max(0,dom_lo[0]-sbx_lo[0]);
  const int nbot = std::max(0,dom_lo[1]-sbx_lo[1]);
  const int ndwn = std::max(0,dom_lo[2]-sbx_lo[2]);

  const int nrgt = std::max(0,sbx_hi[0]-dom_hi[0]);
  const int ntop = std::max(0,sbx_hi[1]-dom_hi[1]);
  const int nup  = std::max(0,sbx_hi[2]-dom_hi[2]);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(sbx_hi), bx_xz_lo_hi_3D(sbx_hi), bx_xy_lo_hi_3D(sbx_hi);
  IntVect bx_yz_hi_lo_3D(sbx_lo), bx_xz_hi_lo_3D(sbx_lo), bx_xy_hi_lo_3D(sbx_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(sbx_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, sbx_hi);

  const Box bx_xz_lo_3D(sbx_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, sbx_hi);

  const Box bx_xy_lo_3D(sbx_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, sbx_hi);

  if (nlft > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_lo_3D, i, j, k,
    {
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        bc_ro_g = ro_g0;

        if (is_undefined_db(mu_g0))
          bc_mu_g = sutherland(get_bc_t_g(bcv));
        else
          bc_mu_g = mu_g0;

        ep_g(i,j,k) = get_bc_ep_g(bcv);
        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });
  }
  
  if (nrgt > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_hi_3D, i, j, k,
    {
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        bc_ro_g = ro_g0;

        if (is_undefined_db(mu_g0))
          bc_mu_g = sutherland(get_bc_t_g(bcv));
        else
          bc_mu_g = mu_g0;

        ep_g(i,j,k) = get_bc_ep_g(bcv);
        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });
  }

  if (nbot > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xz_lo_3D, i, j, k,
    {
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        bc_ro_g = ro_g0;

        if (is_undefined_db(mu_g0))
           bc_mu_g = sutherland(get_bc_t_g(bcv));
        else
           bc_mu_g = mu_g0;

        ep_g(i,j,k) = get_bc_ep_g(bcv);
        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });
  }

  if (ntop > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xz_hi_3D, i, j, k,
    {
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        bc_ro_g = ro_g0;

        if (is_undefined_db(mu_g0))
           bc_mu_g = sutherland(get_bc_t_g(bcv));
        else
           bc_mu_g = mu_g0;

        ep_g(i,j,k) = get_bc_ep_g(bcv);
        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });
  }

  if (ndwn > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_xy_lo_3D, i, j, k,
    {
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        bc_ro_g = ro_g0;

        if (is_undefined_db(mu_g0))
           bc_mu_g = sutherland(get_bc_t_g(bcv));
        else
           bc_mu_g = mu_g0;

        ep_g(i,j,k) = get_bc_ep_g(bcv);
        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });
  }

  if (nup > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(bx_yz_hi_3D, i, j, k,
    {
      const int bcv = bct_ihi(i,j,dom_hi[2]+1,1);
      const int bct = bct_ihi(i,j,dom_hi[2]+1,0);

      if((bct == bc_list.pinf) or (bct == bc_list.pout) or (bct == bc_list.minf))
      {
        bc_ro_g = ro_g0;

        if (is_undefined_db(mu_g0))
           bc_mu_g = sutherland(get_bc_t_g(bcv));
        else
           bc_mu_g = mu_g0;

        ep_g(i,j,k) = get_bc_ep_g(bcv);
        ro_g(i,j,k) = bc_ro_g;
        mu_g(i,j,k) = bc_mu_g;
      }
    });
  }
}
