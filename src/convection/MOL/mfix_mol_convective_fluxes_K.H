#ifndef MFIX_MOL_CONVECTIVE_FLUXES_K_H_
#define MFIX_MOL_CONVECTIVE_FLUXES_K_H_

#include <mfix.H>

#include <mfix_upwind.H>
#include <mfix_algorithm.H>

using namespace amrex;
using namespace aux;

namespace mol {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mfix_compute_fluxes_on_box (const int lev, Box& bx,
                                 Box const& ubx,
                                 Box const& vbx,
                                 Box const& wbx,
                                 FArrayBox& a_fx,
                                 FArrayBox& a_fy,
                                 FArrayBox& a_fz,
                                 const FArrayBox& state_in,
                                 const int state_comp, const int ncomp,
                                 const FArrayBox& xslopes_in,
                                 const FArrayBox& yslopes_in,
                                 const FArrayBox& zslopes_in,
                                 const int slopes_comp,
                                 const FArrayBox& ep_u_mac,
                                 const FArrayBox& ep_v_mac,
                                 const FArrayBox& ep_w_mac,
                                 const GpuArray<int, 3> bc_types,
                                 Array4<int const> const& bct_ilo,
                                 Array4<int const> const& bct_ihi,
                                 Array4<int const> const& bct_jlo,
                                 Array4<int const> const& bct_jhi,
                                 Array4<int const> const& bct_klo,
                                 Array4<int const> const& bct_khi,
                                 Vector<Geometry> geom)
{
  Box domain(geom[lev].Domain());

  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<const Real> const& state = state_in.array();

  Array4<Real> const& fx = a_fx.array();
  Array4<Real> const& fy = a_fy.array();
  Array4<Real> const& fz = a_fz.array();

  Array4<const Real> const& x_slopes = xslopes_in.array();
  Array4<const Real> const& y_slopes = yslopes_in.array();
  Array4<const Real> const& z_slopes = zslopes_in.array();

  Array4<const Real> const& u = ep_u_mac.array();
  Array4<const Real> const& v = ep_v_mac.array();
  Array4<const Real> const& w = ep_w_mac.array();


  // Vectorize the boundary conditions list in order to use it in lambda
  // functions
  // const GpuArray<int, 3> bc_types =
  //   {bc_list.get_minf(), bc_list.get_pinf(), bc_list.get_pout()};

  const int ubx_npoints = ubx.numPts();
  const auto ubx_lo = amrex::lbound(ubx);
  const auto ubx_len = amrex::length(ubx);

  const int vbx_npoints = vbx.numPts();
  const auto vbx_lo = amrex::lbound(vbx);
  const auto vbx_len = amrex::length(vbx);

  const int wbx_npoints = wbx.numPts();
  const auto wbx_lo = amrex::lbound(wbx);
  const auto wbx_len = amrex::length(wbx);

  const int npoints = amrex::max(ubx_npoints,vbx_npoints,wbx_npoints);

  ParallelFor(npoints, [bc_types,ubx_npoints,vbx_npoints,wbx_npoints,ubx_len,
      vbx_len,wbx_len,ubx_lo,vbx_lo,wbx_lo,bct_ilo,bct_ihi,bct_jlo,bct_jhi,
      bct_klo,bct_khi,dom_low,dom_high,fx,fy,fz,slopes_comp,x_slopes,y_slopes,
      z_slopes,state_comp,state,ncomp,u,v,w]
    AMREX_GPU_DEVICE (int idx) noexcept
  {
    const int* bct_data = bc_types.data();
    const int bct_size = bc_types.size();

    if(idx < ubx_npoints)
    {
      int k = idx / (ubx_len.x*ubx_len.y);
      int j = (idx - k*(ubx_len.x*ubx_len.y)) / (ubx_len.x);
      int i = (idx - k*(ubx_len.x*ubx_len.y)) - j*ubx_len.x;

      i += ubx_lo.x;
      j += ubx_lo.y;
      k += ubx_lo.z;

      const Real u_val = u(i,j,k);

      const int bct_ilo_val = bct_ilo(dom_low.x-1,j,k,0);
      const int bct_ihi_val = bct_ihi(dom_high.x+1,j,k,0);

      for(int n(0); n < ncomp; n++) {
        //
        // West face
        //
        // In the case of MINF       we are using the prescribed Dirichlet value
        // In the case of PINF, POUT we are using the upwind value
        Real state_pls = state(i,j,k,state_comp+n);
        Real state_mns = state(i-1,j,k,state_comp+n);
        const Real x_slopes_pls = x_slopes(i,j,k,slopes_comp+n);
        const Real x_slopes_mns = x_slopes(i-1,j,k,slopes_comp+n);

        Real fx_val(0);

        if ((i == dom_low.x) and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_ilo_val)))
        {
          fx_val = u_val * state_mns;
        }
        else if ((i == dom_high.x+1) and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_ihi_val)))
        {
          fx_val = u_val * state_pls;
        }
        else {
          state_pls -= .5*x_slopes_pls;
          state_mns += .5*x_slopes_mns;

          fx_val = u_val * upwind(state_mns, state_pls, u_val);
        }

        fx(i,j,k,n) = fx_val;
      }
    }

    if(idx < vbx_npoints)
    {
      int k = idx / (vbx_len.x*vbx_len.y);
      int j = (idx - k*(vbx_len.x*vbx_len.y)) / (vbx_len.x);
      int i = (idx - k*(vbx_len.x*vbx_len.y)) - j*vbx_len.x;

      i += vbx_lo.x;
      j += vbx_lo.y;
      k += vbx_lo.z;

      const Real v_val = v(i,j,k);

      const int bct_jlo_val = bct_jlo(i,dom_low.y-1,k,0);
      const int bct_jhi_val = bct_jhi(i,dom_high.y+1,k,0);

      for(int n(0); n < ncomp; n++) {
        //
        // South face
        //
        // In the case of MINF       we are using the prescribed Dirichlet value
        // In the case of PINF, POUT we are using the upwind value
        Real state_pls = state(i,j,k,state_comp+n);
        Real state_mns = state(i,j-1,k,state_comp+n);
        const Real y_slopes_pls = y_slopes(i,j,k,slopes_comp+n);
        const Real y_slopes_mns = y_slopes(i,j-1,k,slopes_comp+n);

        Real fy_val(0);

        if((j == dom_low.y) and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_jlo_val)))
        {
          fy_val = v_val * state_mns;
        }
        else if ((j == dom_high.y+1) and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_jhi_val)))
        {
          fy_val = v_val * state_pls;
        }
        else {
          state_pls -= .5*y_slopes_pls;
          state_mns += .5*y_slopes_mns;

          fy_val = v_val * upwind(state_mns, state_pls, v_val);
        }

        fy(i,j,k,n) = fy_val;
      }
    }

    if(idx < wbx_npoints)
    {
      int k = idx / (wbx_len.x*wbx_len.y);
      int j = (idx - k*(wbx_len.x*wbx_len.y)) / (wbx_len.x);
      int i = (idx - k*(wbx_len.x*wbx_len.y)) - j*wbx_len.x;

      i += wbx_lo.x;
      j += wbx_lo.y;
      k += wbx_lo.z;

      const Real w_val = w(i,j,k);

      const int bct_klo_val = bct_klo(i,j,dom_low.z-1,0);
      const int bct_khi_val = bct_khi(i,j,dom_high.z+1,0);

      for(int n(0); n < ncomp; n++) {
        //
        // Bottom face
        //
        // In the case of MINF       we are using the prescribed Dirichlet value
        // In the case of PINF, POUT we are using the upwind value
        Real state_pls = state(i,j,k,state_comp+n);
        Real state_mns = state(i,j,k-1,state_comp+n);
        const Real z_slopes_pls = z_slopes(i,j,k,slopes_comp+n);
        const Real z_slopes_mns = z_slopes(i,j,k-1,slopes_comp+n);

        Real fz_val(0);

        if((k == dom_low.z) and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_klo_val)))
        {
          fz_val = w_val * state_mns;
        }
        else if ((k == dom_high.z+1) and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_khi_val)))
        {
          fz_val = w_val * state_pls;
        }
        else {
          state_pls -= .5*z_slopes_pls;
          state_mns += .5*z_slopes_mns;

          fz_val = w_val * upwind(state_mns, state_pls, w_val);
        }

        fz(i,j,k,n) = fz_val;
      }
    }
  });
}



}//end mol namespace

#endif
