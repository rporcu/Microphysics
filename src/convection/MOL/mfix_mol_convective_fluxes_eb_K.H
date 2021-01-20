#ifndef MFIX_MOL_CONVECTIVE_FLUXES_EB_K_H_
#define MFIX_MOL_CONVECTIVE_FLUXES_EB_K_H_

#include <mfix.H>
#include <mfix_upwind.H>
#include <mfix_algorithm.H>

using namespace aux;
using namespace amrex;

namespace mol {

//
// Compute the three components of the convection term when we have embedded
// boundaries
//
void mfix_compute_eb_fluxes_on_box (const int lev, Box& bx,
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
                                    const FArrayBox& afrac_x_fab,
                                    const FArrayBox& afrac_y_fab,
                                    const FArrayBox& afrac_z_fab,
                                    const FArrayBox& face_centroid_x,
                                    const FArrayBox& face_centroid_y,
                                    const FArrayBox& face_centroid_z,
                                    const FArrayBox& cell_centroid,
                                    const FArrayBox& volfrac,
                                    const FArrayBox& bndry_centroid,
                                    const EBCellFlagFab& flags,
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

  Array4<Real> const& fx = a_fx.array();
  Array4<Real> const& fy = a_fy.array();
  Array4<Real> const& fz = a_fz.array();

  Array4<const Real> const& state = state_in.array();

  Array4<const Real> const& areafrac_x = afrac_x_fab.array();
  Array4<const Real> const& areafrac_y = afrac_y_fab.array();
  Array4<const Real> const& areafrac_z = afrac_z_fab.array();

  Array4<const Real> const& u = ep_u_mac.array();
  Array4<const Real> const& v = ep_v_mac.array();
  Array4<const Real> const& w = ep_w_mac.array();

  Array4<const Real> const& x_slopes = xslopes_in.array();
  Array4<const Real> const& y_slopes = yslopes_in.array();
  Array4<const Real> const& z_slopes = zslopes_in.array();

  const Box ubx_grown = amrex::surroundingNodes(amrex::grow(bx,1),0);
  const Box vbx_grown = amrex::surroundingNodes(amrex::grow(bx,1),1);
  const Box wbx_grown = amrex::surroundingNodes(amrex::grow(bx,1),2);

  FArrayBox s_on_x_face(ubx_grown, ncomp);
  FArrayBox s_on_y_face(vbx_grown, ncomp);
  FArrayBox s_on_z_face(wbx_grown, ncomp);

  s_on_x_face.prefetchToDevice();
  s_on_y_face.prefetchToDevice();
  s_on_z_face.prefetchToDevice();

  // These lines ensure that the temporary Fabs above aren't destroyed
  //   before we're done with them when running with GPUs
  Elixir eli_x = s_on_x_face.elixir();
  Elixir eli_y = s_on_y_face.elixir();
  Elixir eli_z = s_on_z_face.elixir();

  Array4<Real> const& sx = s_on_x_face.array();
  Array4<Real> const& sy = s_on_y_face.array();
  Array4<Real> const& sz = s_on_z_face.array();

  // Face centroids
  const auto& fcx_fab = face_centroid_x.array();
  const auto& fcy_fab = face_centroid_y.array();
  const auto& fcz_fab = face_centroid_z.array();

  // Cell centroid
  const auto& ccc_fab = cell_centroid.array();

  // const GpuArray<int, 3> bc_types =
  //   {bc_list.get_minf(), bc_list.get_pinf(), bc_list.get_pout()};

  const Real my_huge = 1.e200;

  //
  // First compute the convective fluxes at the face center
  // Do this on ALL faces on the tile, i.e. INCLUDE as many ghost faces as
  // possible
  //

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
      bct_klo,bct_khi,state_comp,dom_low,dom_high,state,x_slopes,y_slopes,z_slopes,
      ccc_fab,fcx_fab,fcy_fab,fcz_fab,areafrac_x,areafrac_y,areafrac_z,u,v,w,sx,
      sy,sz,fx,fy,fz,my_huge,ncomp,slopes_comp]
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

      const Real afrac_x = areafrac_x(i,j,k);

      const Real fcx_fab_x = fcx_fab(i,j,k,0);
      const Real fcx_fab_y = fcx_fab(i,j,k,1);

      const Real ccc_fab_x = ccc_fab(i,j,k,0);
      const Real ccc_fab_y = ccc_fab(i,j,k,1);
      const Real ccc_fab_z = ccc_fab(i,j,k,2);

      const Real ccc_fab_mns_x = ccc_fab(i-1,j,k,0);
      const Real ccc_fab_mns_y = ccc_fab(i-1,j,k,1);
      const Real ccc_fab_mns_z = ccc_fab(i-1,j,k,2);

      for(int n(0); n < ncomp; n++) {
        Real sx_ijkn(0);

        if( afrac_x > 0 )
        {
          if(i <= dom_low.x and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_ilo_val)))
          {
            sx_ijkn = state(dom_low.x-1,j,k,state_comp+n);
          }
          else if(i >= dom_high.x+1 and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_ihi_val)))
          {
            sx_ijkn = state(dom_high.x+1,j,k,state_comp+n);
          }
          else
          {
            Real yf = fcx_fab_x; // local (y,z) of centroid of x-face we are extrapolating to
            Real zf = fcx_fab_y;

            Real delta_x = .5 + ccc_fab_x;
            Real delta_y = yf - ccc_fab_y;
            Real delta_z = zf - ccc_fab_z;

            Real state_pls = state(i,j,k,state_comp+n);
            Real state_mns = state(i-1,j,k,state_comp+n);

            Real cc_umax = amrex::max(state_pls, state_mns);
            Real cc_umin = amrex::min(state_pls, state_mns);

            Real upls = state_pls - delta_x * x_slopes(i,j,k,slopes_comp+n)
                                  + delta_y * y_slopes(i,j,k,slopes_comp+n)
                                  + delta_z * z_slopes(i,j,k,slopes_comp+n);

            upls = amrex::max( amrex::min(upls, cc_umax), cc_umin );

            delta_x = .5 - ccc_fab_mns_x;
            delta_y = yf - ccc_fab_mns_y;
            delta_z = zf - ccc_fab_mns_z;

            Real umns = state_mns + delta_x * x_slopes(i-1,j,k,slopes_comp+n)
                                  + delta_y * y_slopes(i-1,j,k,slopes_comp+n)
                                  + delta_z * z_slopes(i-1,j,k,slopes_comp+n);

            umns = amrex::max( amrex::min(umns, cc_umax), cc_umin );

            sx_ijkn = upwind(umns, upls, u_val);
          }
        }
        else {
          sx_ijkn = my_huge;
        }

        sx(i,j,k,n) = sx_ijkn;
        fx(i,j,k,n) = u_val * sx_ijkn;
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

      const Real afrac_y = areafrac_y(i,j,k);

      const Real fcy_fab_x = fcy_fab(i,j,k,0);
      const Real fcy_fab_y = fcy_fab(i,j,k,1);

      const Real ccc_fab_x = ccc_fab(i,j,k,0);
      const Real ccc_fab_y = ccc_fab(i,j,k,1);
      const Real ccc_fab_z = ccc_fab(i,j,k,2);

      const Real ccc_fab_mns_x = ccc_fab(i,j-1,k,0);
      const Real ccc_fab_mns_y = ccc_fab(i,j-1,k,1);
      const Real ccc_fab_mns_z = ccc_fab(i,j-1,k,2);

      for (int n(0); n < ncomp; n++) {
        Real sy_ijkn(0);

        if( afrac_y > 0 ) {
          if( j <= dom_low.y and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_jlo_val)))
          {
            sy_ijkn = state(i,dom_low.y-1,k,state_comp+n);
          }
          else if( j >= dom_high.y+1 and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_jhi_val)))
          {
            sy_ijkn = state(i,dom_high.y+1,k,state_comp+n);
          }
          else
          {
            Real xf = fcy_fab_x; // local (x,z) of centroid of y-face we are extrapolating to
            Real zf = fcy_fab_y;

            Real delta_x = xf  - ccc_fab_x;
            Real delta_y = 0.5 + ccc_fab_y;
            Real delta_z = zf  - ccc_fab_z;

            Real state_pls = state(i,j  ,k,state_comp+n);
            Real state_mns = state(i,j-1,k,state_comp+n);

            Real cc_umax = amrex::max(state_pls, state_mns);
            Real cc_umin = amrex::min(state_pls, state_mns);

            Real vpls = state_pls - delta_y * y_slopes(i,j,k,slopes_comp+n)
                                  + delta_x * x_slopes(i,j,k,slopes_comp+n)
                                  + delta_z * z_slopes(i,j,k,slopes_comp+n);

            vpls = amrex::max( amrex::min(vpls, cc_umax), cc_umin );

            delta_x = xf  - ccc_fab_mns_x;
            delta_y = 0.5 - ccc_fab_mns_y;
            delta_z = zf  - ccc_fab_mns_z;

            Real vmns = state_mns + delta_y * y_slopes(i,j-1,k,slopes_comp+n)
                                  + delta_x * x_slopes(i,j-1,k,slopes_comp+n)
                                  + delta_z * z_slopes(i,j-1,k,slopes_comp+n);

            vmns = amrex::max( amrex::min(vmns, cc_umax), cc_umin );

            sy_ijkn = upwind(vmns, vpls, v_val);
          }
        }
        else {
          sy_ijkn = my_huge;
        }

        sy(i,j,k,n) = sy_ijkn;
        fy(i,j,k,n) = v_val * sy_ijkn;
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

      const Real afrac_z = areafrac_z(i,j,k);

      const Real fcz_fab_x = fcz_fab(i,j,k,0);
      const Real fcz_fab_y = fcz_fab(i,j,k,1);

      const Real ccc_fab_x = ccc_fab(i,j,k,0);
      const Real ccc_fab_y = ccc_fab(i,j,k,1);
      const Real ccc_fab_z = ccc_fab(i,j,k,2);

      const Real ccc_fab_mns_x = ccc_fab(i,j,k-1,0);
      const Real ccc_fab_mns_y = ccc_fab(i,j,k-1,1);
      const Real ccc_fab_mns_z = ccc_fab(i,j,k-1,2);

      for(int n(0); n < ncomp; n++) {
        Real sz_ijkn(0);

        if( afrac_z > 0 ) {
          if( k <= dom_low.z and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_klo_val)))
          {
            sz_ijkn = state(i, j,dom_low.z-1,state_comp+n);
          }
          else if( k >= dom_high.z+1 and
          aux::any_of(&bct_data[0], &bct_data[bct_size], aux::is_equal<int>(bct_khi_val)))
          {
            sz_ijkn = state(i, j,dom_high.z+1,state_comp+n);
          }
          else
          {
            Real xf = fcz_fab_x; // local (x,y) of centroid of z-face we are extrapolating to
            Real yf = fcz_fab_y;

            Real delta_x = xf - ccc_fab_x;
            Real delta_y = yf - ccc_fab_y;
            Real delta_z = .5 + ccc_fab_z;

            Real state_pls = state(i,j,k  ,state_comp+n);
            Real state_mns = state(i,j,k-1,state_comp+n);

            Real cc_umax = amrex::max(state_pls, state_mns);
            Real cc_umin = amrex::min(state_pls, state_mns);

            Real wpls = state_pls - delta_z * z_slopes(i,j,k,slopes_comp+n)
                                  + delta_x * x_slopes(i,j,k,slopes_comp+n)
                                  + delta_y * y_slopes(i,j,k,slopes_comp+n);

            wpls = amrex::max( amrex::min(wpls, cc_umax), cc_umin );

            delta_x = xf - ccc_fab_mns_x;
            delta_y = yf - ccc_fab_mns_y;
            delta_z = .5 - ccc_fab_mns_z;

            Real wmns = state_mns + delta_z * z_slopes(i,j,k-1,slopes_comp+n)
                                  + delta_x * x_slopes(i,j,k-1,slopes_comp+n)
                                  + delta_y * y_slopes(i,j,k-1,slopes_comp+n);

            wmns = amrex::max( amrex::min(wmns, cc_umax), cc_umin );

            sz_ijkn = upwind(wmns, wpls, w_val);
          }
        }
        else {
          sz_ijkn = my_huge;
        }

        sz(i,j,k,n) = sz_ijkn;
        fz(i,j,k,n) = w_val * sz_ijkn;
      }
    }
  });
}

}//end mol namespace

#endif
