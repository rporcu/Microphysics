#include <mfix.H>
#include <param_mod_F.H>

namespace aux {

struct is_equal {
  AMREX_GPU_HOST_DEVICE
  is_equal () {}

  template <class T>
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  bool operator() (const T& x, const T& y) const {return x == y;}
};

template <class Operator> bool
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
any (const int bc, const int* bc_types, const int size, Operator op)
{
  for(int i(0); i < size; ++i)
    if(op(bc, bc_types[i]))
      return true;

  return false;
}

//
// Compute upwind non-normal velocity
//
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real
upwind (const Real velocity_minus, const Real velocity_plus, const Real u_edge)
{
  // Small value to protect against tiny velocities used in upwinding
  const Real small_velocity(1.e-10);

  if(std::abs(u_edge) < small_velocity)
    return .5*(velocity_minus+velocity_plus);

  return u_edge > 0 ? velocity_minus : velocity_plus;
}

} // end aux namespace

using namespace aux;

//
// Compute the three components of the convection term
//
void
mfix::mfix_compute_fluxes (int lev,
                           Vector< MultiFab* >& a_fx,
                           Vector< MultiFab* >& a_fy,
                           Vector< MultiFab* >& a_fz,
                           Vector< MultiFab* > const& state_in,
                           const int state_comp, const int ncomp,
                           Vector< MultiFab* > const& xslopes_in,
                           Vector< MultiFab* > const& yslopes_in,
                           Vector< MultiFab* > const& zslopes_in,
                           const int slopes_comp,
                           Vector< MultiFab* > const& ep_u_mac,
                           Vector< MultiFab* > const& ep_v_mac,
                           Vector< MultiFab* > const& ep_w_mac)
{
  // Get EB geometric info
  Array< const MultiCutFab*,3> areafrac;
  Array< const MultiCutFab*,3> facecent;
  const amrex::MultiFab*    volfrac;
  const amrex::MultiCutFab* bndrycent;

  areafrac  =   ebfactory[lev]->getAreaFrac();
  facecent  =   ebfactory[lev]->getFaceCent();
  volfrac   = &(ebfactory[lev]->getVolFrac());
  bndrycent = &(ebfactory[lev]->getBndryCent());

  const auto& cellcent = ebfactory[lev]->getCentroid();

  // Create cc_mask
  iMultiFab cc_mask(grids[lev], dmap[lev], 1, 1);

  const int covered_value = 1;
  const int notcovered_value = 0;
  const int physical_boundaries_value = 0;
  const int interior_value = 1;

  cc_mask.BuildMask(geom[lev].Domain(), geom[lev].periodicity(),
                    covered_value, notcovered_value,
                    physical_boundaries_value, interior_value);

  // We do this here to avoid any confusion about the FAB setVal.
  a_fx[lev]->setVal(covered_val);
  a_fy[lev]->setVal(covered_val);
  a_fz[lev]->setVal(covered_val);

  for (MFIter mfi(*state_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    // Tilebox
    Box bx = mfi.tilebox ();

    // this is to check efficiently if this tile contains any eb stuff
    const EBFArrayBox& state_fab = static_cast<EBFArrayBox const&>((*state_in[lev])[mfi]);
    const EBCellFlagFab&  flags = state_fab.getEBCellFlagFab();

    if (flags.getType(amrex::grow(bx,0)) != FabType::covered )
    {
      // No cut cells in tile + nghost-cell width halo -> use non-eb routine
      if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
      {
        mfix_compute_fluxes_on_box(
              lev, bx, (*a_fx[lev])[mfi], (*a_fy[lev])[mfi], (*a_fz[lev])[mfi], 
              (*state_in[lev])[mfi], state_comp, ncomp,
              (*xslopes_in[lev])[mfi], (*yslopes_in[lev])[mfi], (*zslopes_in[lev])[mfi], slopes_comp,
              (*ep_u_mac[lev])[mfi], (*ep_v_mac[lev])[mfi], (*ep_w_mac[lev])[mfi]);
      }
      else
      {
        mfix_compute_eb_fluxes_on_box(
              lev, bx, (*a_fx[lev])[mfi], (*a_fy[lev])[mfi], (*a_fz[lev])[mfi], 
              (*state_in[lev])[mfi], state_comp, ncomp,
              (*xslopes_in[lev])[mfi], (*yslopes_in[lev])[mfi], (*zslopes_in[lev])[mfi], slopes_comp,
              (*ep_u_mac[lev])[mfi], (*ep_v_mac[lev])[mfi], (*ep_w_mac[lev])[mfi],
              (*areafrac[0])[mfi], (*areafrac[1])[mfi], (*areafrac[2])[mfi], 
              (*facecent[0])[mfi], (*facecent[1])[mfi], (*facecent[2])[mfi], 
              cellcent[mfi], (*volfrac)[mfi], (*bndrycent)[mfi], cc_mask[mfi], flags);
      }
    }
  } // MFIter
}

void
mfix::mfix_compute_fluxes_on_box (const int lev, Box& bx,
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
                                  const FArrayBox& ep_w_mac)
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

  Array4<int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<int> const& bct_klo = bc_klo[lev]->array();
  Array4<int> const& bct_khi = bc_khi[lev]->array();

  const Box ubx = amrex::surroundingNodes(bx,0);
  const Box vbx = amrex::surroundingNodes(bx,1);
  const Box wbx = amrex::surroundingNodes(bx,2);

  // Vectorize the boundary conditions list in order to use it in lambda
  // functions
  const GpuArray<int, 3> bc_types =
    {bc_list.get_minf(), bc_list.get_pinf(), bc_list.get_pout()};

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

  ParallelFor(npoints, [=] AMREX_GPU_DEVICE (int idx) noexcept
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
          any(bct_ilo_val, bct_data, bct_size, aux::is_equal()))
        {
          fx_val = u_val * state_mns;
        }
        else if ((i == dom_high.x+1) and
            any(bct_ihi_val, bct_data, bct_size, aux::is_equal()))
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
          any(bct_jlo_val, bct_data, bct_size, aux::is_equal()))
        {
          fy_val = v_val * state_mns;
        }
        else if ((j == dom_high.y+1) and
            any(bct_jhi_val, bct_data, bct_size, aux::is_equal()))
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
          any(bct_klo_val, bct_data, bct_size, aux::is_equal()))
        {
          fz_val = w_val * state_mns;
        }
        else if ((k == dom_high.z+1) and
            any(bct_khi_val, bct_data, bct_size, aux::is_equal()))
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

//
// Compute the three components of the convection term when we have embedded
// boundaries
//
void
mfix::mfix_compute_eb_fluxes_on_box (const int lev, Box& bx,
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
                                     const IArrayBox& cc_mask, 
                                     const EBCellFlagFab& flags)
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

  Array4<int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<int> const& bct_klo = bc_klo[lev]->array();
  Array4<int> const& bct_khi = bc_khi[lev]->array();

  const Box ubx = amrex::surroundingNodes(bx,0);
  const Box vbx = amrex::surroundingNodes(bx,1);
  const Box wbx = amrex::surroundingNodes(bx,2);

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

  const GpuArray<int, 3> bc_types =
    {bc_list.get_minf(), bc_list.get_pinf(), bc_list.get_pout()};

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

  ParallelFor(npoints, [=] AMREX_GPU_DEVICE (int idx) noexcept
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
            any(bct_ilo_val, bct_data, bct_size, aux::is_equal()))
          {
            sx_ijkn = state(dom_low.x-1,j,k,state_comp+n);
          }
          else if(i >= dom_high.x+1 and
            any(bct_ihi_val, bct_data, bct_size, aux::is_equal()))
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
           any(bct_jlo_val, bct_data, bct_size, aux::is_equal()))
          {
            sy_ijkn = state(i,dom_low.y-1,k,state_comp+n);
          }
          else if( j >= dom_high.y+1 and
           any(bct_jhi_val, bct_data, bct_size, aux::is_equal()))
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
           any(bct_klo_val, bct_data, bct_size, aux::is_equal()))
          {                    
            sz_ijkn = state(i, j,dom_low.z-1,state_comp+n);
          }                    
          else if( k >= dom_high.z+1 and
           any(bct_khi_val, bct_data, bct_size, aux::is_equal()))
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
