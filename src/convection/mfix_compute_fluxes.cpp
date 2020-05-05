#include <mfix.H>
#include <param_mod_F.H>

namespace ugradu_aux {

//
// Compute upwind non-normal velocity
//
AMREX_GPU_HOST_DEVICE
Real
upwind (const Real velocity_minus, const Real velocity_plus, const Real u_edge)
{
  // Small value to protect against tiny velocities used in upwinding
  const Real small_velocity(1.e-10);

  if(std::abs(u_edge) < small_velocity)
    return .5*(velocity_minus+velocity_plus);

  return u_edge > 0 ? velocity_minus : velocity_plus;
}

AMREX_GPU_HOST_DEVICE
bool
is_equal_to_any (const int bc, const int* bc_types, const int size)
{
  for(int i(0); i < size; ++i)
    if(bc == bc_types[i])
      return true;

  return false;
}

} // end namespace ugradu_aux

using namespace ugradu_aux;

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

  amrex::ParallelFor(ubx,ncomp,
    [slopes_comp,state_comp,dom_low,dom_high,bct_ilo,bct_ihi,bc_types,state,x_slopes,u,fx]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    //
    // West face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if ((i == dom_low.x) and
     ugradu_aux::is_equal_to_any(bct_ilo(dom_low.x-1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      fx(i,j,k,n) = u(i,j,k) * state(i-1,j,k,state_comp+n);
    }
    else if ((i == dom_high.x+1) and
     ugradu_aux::is_equal_to_any(bct_ihi(dom_high.x+1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      fx(i,j,k,n) = u(i,j,k) * state(i,j,k,state_comp+n);
    }
    else {
      const Real state_pls = state(i  ,j,k,state_comp+n) - .5*x_slopes(i  ,j,k,slopes_comp+n);
      const Real state_mns = state(i-1,j,k,state_comp+n) + .5*x_slopes(i-1,j,k,slopes_comp+n);
      const Real u_ijk = u(i,j,k);

      fx(i,j,k,n) = u_ijk * upwind(state_mns, state_pls, u_ijk);
    }
  });

  amrex::ParallelFor(vbx,ncomp,
    [slopes_comp,state_comp,dom_low,dom_high,bct_jlo,bct_jhi,bc_types,state,y_slopes,v,fy]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    //
    // South face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((j == dom_low.y) and
     ugradu_aux::is_equal_to_any(bct_jlo(i,dom_low.y-1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      fy(i,j,k,n) = v(i,j,k) * state(i,j-1,k,state_comp+n);
    }
    else if ((j == dom_high.y+1) and
     ugradu_aux::is_equal_to_any(bct_jhi(i,dom_high.y+1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      fy(i,j,k,n) = v(i,j,k) * state(i,j,k,state_comp+n);
    }
    else {
      const Real state_pls = state(i,j  ,k,state_comp+n) - .5*y_slopes(i,j  ,k,slopes_comp+n);
      const Real state_mns = state(i,j-1,k,state_comp+n) + .5*y_slopes(i,j-1,k,slopes_comp+n);
      const Real v_ijk = v(i,j,k);

      fy(i,j,k,n) = v_ijk * upwind(state_mns, state_pls, v_ijk);
    }
  });

  amrex::ParallelFor(wbx, ncomp,
    [slopes_comp,state_comp,dom_low,dom_high,bct_klo,bct_khi,bc_types,state,z_slopes,w,fz]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    //
    // Bottom face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((k == dom_low.z) and
     ugradu_aux::is_equal_to_any(bct_klo(i,j,dom_low.z-1,0),
                                 bc_types.data(), bc_types.size()))
    {
      fz(i,j,k,n) = w(i,j,k) * state(i,j,k-1,state_comp+n);
    }
    else if ((k == dom_high.z+1) and
     ugradu_aux::is_equal_to_any(bct_khi(i,j,dom_high.z+1,0),
                                 bc_types.data(), bc_types.size()))
    {
      fz(i,j,k,n) = w(i,j,k) * state(i,j,k,state_comp+n);
    }
    else {
      const Real state_pls = state(i,j,k  ,state_comp+n) - .5*z_slopes(i,j,k  ,slopes_comp+n);
      const Real state_mns = state(i,j,k-1,state_comp+n) + .5*z_slopes(i,j,k-1,slopes_comp+n);
      const Real w_ijk = w(i,j,k);

      fz(i,j,k,n) = w_ijk * upwind(state_mns, state_pls, w_ijk);
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

  const Real my_huge = get_my_huge();

  //
  // First compute the convective fluxes at the face center
  // Do this on ALL faces on the tile, i.e. INCLUDE as many ghost faces as
  // possible
  //

  //
  // ===================== X =====================
  //
  amrex::ParallelFor(ubx,ncomp,
    [my_huge,slopes_comp,state_comp,dom_low,dom_high,bct_ilo,bct_ihi,bc_types,areafrac_x,fcx_fab,ccc_fab,
     x_slopes,y_slopes,z_slopes,state,u,sx,fx]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    if( areafrac_x(i,j,k) > 0 )
    {
      if(i <= dom_low.x and
        ugradu_aux::is_equal_to_any(bct_ilo(dom_low.x-1,j,k,0),
                                    bc_types.data(), bc_types.size()))
      {
        Real sx_ijkn = state(dom_low.x-1,j,k,state_comp+n);

        sx(i,j,k,n) = sx_ijkn;
        fx(i,j,k,n) = u(i,j,k) * sx_ijkn;
      }
      else if(i >= dom_high.x+1 and
        ugradu_aux::is_equal_to_any(bct_ihi(dom_high.x+1,j,k,0),
                                    bc_types.data(), bc_types.size()))
      {
        Real sx_ijkn = state(dom_high.x+1,j,k,state_comp+n);

        sx(i,j,k,n) = sx_ijkn;
        fx(i,j,k,n) = u(i,j,k) * sx_ijkn;
      }
      else 
      {
        Real yf = fcx_fab(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
        Real zf = fcx_fab(i,j,k,1);

        Real delta_x = .5 + ccc_fab(i,j,k,0);
        Real delta_y = yf - ccc_fab(i,j,k,1);
        Real delta_z = zf - ccc_fab(i,j,k,2);

        Real state_pls = state(i,j,k,state_comp+n);
        Real state_mns = state(i-1,j,k,state_comp+n);

        Real cc_umax = std::max(state_pls, state_mns);
        Real cc_umin = std::min(state_pls, state_mns);

        Real upls = state_pls - delta_x * x_slopes(i,j,k,slopes_comp+n) 
                              + delta_y * y_slopes(i,j,k,slopes_comp+n) 
                              + delta_z * z_slopes(i,j,k,slopes_comp+n);

        upls = std::max( std::min(upls, cc_umax), cc_umin );

        delta_x = .5 - ccc_fab(i-1,j,k,0);
        delta_y = yf - ccc_fab(i-1,j,k,1);
        delta_z = zf - ccc_fab(i-1,j,k,2);

        Real umns = state_mns + delta_x * x_slopes(i-1,j,k,slopes_comp+n) 
                              + delta_y * y_slopes(i-1,j,k,slopes_comp+n) 
                              + delta_z * z_slopes(i-1,j,k,slopes_comp+n);

        umns = std::max( std::min(umns, cc_umax), cc_umin );

        const Real u_ijk = u(i,j,k);
        const Real sx_ijkn = upwind(umns, upls, u_ijk);

        sx(i,j,k,n) = sx_ijkn;
        fx(i,j,k,n) = u_ijk * sx_ijkn;
      }
    }
    else {
      const Real sx_ijkn = my_huge;

      sx(i,j,k,n) = sx_ijkn;
      fx(i,j,k,n) = u(i,j,k) * sx_ijkn;
    }
  });

  //
  // ===================== Y =====================
  //
  amrex::ParallelFor(vbx, ncomp,
    [my_huge,slopes_comp,state_comp,dom_low,dom_high,bct_jlo,bct_jhi,bc_types,areafrac_y,fcy_fab,ccc_fab,
     x_slopes,y_slopes,z_slopes,state,v,sy,fy]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    if( areafrac_y(i,j,k) > 0 ) {
      if( j <= dom_low.y and
       ugradu_aux::is_equal_to_any(bct_jlo(i,dom_low.y-1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        Real sy_ijkn = state(i,dom_low.y-1,k,state_comp+n);

        sy(i,j,k,n) = sy_ijkn;
        fy(i,j,k,n) = v(i,j,k) * sy_ijkn;
      }
      else if( j >= dom_high.y+1 and
       ugradu_aux::is_equal_to_any(bct_jhi(i,dom_high.y+1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        Real sy_ijkn = state(i,dom_high.y+1,k,state_comp+n);

        sy(i,j,k,n) = sy_ijkn;
        fy(i,j,k,n) = v(i,j,k) * sy_ijkn;
      }
      else 
      {
        Real xf = fcy_fab(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
        Real zf = fcy_fab(i,j,k,1);

        Real delta_x = xf  - ccc_fab(i,j,k,0);
        Real delta_y = 0.5 + ccc_fab(i,j,k,1);
        Real delta_z = zf  - ccc_fab(i,j,k,2);

        Real state_pls = state(i,j  ,k,state_comp+n);
        Real state_mns = state(i,j-1,k,state_comp+n);

        Real cc_umax = std::max(state_pls, state_mns);
        Real cc_umin = std::min(state_pls, state_mns);

        Real vpls = state_pls - delta_y * y_slopes(i,j,k,slopes_comp+n) 
                              + delta_x * x_slopes(i,j,k,slopes_comp+n) 
                              + delta_z * z_slopes(i,j,k,slopes_comp+n);

        vpls = std::max( std::min(vpls, cc_umax), cc_umin );

        delta_x = xf  - ccc_fab(i,j-1,k,0);
        delta_y = 0.5 - ccc_fab(i,j-1,k,1);
        delta_z = zf  - ccc_fab(i,j-1,k,2);

        Real vmns = state_mns + delta_y * y_slopes(i,j-1,k,slopes_comp+n) 
                              + delta_x * x_slopes(i,j-1,k,slopes_comp+n) 
                              + delta_z * z_slopes(i,j-1,k,slopes_comp+n);

        vmns = std::max( std::min(vmns, cc_umax), cc_umin );

        const Real v_ijk = v(i,j,k);
        const Real sy_ijkn = upwind(vmns, vpls, v_ijk);

        sy(i,j,k,n) = sy_ijkn;
        fy(i,j,k,n) = v_ijk * sy_ijkn;
      }
    }
    else {
      const Real sy_ijkn = my_huge; 

      sy(i,j,k,n) = sy_ijkn;
      fy(i,j,k,n) = v(i,j,k) * sy_ijkn;
    }                      
  });                      
                           
  //                       
  // =====================  Z =====================
  //                       
  amrex::ParallelFor(wbx, ncomp,
    [my_huge,slopes_comp,state_comp,dom_low,dom_high,bct_klo,bct_khi,bc_types,areafrac_z,fcz_fab,ccc_fab,
     x_slopes,y_slopes,z_slopes,state,w,sz,fz]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {                        
    if( areafrac_z(i,j,k) > 0 ) {
      if( k <= dom_low.z and
       ugradu_aux::is_equal_to_any(bct_klo(i,j,dom_low.z-1,0),
                                   bc_types.data(), bc_types.size()))
      {                    
        const Real sz_ijkn = state(i, j,dom_low.z-1,state_comp+n);

        sz(i,j,k,n) = sz_ijkn;
        fz(i,j,k,n) = w(i,j,k) * sz_ijkn;
      }                    
      else if( k >= dom_high.z+1 and
       ugradu_aux::is_equal_to_any(bct_khi(i,j,dom_high.z+1,0),
                                   bc_types.data(), bc_types.size()))
      {                    
        const Real sz_ijkn = state(i, j,dom_high.z+1,state_comp+n);

        sz(i,j,k,n) = sz_ijkn;
        fz(i,j,k,n) = w(i,j,k) * sz_ijkn;
      }                    
      else                 
      {                    
        Real xf = fcz_fab(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
        Real yf = fcz_fab(i,j,k,1);

        Real delta_x = xf - ccc_fab(i,j,k,0);
        Real delta_y = yf - ccc_fab(i,j,k,1);
        Real delta_z = .5 + ccc_fab(i,j,k,2);

        Real state_pls = state(i,j,k  ,state_comp+n);
        Real state_mns = state(i,j,k-1,state_comp+n);

        Real cc_umax = std::max(state_pls, state_mns);
        Real cc_umin = std::min(state_pls, state_mns);

        Real wpls = state_pls - delta_z * z_slopes(i,j,k,slopes_comp+n) 
                              + delta_x * x_slopes(i,j,k,slopes_comp+n) 
                              + delta_y * y_slopes(i,j,k,slopes_comp+n);

        wpls = std::max( std::min(wpls, cc_umax), cc_umin );

        delta_x = xf - ccc_fab(i,j,k-1,0);
        delta_y = yf - ccc_fab(i,j,k-1,1);
        delta_z = .5 - ccc_fab(i,j,k-1,2);

        Real wmns = state_mns + delta_z * z_slopes(i,j,k-1,slopes_comp+n) 
                              + delta_x * x_slopes(i,j,k-1,slopes_comp+n) 
                              + delta_y * y_slopes(i,j,k-1,slopes_comp+n);

        wmns = std::max( std::min(wmns, cc_umax), cc_umin );

        const Real w_ijk = w(i,j,k);
        const Real sz_ijkn = upwind(wmns, wpls, w_ijk);

        sz(i,j,k,n) = sz_ijkn;
        fz(i,j,k,n) = w_ijk * sz_ijkn;
      }
    }
    else {
      const Real sz_ijkn = my_huge;

      sz(i,j,k,n) = sz_ijkn;
      fz(i,j,k,n) = w(i,j,k) * sz_ijkn;
    }
  });
}
