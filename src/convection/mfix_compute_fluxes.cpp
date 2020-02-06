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
                           Vector< std::unique_ptr<MultiFab> >& a_fx,
                           Vector< std::unique_ptr<MultiFab> >& a_fy,
                           Vector< std::unique_ptr<MultiFab> >& a_fz,
                           Vector< std::unique_ptr<MultiFab> >& state_in,
                           const int state_comp, const int ncomp,
                           Vector< std::unique_ptr<MultiFab> >& xslopes_in,
                           Vector< std::unique_ptr<MultiFab> >& yslopes_in,
                           Vector< std::unique_ptr<MultiFab> >& zslopes_in,
                           const int slopes_comp,
                           Vector< std::unique_ptr<MultiFab> >& ep_u_mac,
                           Vector< std::unique_ptr<MultiFab> >& ep_v_mac,
                           Vector< std::unique_ptr<MultiFab> >& ep_w_mac)
{
        // Get EB geometric info
        Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
        Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;
        const amrex::MultiFab*                    volfrac;
        const amrex::MultiCutFab*                 bndrycent;


        areafrac  =   ebfactory[lev] -> getAreaFrac();
        facecent  =   ebfactory[lev] -> getFaceCent();
        volfrac   = &(ebfactory[lev] -> getVolFrac());
        bndrycent = &(ebfactory[lev] -> getBndryCent());

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
                // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
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
    Real state_w(0)  ;
    //Real state_e(0); DECLARED_BUT_NEVER_REFERENCED
    Real state_mns(0); Real state_pls(0);

    //
    // West face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if ((i == dom_low.x) and
     ugradu_aux::is_equal_to_any(bct_ilo(dom_low.x-1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_w = state(i-1,j,k,state_comp+n);

    } else if ((i == dom_high.x+1) and
     ugradu_aux::is_equal_to_any(bct_ihi(dom_high.x+1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_w = state(i,j,k,state_comp+n);
    } else {
      state_pls = state(i  ,j,k,state_comp+n) - .5*x_slopes(i  ,j,k,slopes_comp+n);
      state_mns = state(i-1,j,k,state_comp+n) + .5*x_slopes(i-1,j,k,slopes_comp+n);
      state_w = upwind( state_mns, state_pls, u(i,j,k) );
    }
    fx(i,j,k,n) = u(i,j,k) * state_w;
  });

  amrex::ParallelFor(vbx,ncomp,
    [slopes_comp,state_comp,dom_low,dom_high,bct_jlo,bct_jhi,bc_types,state,y_slopes,v,fy]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real state_s(0)  ;
    //Real state_n(0); DECLARED_BUT_NEVER_REFERENCED
    Real state_mns(0); Real state_pls(0);

    //
    // South face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((j == dom_low.y) and
     ugradu_aux::is_equal_to_any(bct_jlo(i,dom_low.y-1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_s = state(i,j-1,k,state_comp+n);
    } else if ((j == dom_high.y+1) and
     ugradu_aux::is_equal_to_any(bct_jhi(i,dom_high.y+1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_s = state(i,j,k,state_comp+n);
    } else {
      state_pls = state(i,j  ,k,state_comp+n) - .5*y_slopes(i,j  ,k,slopes_comp+n);
      state_mns = state(i,j-1,k,state_comp+n) + .5*y_slopes(i,j-1,k,slopes_comp+n);

      state_s = upwind( state_mns, state_pls, v(i,j,k) );
    }
    fy(i,j,k,n) = v(i,j,k) * state_s;
  });

  amrex::ParallelFor(wbx, ncomp,
    [slopes_comp,state_comp,dom_low,dom_high,bct_klo,bct_khi,bc_types,state,z_slopes,w,fz]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real state_b(0)  ;
    Real state_mns(0); Real state_pls(0);
    
    //
    // Bottom face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((k == dom_low.z) and
     ugradu_aux::is_equal_to_any(bct_klo(i,j,dom_low.z-1,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_b = state(i,j,k-1,state_comp+n);
    } else if ((k == dom_high.z+1) and
     ugradu_aux::is_equal_to_any(bct_khi(i,j,dom_high.z+1,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_b = state(i,j,k,state_comp+n);
    } else {
      state_pls = state(i,j,k  ,state_comp+n) - .5*z_slopes(i,j,k  ,slopes_comp+n);
      state_mns = state(i,j,k-1,state_comp+n) + .5*z_slopes(i,j,k-1,slopes_comp+n);

      state_b = upwind(state_mns, state_pls, w(i,j,k));
    }
    fz(i,j,k,n) = w(i,j,k) * state_b;
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

  const auto& ccm_fab = cc_mask.const_array();

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
    Real upls(0); Real umns(0);

    if( areafrac_x(i,j,k) > 0 ) {
      if( i <= dom_low.x and
           ugradu_aux::is_equal_to_any(bct_ilo(dom_low.x-1,j,k,0),
                                       bc_types.data(), bc_types.size()))
        {
         sx(i,j,k,n) = state(dom_low.x-1,j,k,state_comp+n);
      }
      else if( i >= dom_high.x+1 and
           ugradu_aux::is_equal_to_any(bct_ihi(dom_high.x+1,j,k,0),
                                       bc_types.data(), bc_types.size()))
      {
         sx(i,j,k,n) = state(dom_high.x+1,j,k,state_comp+n);
      }
      else 
      {
         Real yf = fcx_fab(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
         Real zf = fcx_fab(i,j,k,1);

         Real xc = ccc_fab(i,j,k,0); // centroid of cell (i,j,k)
         Real yc = ccc_fab(i,j,k,1);
         Real zc = ccc_fab(i,j,k,2);

         Real delta_x = 0.5 + xc;
         Real delta_y = yf  - yc;
         Real delta_z = zf  - zc;

         Real cc_umax = std::max(state(i,j,k,state_comp+n), state(i-1,j,k,state_comp+n));
         Real cc_umin = std::min(state(i,j,k,state_comp+n), state(i-1,j,k,state_comp+n));

         upls = state(i  ,j,k,state_comp+n) - delta_x * x_slopes(i,j,k,slopes_comp+n) 
                                            + delta_y * y_slopes(i,j,k,slopes_comp+n) 
                                            + delta_z * z_slopes(i,j,k,slopes_comp+n);

         upls = std::max( std::min(upls, cc_umax), cc_umin );

              xc = ccc_fab(i-1,j,k,0); // centroid of cell (i,j,k)
              yc = ccc_fab(i-1,j,k,1);
              zc = ccc_fab(i-1,j,k,2);

              delta_x = 0.5 - xc;
              delta_y = yf  - yc;
              delta_z = zf  - zc;

         umns = state(i-1,j,k,state_comp+n) + delta_x * x_slopes(i-1,j,k,slopes_comp+n) 
                                            + delta_y * y_slopes(i-1,j,k,slopes_comp+n) 
                                            + delta_z * z_slopes(i-1,j,k,slopes_comp+n);

         umns = std::max( std::min(umns, cc_umax), cc_umin );

         sx(i,j,k,n) = upwind( umns, upls, u(i,j,k) );
      }
    } else {
        sx(i,j,k,n) = my_huge;
    }
    fx(i,j,k,n) = u(i,j,k) * sx(i,j,k,n);
  });

  //
  // ===================== Y =====================
  //
  amrex::ParallelFor(vbx, ncomp,
    [my_huge,slopes_comp,state_comp,dom_low,dom_high,bct_jlo,bct_jhi,bc_types,areafrac_y,fcy_fab,ccc_fab,
     x_slopes,y_slopes,z_slopes,state,v,sy,fy]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real vpls(0); Real vmns(0);

    if( areafrac_y(i,j,k) > 0 ) {
      if( j <= dom_low.y and
       ugradu_aux::is_equal_to_any(bct_jlo(i,dom_low.y-1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sy(i,j,k,n) = state(i,dom_low.y-1,k,state_comp+n);
      }
      else if( j >= dom_high.y+1 and
       ugradu_aux::is_equal_to_any(bct_jhi(i,dom_high.y+1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sy(i,j,k,n) = state(i,dom_high.y+1,k,state_comp+n);
      }
      else 
      {
         Real xf = fcy_fab(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
         Real zf = fcy_fab(i,j,k,1);

         Real xc = ccc_fab(i,j,k,0); // centroid of cell (i,j,k)
         Real yc = ccc_fab(i,j,k,1);
         Real zc = ccc_fab(i,j,k,2);

         Real delta_x = xf  - xc;
         Real delta_y = 0.5 + yc;
         Real delta_z = zf  - zc;

         Real cc_umax = std::max(state(i,j,k,state_comp+n), state(i,j-1,k,state_comp+n));
         Real cc_umin = std::min(state(i,j,k,state_comp+n), state(i,j-1,k,state_comp+n));

         vpls = state(i,j  ,k,state_comp+n) - delta_y * y_slopes(i,j,k,slopes_comp+n) 
                                            + delta_x * x_slopes(i,j,k,slopes_comp+n) 
                                            + delta_z * z_slopes(i,j,k,slopes_comp+n);
         vpls = std::max( std::min(vpls, cc_umax), cc_umin );

              xc = ccc_fab(i,j-1,k,0); // centroid of cell (i,j-1,k)
              yc = ccc_fab(i,j-1,k,1);
              zc = ccc_fab(i,j-1,k,2);

              delta_x = xf  - xc;
              delta_y = 0.5 - yc;
              delta_z = zf  - zc;

         vmns = state(i,j-1,k,state_comp+n) + delta_y * y_slopes(i,j-1,k,slopes_comp+n) 
                                            + delta_x * x_slopes(i,j-1,k,slopes_comp+n) 
                                            + delta_z * z_slopes(i,j-1,k,slopes_comp+n);

         vmns = std::max( std::min(vmns, cc_umax), cc_umin );

         sy(i,j,k,n) = upwind( vmns, vpls, v(i,j,k) );
      }
    }
    else {
        sy(i,j,k,n) = my_huge;
    }
    fy(i,j,k,n) = v(i,j,k) * sy(i,j,k,n);
  });

  //
  // ===================== Z =====================
  //
  amrex::ParallelFor(wbx,ncomp,
    [my_huge,slopes_comp,state_comp,dom_low,dom_high,bct_klo,bct_khi,bc_types,areafrac_z,fcz_fab,ccc_fab,
     x_slopes,y_slopes,z_slopes,state,w,sz,fz]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real wpls(0); Real wmns(0);

    if( areafrac_z(i,j,k) > 0 ) {
      if( k <= dom_low.z and
       ugradu_aux::is_equal_to_any(bct_klo(i,j,dom_low.z-1,0),
                                   bc_types.data(), bc_types.size()))
      {
        sz(i,j,k,n) = state(i,j,dom_low.z-1,state_comp+n);
      }
      else if( k >= dom_high.z+1 and
       ugradu_aux::is_equal_to_any(bct_khi(i,j,dom_high.z+1,0),
                                   bc_types.data(), bc_types.size()))
      {
        sz(i,j,k,n) = state(i,j,dom_high.z+1,state_comp+n);
      }
      else 
      {
         Real xf = fcz_fab(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
         Real yf = fcz_fab(i,j,k,1);

         Real xc = ccc_fab(i,j,k,0); // centroid of cell (i,j,k)
         Real yc = ccc_fab(i,j,k,1);
         Real zc = ccc_fab(i,j,k,2);

         Real delta_x = xf  - xc;
         Real delta_y = yf  - yc;
         Real delta_z = 0.5 + zc;

         Real cc_umax = std::max(state(i,j,k,state_comp+n), state(i,j,k-1,state_comp+n));
         Real cc_umin = std::min(state(i,j,k,state_comp+n), state(i,j,k-1,state_comp+n));

         wpls = state(i,j,k  ,state_comp+n) - delta_z * z_slopes(i,j,k,slopes_comp+n) 
                                            + delta_x * x_slopes(i,j,k,slopes_comp+n) 
                                            + delta_y * y_slopes(i,j,k,slopes_comp+n);

         wpls = std::max( std::min(wpls, cc_umax), cc_umin );

              xc = ccc_fab(i,j,k-1,0); // centroid of cell (i,j,k-1)
              yc = ccc_fab(i,j,k-1,1);
              zc = ccc_fab(i,j,k-1,2);

              delta_x = xf  - xc;
              delta_y = yf  - yc;
              delta_z = 0.5 - zc;

         wmns = state(i,j,k-1,state_comp+n) + delta_z * z_slopes(i,j,k-1,slopes_comp+n) 
                                            + delta_x * x_slopes(i,j,k-1,slopes_comp+n) 
                                            + delta_y * y_slopes(i,j,k-1,slopes_comp+n);
         wmns = std::max( std::min(wmns, cc_umax), cc_umin );

         sz(i,j,k,n) = upwind( wmns, wpls, w(i,j,k) );
      }
    }
    else {
        sz(i,j,k,n) = my_huge;
    }
    fz(i,j,k,n) = w(i,j,k) * sz(i,j,k,n);
  });
}
