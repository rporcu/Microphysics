#include <mfix.H>
#include <mfix_divop_conv.hpp>
#include <param_mod_F.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>

namespace ugradu_aux {

//
// Compute upwind non-normal velocity
//
AMREX_GPU_HOST_DEVICE
Real
upwind(const Real velocity_minus,
       const Real velocity_plus,
       const Real u_edge)
{
  // Small value to protect against tiny velocities used in upwinding
  const Real small_velocity(1.e-10);

  if(std::abs(u_edge) < small_velocity)
    return .5*(velocity_minus+velocity_plus);

  return u_edge > 0 ? velocity_minus : velocity_plus;
}

AMREX_GPU_HOST_DEVICE
bool
is_equal_to_any(const int bc,
                const int* bc_types,
                const int size)
{
  for(int i(0); i < size; ++i)
  {
    if(bc == bc_types[i])
      return true;
  }
  return false;
}

} // end namespace ugradu_aux

using namespace ugradu_aux;

//
// Compute the three components of the convection term
//
void 
mfix::mfix_compute_fluxes(int lev,
                          Vector< std::unique_ptr<MultiFab> >& conv_in, 
                          const int conv_comp, 
                          Vector< std::unique_ptr<MultiFab> >& state_in,
                          const int state_comp, const int ncomp,
                          Vector< std::unique_ptr<MultiFab> >& xslopes_in,
                          Vector< std::unique_ptr<MultiFab> >& yslopes_in,
                          Vector< std::unique_ptr<MultiFab> >& zslopes_in,
                          const int slopes_comp,
                          Vector< std::unique_ptr<MultiFab> >& u_mac,
                          Vector< std::unique_ptr<MultiFab> >& v_mac,
                          Vector< std::unique_ptr<MultiFab> >& w_mac,
                          const bool is_conservative)
{
        Box domain(geom[lev].Domain());

        // Get EB geometric info
        Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
        Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;
        const amrex::MultiFab*                    volfrac;
        const amrex::MultiCutFab*                 bndrycent;

        areafrac  =   ebfactory[lev] -> getAreaFrac();
        facecent  =   ebfactory[lev] -> getFaceCent();
        volfrac   = &(ebfactory[lev] -> getVolFrac());
        bndrycent = &(ebfactory[lev] -> getBndryCent());

        // Create cc_mask
        iMultiFab cc_mask(grids[lev], dmap[lev], 1, 1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       {
           std::vector< std::pair<int,Box> > isects;
           const std::vector<IntVect>& pshifts = geom[lev].periodicity().shiftIntVect();
           const BoxArray& ba = cc_mask.boxArray();
           for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
           {
               Array4<int> const& fab = cc_mask.array(mfi);
               const Box& bx = mfi.fabbox();
               for (const auto& iv : pshifts)
               {
                   ba.intersections(bx+iv, isects);
                   for (const auto& is : isects)
                   {
                       const Box& b = is.second-iv;
                       AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( b, i, j, k,
                       {
                           fab(i,j,k) = 1;
                       });
                   }
               }
           }
        }

        for (MFIter mfi(*state_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox ();

            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox&  conv_fab = static_cast<EBFArrayBox const&>((*conv_in[lev])[mfi]);
            const EBCellFlagFab&  flags = conv_fab.getEBCellFlagFab();

            if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
            {
                // If tile is completely covered by EB geometry, set slopes
                // value to some very large number so we know if
                // we accidentally use these covered slopes later in calculations
                (*conv_in[lev])[mfi].setVal( 1.2345e300, bx, 0, conv_in[lev]->nComp());
            }
            else
            {
                // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
                if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
                {
                    mfix_compute_ugradu(bx, conv_in, conv_comp, state_in, state_comp, ncomp,
                                        xslopes_u, yslopes_u, zslopes_u, slopes_comp,
                                        u_mac, v_mac, w_mac, &mfi, domain, lev, false);
                }
                else
                {
                    mfix_compute_ugradu_eb(bx, conv_in, conv_comp, state_in, state_comp, ncomp,
                                           xslopes_u, yslopes_u, zslopes_u, slopes_comp,
                                           u_mac, v_mac, w_mac, &mfi, areafrac, facecent,
                                           volfrac, bndrycent, &cc_mask, domain, flags, lev, false);
                }
            }
        } // MFIter
}

void
mfix::mfix_compute_ugradu( Box& bx,
                           Vector< std::unique_ptr<MultiFab> >& conv_in, const int conv_comp, 
                           Vector< std::unique_ptr<MultiFab> >& state_in,
                           const int state_comp, const int ncomp,
                           Vector< std::unique_ptr<MultiFab> >& xslopes_in,
                           Vector< std::unique_ptr<MultiFab> >& yslopes_in,
                           Vector< std::unique_ptr<MultiFab> >& zslopes_in,
                           const int slopes_comp,
                           Vector< std::unique_ptr<MultiFab> >& u_mac,
                           Vector< std::unique_ptr<MultiFab> >& v_mac,
                           Vector< std::unique_ptr<MultiFab> >& w_mac,
                           MFIter* mfi,
                           Box& domain,
                           const int lev,
                           const bool is_conservative)
{
  const Real* dx = geom[lev].CellSize();
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& ugradu = conv_in[lev]->array(*mfi); 
  
  Array4<Real> const& state     = state_in[lev]->array(*mfi);
  Array4<Real> const& epsilon_g = ep_g[lev]->array(*mfi);
  
  Array4<Real> const& u = u_mac[lev]->array(*mfi);
  Array4<Real> const& v = v_mac[lev]->array(*mfi);
  Array4<Real> const& w = w_mac[lev]->array(*mfi);

  Array4<Real> const& x_slopes = xslopes_in[lev]->array(*mfi);
  Array4<Real> const& y_slopes = yslopes_in[lev]->array(*mfi);
  Array4<Real> const& z_slopes = zslopes_in[lev]->array(*mfi);

  Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
  Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();
  Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
  Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();
  Array4<int> const& bc_klo_type = bc_klo[lev]->array();
  Array4<int> const& bc_khi_type = bc_khi[lev]->array();

  const Real i_dx(1/dx[0]), i_dy(1/dx[1]), i_dz(1/dx[2]);

  // Vectorize the boundary conditions list in order to use it in lambda
  // functions
  const GpuArray<int, 3> bc_types =
    {bc_list.get_minf(), bc_list.get_pinf(), bc_list.get_pout()};

  AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
  {
    Real state_w(0)  ; Real state_e(0);
    Real state_s(0)  ; Real state_n(0);
    Real state_b(0)  ; Real state_t(0);
    Real state_mns(0); Real state_pls(0);

    //
    // West face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((i == dom_low.x) and
     ugradu_aux::is_equal_to_any(bc_ilo_type(dom_low.x-1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_w = state(i-1,j,k,state_comp+n);
    } else {
      state_pls = state(i  ,j,k,state_comp+n) - .5*x_slopes(i  ,j,k,slopes_comp+n);
      state_mns = state(i-1,j,k,state_comp+n) + .5*x_slopes(i-1,j,k,slopes_comp+n);
      state_w = upwind( state_mns, state_pls, u(i,j,k) );
    }

    //
    // East face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((i == dom_high.x) and
     ugradu_aux::is_equal_to_any(bc_ihi_type(dom_high.x+1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_e = state(i+1,j,k,state_comp+n);
    } else {
      state_pls = state(i+1,j,k,state_comp+n) - .5*x_slopes(i+1,j,k,slopes_comp+n);
      state_mns = state(i  ,j,k,state_comp+n) + .5*x_slopes(i  ,j,k,slopes_comp+n);

      state_e = upwind( state_mns, state_pls, u(i+1,j,k) );
    }

    //
    // South face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((j == dom_low.y) and
     ugradu_aux::is_equal_to_any(bc_jlo_type(i,dom_low.y-1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_s = state(i,j-1,k,state_comp+n);
    } else {
      state_pls = state(i,j  ,k,state_comp+n) - .5*y_slopes(i,j  ,k,slopes_comp+n);
      state_mns = state(i,j-1,k,state_comp+n) + .5*y_slopes(i,j-1,k,slopes_comp+n);

      state_s = upwind( state_mns, state_pls, v(i,j,k) );
    }

    //
    // North face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((j == dom_high.y) and
     ugradu_aux::is_equal_to_any(bc_jhi_type(i,dom_high.y+1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_n = state(i,j+1,k,state_comp+n);
    } else {
      state_pls = state(i,j+1,k,state_comp+n) - .5*y_slopes(i,j+1,k,slopes_comp+n);
      state_mns = state(i,j  ,k,state_comp+n) + .5*y_slopes(i,j  ,k,slopes_comp+n);

      state_n = upwind( state_mns, state_pls, v(i,j+1,k) );
    }

    //
    // Bottom face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((k == dom_low.z) and
     ugradu_aux::is_equal_to_any(bc_klo_type(i,j,dom_low.z-1,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_b = state(i,j,k-1,state_comp+n);
    } else {
      state_pls = state(i,j,k  ,state_comp+n) - .5*z_slopes(i,j,k  ,slopes_comp+n);
      state_mns = state(i,j,k-1,state_comp+n) + .5*z_slopes(i,j,k-1,slopes_comp+n);

      state_b = upwind( state_mns, state_pls, w(i,j,k) );
    }

    //
    // Top face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((k == dom_high.z) and
     ugradu_aux::is_equal_to_any(bc_khi_type(i,j,dom_high.z+1,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_t = state(i,j,k+1,state_comp+n);
    } else {
      state_pls = state(i,j,k+1,state_comp+n) - .5*z_slopes(i,j,k+1,slopes_comp+n);
      state_mns = state(i,j,k  ,state_comp+n) + .5*z_slopes(i,j,k  ,slopes_comp+n);

      state_t = upwind( state_mns, state_pls, w(i,j,k+1) );
    }

    Real epu_hi_x = u(i+1,j,k);
    Real epu_lo_x = u(i  ,j,k);
    Real epv_hi_y = v(i,j+1,k);
    Real epv_lo_y = v(i,j  ,k);
    Real epw_hi_z = w(i,j,k+1);
    Real epw_lo_z = w(i,j,k  );

    ugradu(i,j,k,conv_comp+n) = (epu_hi_x*state_e - epu_lo_x*state_w) * i_dx +
                                (epv_hi_y*state_n - epv_lo_y*state_s) * i_dy +
                                (epw_hi_z*state_t - epw_lo_z*state_b) * i_dz;
    if (!is_conservative)
    {
       Real diveumac = (epu_hi_x - epu_lo_x) * i_dx +
                       (epv_hi_y - epv_lo_y) * i_dy + 
                       (epw_hi_z - epw_lo_z) * i_dz;
       ugradu(i,j,k,conv_comp+n) = ugradu(i,j,k,conv_comp+n) - state(i,j,k,state_comp+n)*diveumac;
    }

    //
    // Return the negative
    //
    const Real coefficient(-1/epsilon_g(i,j,k));
    ugradu(i,j,k,conv_comp+n) *= coefficient; 
  });

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif
}


//
// Compute the three components of the convection term when we have embedded
// boundaries
//
void
mfix::mfix_compute_ugradu_eb(Box& bx,
                             Vector< std::unique_ptr<MultiFab> >& conv_in, 
                             const int conv_comp, 
                             Vector< std::unique_ptr<MultiFab> >& state_in,
                             const int state_comp, const int ncomp,
                             Vector< std::unique_ptr<MultiFab> >& xslopes_in,
                             Vector< std::unique_ptr<MultiFab> >& yslopes_in,
                             Vector< std::unique_ptr<MultiFab> >& zslopes_in,
                             const int slopes_comp,
                             Vector< std::unique_ptr<MultiFab> >& u_mac,
                             Vector< std::unique_ptr<MultiFab> >& v_mac,
                             Vector< std::unique_ptr<MultiFab> >& w_mac,
                             MFIter* mfi,
                             Array<const MultiCutFab*,AMREX_SPACEDIM>& areafrac,
                             Array<const MultiCutFab*,AMREX_SPACEDIM>& facecent,
                             const MultiFab* volfrac,
                             const MultiCutFab* bndrycent,
                             const iMultiFab* cc_mask,
                             Box& domain,
                             const EBCellFlagFab& flags,
                             const int lev,
                             const bool is_conservative)
{
  AMREX_ASSERT_WITH_MESSAGE(nghost >= 4, "Compute divop_conv(): ng must be >= 4");

  const Real* dx = geom[lev].CellSize();
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& ugradu = conv_in[lev]->array(*mfi);

  Array4<Real> const& state     = state_in[lev]->array(*mfi);
  Array4<Real> const& epsilon_g =    ep_g[lev]->array(*mfi);

  Array4<const Real> const& areafrac_x = areafrac[0]->array(*mfi);
  Array4<const Real> const& areafrac_y = areafrac[1]->array(*mfi);
  Array4<const Real> const& areafrac_z = areafrac[2]->array(*mfi);

  Array4<Real> const& u = u_mac[lev]->array(*mfi);
  Array4<Real> const& v = v_mac[lev]->array(*mfi);
  Array4<Real> const& w = w_mac[lev]->array(*mfi);

  Array4<Real> const& x_slopes = xslopes_in[lev]->array(*mfi);
  Array4<Real> const& y_slopes = yslopes_in[lev]->array(*mfi);
  Array4<Real> const& z_slopes = zslopes_in[lev]->array(*mfi);

  Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
  Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();
  Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
  Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();
  Array4<int> const& bc_klo_type = bc_klo[lev]->array();
  Array4<int> const& bc_khi_type = bc_khi[lev]->array();

  // Number of Halo layers
  const int nh(3);

  const Box ubx = amrex::surroundingNodes(amrex::grow(bx,nh),0);
  const Box vbx = amrex::surroundingNodes(amrex::grow(bx,nh),1);
  const Box wbx = amrex::surroundingNodes(amrex::grow(bx,nh),2);

  const Box ubx_inner = amrex::surroundingNodes(amrex::grow(bx,nh-1),0);
  const Box vbx_inner = amrex::surroundingNodes(amrex::grow(bx,nh-1),1);
  const Box wbx_inner = amrex::surroundingNodes(amrex::grow(bx,nh-1),2);
  
  FArrayBox fxfab(ubx, ncomp);
  FArrayBox fyfab(vbx, ncomp);
  FArrayBox fzfab(wbx, ncomp);

  Array4<Real> const& fx = fxfab.array();
  Array4<Real> const& fy = fyfab.array();
  Array4<Real> const& fz = fzfab.array();
  
  FArrayBox s_on_x_face(ubx, ncomp);
  FArrayBox s_on_y_face(vbx, ncomp);
  FArrayBox s_on_z_face(wbx, ncomp);

  Array4<Real> const& sx = s_on_x_face.array();
  Array4<Real> const& sy = s_on_y_face.array();
  Array4<Real> const& sz = s_on_z_face.array();

  // Face centroids
  const auto& fcx_fab = facecent[0]->array(*mfi);
  const auto& fcy_fab = facecent[1]->array(*mfi);
  const auto& fcz_fab = facecent[2]->array(*mfi);

  const auto& ccm_fab = cc_mask->const_array(*mfi);

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
  AMREX_HOST_DEVICE_FOR_4D(ubx, ncomp, i, j, k, n,
  {
    Real upls(0); Real umns(0);

    if( areafrac_x(i,j,k) > 0 ) {
      if( i <= dom_low.x and
       ugradu_aux::is_equal_to_any(bc_ilo_type(dom_low.x-1,j,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sx(i,j,k,n) = state(dom_low.x-1,j,k,state_comp+n);
      }
      else if( i >= dom_high.x+1 and
       ugradu_aux::is_equal_to_any(bc_ihi_type(dom_high.x+1,j,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sx(i,j,k,n) = state(dom_high.x+1,j,k,state_comp+n);
      }
      else {
        upls = state(i  ,j,k,state_comp+n) - .5*x_slopes(i  ,j,k,slopes_comp+n);
        umns = state(i-1,j,k,state_comp+n) + .5*x_slopes(i-1,j,k,slopes_comp+n);

        sx(i,j,k,n) = upwind( umns, upls, u(i,j,k) );
      }
    }
    else {
        sx(i,j,k,n) = my_huge; 
    }
  });

  AMREX_HOST_DEVICE_FOR_4D(ubx, ncomp, i, j, k, n,
  {
    int jj = j + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,0)));
    int kk = k + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,1)));

    Real fracy = (ccm_fab(i-1,jj,k) || ccm_fab(i,jj,k)) ? std::abs(fcx_fab(i,j,k,0)) : 0.0;
    Real fracz = (ccm_fab(i-1,j,kk) || ccm_fab(i,j,kk)) ? std::abs(fcx_fab(i,j,k,1)) : 0.0;

    Real s_on_x_centroid = (1.0-fracy)*(1.0-fracz)*sx(i, j,k ,n)+
                                fracy *(1.0-fracz)*sx(i,jj,k ,n)+
                                fracz *(1.0-fracy)*sx(i, j,kk,n)+
                                fracy *     fracz *sx(i,jj,kk,n);

    fx(i,j,k,n) = u(i,j,k) * s_on_x_centroid;
  });

  //
  // ===================== Y =====================
  //
  AMREX_HOST_DEVICE_FOR_4D(vbx, ncomp, i, j, k, n,
  {
    Real vpls(0); Real vmns(0);

    if( areafrac_y(i,j,k) > 0 ) {
      if( j <= dom_low.y and
       ugradu_aux::is_equal_to_any(bc_jlo_type(i,dom_low.y-1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sy(i,j,k,n) = state(i,dom_low.y-1,k,state_comp+n);
      }
      else if( j >= dom_high.y+1 and
       ugradu_aux::is_equal_to_any(bc_jhi_type(i,dom_high.y+1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sy(i,j,k,n) = state(i,dom_high.y+1,k,state_comp+n);
      }
      else {
        vpls = state(i,j  ,k,state_comp+n) - .5*y_slopes(i,j  ,k,slopes_comp+n);
        vmns = state(i,j-1,k,state_comp+n) + .5*y_slopes(i,j-1,k,slopes_comp+n);

        sy(i,j,k,n) = upwind( vmns, vpls, v(i,j,k) );
      }
    }
    else {
        sy(i,j,k,n) = my_huge;
    }
  });

  AMREX_HOST_DEVICE_FOR_4D(vbx, ncomp, i, j, k, n,
  {
    int ii = i + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,0)));
    int kk = k + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,1)));

    Real fracx = (ccm_fab(ii,j-1,k) || ccm_fab(ii,j,k)) ? std::abs(fcy_fab(i,j,k,0)) : 0.0;
    Real fracz = (ccm_fab(i,j-1,kk) || ccm_fab(i,j,kk)) ? std::abs(fcy_fab(i,j,k,1)) : 0.0;

    Real s_on_y_centroid = (1.0-fracx)*(1.0-fracz)*sy(i ,j,k ,n)+
                                fracx *(1.0-fracz)*sy(ii,j,k ,n)+
                                fracz *(1.0-fracx)*sy(i ,j,kk,n)+
                                fracx *     fracz *sy(ii,j,kk,n);

    fy(i,j,k,n) = v(i,j,k) * s_on_y_centroid;

  });

  //
  // ===================== Z =====================
  //
  AMREX_HOST_DEVICE_FOR_4D(wbx, ncomp, i, j, k, n,
  {
    Real w_face(0);
    Real wpls(0); Real wmns(0);

    if( areafrac_z(i,j,k) > 0 ) {
      if( k <= dom_low.z and
       ugradu_aux::is_equal_to_any(bc_klo_type(i,j,dom_low.z-1,0),
                                   bc_types.data(), bc_types.size()))
      {
        sz(i,j,k,n) = state(i,j,dom_low.z-1,state_comp+n);
      }
      else if( k >= dom_high.z+1 and
       ugradu_aux::is_equal_to_any(bc_khi_type(i,j,dom_high.z+1,0),
                                   bc_types.data(), bc_types.size()))
      {
        sz(i,j,k,n) = state(i,j,dom_high.z+1,state_comp+n);
      }
      else {
        wpls = state(i,j,k  ,state_comp+n) - .5*z_slopes(i,j,k  ,slopes_comp+n);
        wmns = state(i,j,k-1,state_comp+n) + .5*z_slopes(i,j,k-1,slopes_comp+n);

        sz(i,j,k,n) = upwind( wmns, wpls, w(i,j,k) );
      }
    }
    else {
        sz(i,j,k,n) = my_huge;
    }
  });

  AMREX_HOST_DEVICE_FOR_4D(wbx, ncomp, i, j, k, n,
  {

    int ii = i + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,0)));
    int jj = j + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,1)));

    Real fracx = (ccm_fab(ii,j,k-1) || ccm_fab(ii,j,k)) ? std::abs(fcz_fab(i,j,k,0)) : 0.0;
    Real fracy = (ccm_fab(i,jj,k-1) || ccm_fab(i,jj,k)) ? std::abs(fcz_fab(i,j,k,1)) : 0.0;

    Real s_on_z_centroid = (1.0-fracx)*(1.0-fracy)*sz(i ,j ,k,n)+
                                fracx *(1.0-fracy)*sz(ii,j ,k,n)+
                                fracy *(1.0-fracx)*sz(i ,jj,k,n)+
                                fracx *     fracy *sz(ii,jj,k,n);

    fz(i,j,k,n) = w(i,j,k) * s_on_z_centroid;
  });

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  const int cyclic_x = geom[0].isPeriodic(0) ? 1 : 0;
  const int cyclic_y = geom[0].isPeriodic(1) ? 1 : 0;
  const int cyclic_z = geom[0].isPeriodic(2) ? 1 : 0;

  // Compute div(tau) with EB algorithm
  compute_divop_conv(bx, *conv_in[lev], *ep_g[lev], conv_comp, ncomp, mfi, fxfab, fyfab, fzfab, 
                     areafrac, facecent, flags, volfrac, bndrycent, domain,
                     cyclic_x, cyclic_y, cyclic_z, dx);

  AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
  {
    const Real coefficient(-1/epsilon_g(i,j,k));
    ugradu(i,j,k,conv_comp+n) *= coefficient; 
  });

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif
}
