#include <mfix.H>
#include <mfix_mac_F.H>
#include <mfix_divop_F.H>

#define MY_HUGE 1.e200

namespace ugradu_auxiliary {

//
// Compute upwind non-normal velocity
//
Real
upwind( const Real velocity_minus,
        const Real velocity_plus,
        const Real u_edge)
{
  // Small value to protect against tiny velocities used in upwinding
  const Real small_velocity(1.e-10);

  if(std::abs(u_edge) < small_velocity)
    return .5*(velocity_minus+velocity_plus);

  return u_edge > 0 ? velocity_minus : velocity_plus;
}

} // end namespace ugradu_auxiliary

using namespace ugradu_auxiliary;

//
// Compute the three components of the convection term
//
void
mfix::mfix_compute_ugradu( Box& bx,
                           Vector< std::unique_ptr<MultiFab> >& conv, 
                           Vector< std::unique_ptr<MultiFab> >& vel,
                           MFIter* mfi,
                           Box& domain,
                           const int lev)
{
  const Real* dx = geom[lev].CellSize();
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& ugradu = conv[lev]->array(*mfi); 
  
  Array4<Real> const& velocity = vel[lev]->array(*mfi);
  Array4<Real> const& epsilon_g = ep_g[lev]->array(*mfi);
  
  Array4<Real> const& u = m_u_mac[lev]->array(*mfi);
  Array4<Real> const& v = m_v_mac[lev]->array(*mfi);
  Array4<Real> const& w = m_w_mac[lev]->array(*mfi);

  Array4<Real> const& x_slopes = xslopes[lev]->array(*mfi);
  Array4<Real> const& y_slopes = yslopes[lev]->array(*mfi);
  Array4<Real> const& z_slopes = zslopes[lev]->array(*mfi);

  Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
  Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();
  Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
  Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();
  Array4<int> const& bc_klo_type = bc_klo[lev]->array();
  Array4<int> const& bc_khi_type = bc_khi[lev]->array();

  const Real i_dx(1/dx[0]), i_dy(1/dx[1]), i_dz(1/dx[2]);

  // Vectorize the boundary conditions list in order to use it in lambda
  // functions
  const Vector<int> bc = {bc_list.minf, bc_list.pinf, bc_list.pout};

  AMREX_CUDA_HOST_DEVICE_FOR_3D(bx, i, j, k,
  {
      Real u_w(0);
      Real v_w(0);
      Real w_w(0);
      Real u_e(0);
      Real v_e(0);
      Real w_e(0);
      Real u_s(0);
      Real v_s(0);
      Real w_s(0);
      Real u_n(0);
      Real v_n(0);
      Real w_n(0);
      Real u_b(0);
      Real v_b(0);
      Real w_b(0);
      Real u_t(0);
      Real v_t(0);
      Real w_t(0);
      Real upls(0);
      Real umns(0);
      Real vpls(0);
      Real vmns(0);
      Real wpls(0);
      Real wmns(0);

    //
    // West face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((i == dom_low.x) and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_ilo_type(dom_low.x-1,j,k,0) == c;}))
    {
      u_w = velocity(i-1,j,k,0);
      v_w = velocity(i-1,j,k,1);
      w_w = velocity(i-1,j,k,2);
    }
    else {
      upls = velocity(i  ,j,k,0) - .5*x_slopes(i  ,j,k,0);
      umns = velocity(i-1,j,k,0) + .5*x_slopes(i-1,j,k,0);
      vpls = velocity(i  ,j,k,1) - .5*x_slopes(i  ,j,k,1);
      vmns = velocity(i-1,j,k,1) + .5*x_slopes(i-1,j,k,1);
      wpls = velocity(i  ,j,k,2) - .5*x_slopes(i  ,j,k,2);
      wmns = velocity(i-1,j,k,2) + .5*x_slopes(i-1,j,k,2);

      u_w = upwind( umns, upls, u(i,j,k) );
      v_w = upwind( vmns, vpls, u(i,j,k) );
      w_w = upwind( wmns, wpls, u(i,j,k) );
    }

    //
    // East face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((i == dom_high.x) and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_ihi_type(dom_high.x+1,j,k,0) == c;}))
    {
      u_e = velocity(i+1,j,k,0);
      v_e = velocity(i+1,j,k,1);
      w_e = velocity(i+1,j,k,2);
    }
    else {
      upls = velocity(i+1,j,k,0) - .5*x_slopes(i+1,j,k,0);
      umns = velocity(i  ,j,k,0) + .5*x_slopes(i  ,j,k,0);
      vpls = velocity(i+1,j,k,1) - .5*x_slopes(i+1,j,k,1);
      vmns = velocity(i  ,j,k,1) + .5*x_slopes(i  ,j,k,1);
      wpls = velocity(i+1,j,k,2) - .5*x_slopes(i+1,j,k,2);
      wmns = velocity(i  ,j,k,2) + .5*x_slopes(i  ,j,k,2);

      u_e = upwind( umns, upls, u(i+1,j,k) );
      v_e = upwind( vmns, vpls, u(i+1,j,k) );
      w_e = upwind( wmns, wpls, u(i+1,j,k) );
    }

    //
    // South face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((j == dom_low.y) and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_jlo_type(i,dom_low.y-1,k,0) == c;}))
    {
      u_s = velocity(i,j-1,k,0);
      v_s = velocity(i,j-1,k,1);
      w_s = velocity(i,j-1,k,2);
    }
    else {
      upls = velocity(i,j  ,k,0) - .5*y_slopes(i,j  ,k,0);
      umns = velocity(i,j-1,k,0) + .5*y_slopes(i,j-1,k,0);
      vpls = velocity(i,j  ,k,1) - .5*y_slopes(i,j  ,k,1);
      vmns = velocity(i,j-1,k,1) + .5*y_slopes(i,j-1,k,1);
      wpls = velocity(i,j  ,k,2) - .5*y_slopes(i,j  ,k,2);
      wmns = velocity(i,j-1,k,2) + .5*y_slopes(i,j-1,k,2);

      u_s = upwind( umns, upls, v(i,j,k) );
      v_s = upwind( vmns, vpls, v(i,j,k) );
      w_s = upwind( wmns, wpls, v(i,j,k) );
    }

    //
    // North face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((j == dom_high.y) and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_jhi_type(i,dom_high.y+1,k,0) == c;}))
    {
      u_n = velocity(i,j+1,k,0);
      v_n = velocity(i,j+1,k,1);
      w_n = velocity(i,j+1,k,2);
    }
    else {
      upls = velocity(i,j+1,k,0) - .5*y_slopes(i,j+1,k,0);
      umns = velocity(i,j  ,k,0) + .5*y_slopes(i,j  ,k,0);
      vpls = velocity(i,j+1,k,1) - .5*y_slopes(i,j+1,k,1);
      vmns = velocity(i,j  ,k,1) + .5*y_slopes(i,j  ,k,1);
      wpls = velocity(i,j+1,k,2) - .5*y_slopes(i,j+1,k,2);
      wmns = velocity(i,j  ,k,2) + .5*y_slopes(i,j  ,k,2);

      u_n = upwind( umns, upls, v(i,j+1,k) );
      v_n = upwind( vmns, vpls, v(i,j+1,k) );
      w_n = upwind( wmns, wpls, v(i,j+1,k) );
    }

    //
    // Bottom face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((k == dom_low.z) and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_klo_type(i,j,dom_low.z-1,0) == c;}))
    {
      u_b = velocity(i,j,k-1,0);
      v_b = velocity(i,j,k-1,1);
      w_b = velocity(i,j,k-1,2);
    }
    else {
      upls = velocity(i,j,k  ,0) - .5*z_slopes(i,j,k  ,0);
      umns = velocity(i,j,k-1,0) + .5*z_slopes(i,j,k-1,0);
      vpls = velocity(i,j,k  ,1) - .5*z_slopes(i,j,k  ,1);
      vmns = velocity(i,j,k-1,1) + .5*z_slopes(i,j,k-1,1);
      wpls = velocity(i,j,k  ,2) - .5*z_slopes(i,j,k  ,2);
      wmns = velocity(i,j,k-1,2) + .5*z_slopes(i,j,k-1,2);

      u_b = upwind( umns, upls, w(i,j,k) );
      v_b = upwind( vmns, vpls, w(i,j,k) );
      w_b = upwind( wmns, wpls, w(i,j,k) );
    }

    //
    // Top face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((k == dom_high.z) and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_khi_type(i,j,dom_high.z+1,0) == c;}))
    {
      u_t = velocity(i,j,k+1,0);
      v_t = velocity(i,j,k+1,1);
      w_t = velocity(i,j,k+1,2);
    }
    else {
      upls = velocity(i,j,k+1,0) - .5*z_slopes(i,j,k+1,0);
      umns = velocity(i,j,k  ,0) + .5*z_slopes(i,j,k  ,0);
      vpls = velocity(i,j,k+1,1) - .5*z_slopes(i,j,k+1,1);
      vmns = velocity(i,j,k  ,1) + .5*z_slopes(i,j,k  ,1);
      wpls = velocity(i,j,k+1,2) - .5*z_slopes(i,j,k+1,2);
      wmns = velocity(i,j,k  ,2) + .5*z_slopes(i,j,k  ,2);

      u_t = upwind( umns, upls, w(i,j,k+1) );
      v_t = upwind( vmns, vpls, w(i,j,k+1) );
      w_t = upwind( wmns, wpls, w(i,j,k+1) );
    }

    // Define the convective terms -- conservatively
    //   ugradu = ( div(ep_g u^MAC u^cc) - u^cc div(ep_g u^MAC) ) / ep_g
    Real epu_hi_x = .5*(epsilon_g(i+1,j,k)+epsilon_g(i,j,k)) * u(i+1,j,k);
    Real epu_lo_x = .5*(epsilon_g(i-1,j,k)+epsilon_g(i,j,k)) * u(i  ,j,k);
    Real epv_hi_y = .5*(epsilon_g(i,j+1,k)+epsilon_g(i,j,k)) * v(i,j+1,k);
    Real epv_lo_y = .5*(epsilon_g(i,j-1,k)+epsilon_g(i,j,k)) * v(i,j  ,k);
    Real epw_hi_z = .5*(epsilon_g(i,j,k+1)+epsilon_g(i,j,k)) * w(i,j,k+1);
    Real epw_lo_z = .5*(epsilon_g(i,j,k-1)+epsilon_g(i,j,k)) * w(i,j,k  );

    Real divumac = (epu_hi_x - epu_lo_x) * i_dx +
                   (epv_hi_y - epv_lo_y) * i_dy + 
                   (epw_hi_z - epw_lo_z) * i_dz;

    ugradu(i,j,k,0) = (epu_hi_x*u_e - epu_lo_x*u_w) * i_dx +
                      (epv_hi_y*u_n - epv_lo_y*u_s) * i_dy +
                      (epw_hi_z*u_t - epw_lo_z*u_b) * i_dz -
                      velocity(i,j,k,0)*divumac;

    ugradu(i,j,k,1) = (epu_hi_x*v_e - epu_lo_x*v_w) * i_dx +
                      (epv_hi_y*v_n - epv_lo_y*v_s) * i_dy +
                      (epw_hi_z*v_t - epw_lo_z*v_b) * i_dz -
                      velocity(i,j,k,1)*divumac;

    ugradu(i,j,k,2) = (epu_hi_x*w_e - epu_lo_x*w_w) * i_dx +
                      (epv_hi_y*w_n - epv_lo_y*w_s) * i_dy +
                      (epw_hi_z*w_t - epw_lo_z*w_b) * i_dz -
                      velocity(i,j,k,2)*divumac;

    //
    // Return the negative
    //
    const Real coefficient(-1/epsilon_g(i,j,k));
    ugradu(i,j,k,0) *= coefficient; 
    ugradu(i,j,k,1) *= coefficient; 
    ugradu(i,j,k,2) *= coefficient; 
  });
}


//
// Compute [TODO description]
//
void
mfix::mfix_compute_ugradu_eb( Box& bx,
                              Vector< std::unique_ptr<MultiFab> >& conv, 
                              Vector< std::unique_ptr<MultiFab> >& vel,
                              MFIter* mfi,
                              Array<const MultiCutFab*,AMREX_SPACEDIM>& areafrac,
                              Array<const MultiCutFab*,AMREX_SPACEDIM>& facecent,
                              const amrex::MultiFab* volfrac,
                              const amrex::MultiCutFab* bndrycent,
                              Box& domain,
                              const EBCellFlagFab& flags,
                              const int lev)
{
  AMREX_ASSERT_WITH_MESSAGE(nghost >= 4, "Compute divop(): ng must be >= 4");

  const Real* dx = geom[lev].CellSize();
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& ugradu = conv[lev]->array(*mfi);

  Array4<Real> const& velocity = vel[lev]->array(*mfi);
  Array4<Real> const& epsilon_g = ep_g[lev]->array(*mfi);

  Array4<const Real> const& areafrac_x = areafrac[0]->array(*mfi);
  Array4<const Real> const& areafrac_y = areafrac[1]->array(*mfi);
  Array4<const Real> const& areafrac_z = areafrac[2]->array(*mfi);

  Array4<Real> const& u = m_u_mac[lev]->array(*mfi);
  Array4<Real> const& v = m_v_mac[lev]->array(*mfi);
  Array4<Real> const& w = m_w_mac[lev]->array(*mfi);

  Array4<Real> const& x_slopes = xslopes[lev]->array(*mfi);
  Array4<Real> const& y_slopes = yslopes[lev]->array(*mfi);
  Array4<Real> const& z_slopes = zslopes[lev]->array(*mfi);

  Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
  Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();
  Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
  Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();
  Array4<int> const& bc_klo_type = bc_klo[lev]->array();
  Array4<int> const& bc_khi_type = bc_khi[lev]->array();

  // Number of Halo layers
  const int nh(3);

  Box ubx = amrex::surroundingNodes(amrex::grow(bx,nh),0);
  Box vbx = amrex::surroundingNodes(amrex::grow(bx,nh),1);
  Box wbx = amrex::surroundingNodes(amrex::grow(bx,nh),2);

  const int ncomp(3);
  
  FArrayBox fxfab(ubx, ncomp);
  FArrayBox fyfab(vbx, ncomp);
  FArrayBox fzfab(wbx, ncomp);

  Array4<Real> const& fx = fxfab.array();
  Array4<Real> const& fy = fyfab.array();
  Array4<Real> const& fz = fzfab.array();

  Vector<int> bc = {bc_list.minf, bc_list.pinf, bc_list.pout};

  //
  // First compute the convective fluxes at the face center
  // Do this on ALL faces on the tile, i.e. INCLUDE as many ghost faces as
  // possible
  //

  //
  // ===================== X =====================
  //
  AMREX_CUDA_HOST_DEVICE_FOR_4D(ubx, ncomp, i, j, k, n,
  {
	  Real u_face(0); 
	  Real upls(0); 
	  Real umns(0); 

    if( areafrac_x(i,j,k) > 0 ) {
      if( i <= dom_low.x and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_ilo_type(dom_low.x-1,j,k,0) == c;}))
      {
        u_face = velocity(dom_low.x-1,j,k,n);
      }
      else if( i >= dom_high.x+1 and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_ihi_type(dom_high.x+1,j,k,0) == c;}))
      {
        u_face = velocity(dom_high.x+1,j,k,n);
      }
      else {
        upls = velocity(i  ,j,k,n) - .5*x_slopes(i  ,j,k,n);
        umns = velocity(i-1,j,k,n) + .5*x_slopes(i-1,j,k,n);

        u_face = upwind( umns, upls, u(i,j,k) );
      }
    }
    else {
      u_face = MY_HUGE; 
    }
    fx(i,j,k,n) = .5*(epsilon_g(i-1,j,k)+epsilon_g(i,j,k)) * u(i,j,k) * u_face;
  });

  //
  // ===================== Y =====================
  //
  AMREX_CUDA_HOST_DEVICE_FOR_4D(vbx, ncomp, i, j, k, n,
  {
	  Real v_face(0); 
	  Real vpls(0); 
	  Real vmns(0); 
    if( areafrac_y(i,j,k) > 0 ) {
      if( j <= dom_low.y and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_jlo_type(i,dom_low.y-1,k,0) == c;}))
      {
        v_face = velocity(i,dom_low.y-1,k,n);
      }
      else if( j >= dom_high.y+1 and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_jhi_type(i,dom_high.y+1,k,0) == c;}))
      {
        v_face = velocity(i,dom_high.y+1,k,n);
      }
      else {
        vpls = velocity(i,j  ,k,n) - .5*y_slopes(i,j  ,k,n);
        vmns = velocity(i,j-1,k,n) + .5*y_slopes(i,j-1,k,n);

        v_face = upwind( vmns, vpls, v(i,j,k) );
      }
    }
    else {
      v_face = MY_HUGE;
    }
    fy(i,j,k,n) = .5*(epsilon_g(i,j-1,k)+epsilon_g(i,j,k)) * v(i,j,k) * v_face;
  });

  //
  // ===================== Z =====================
  //
  AMREX_CUDA_HOST_DEVICE_FOR_4D(wbx, ncomp, i, j, k, n,
  {
	  Real w_face(0); 
	  Real wpls(0); 
	  Real wmns(0); 
    if( areafrac_z(i,j,k) > 0 ) {
      if( k <= dom_low.z and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_klo_type(i,j,dom_low.z-1,0) == c;}))
      {
        w_face = velocity(i,j,dom_low.z-1,n);
      }
      else if( k >= dom_high.z+1 and
       std::any_of(bc.begin(),
                   bc.end(),
                   [&](int c){return bc_khi_type(i,j,dom_high.z+1,0) == c;}))
      {
        w_face = velocity(i,j,dom_high.z+1,n);
      }
      else {
        wpls = velocity(i,j,k  ,n) - .5*z_slopes(i,j,k  ,n);
        wmns = velocity(i,j,k-1,n) + .5*z_slopes(i,j,k-1,n);

        w_face = upwind( wmns, wpls, w(i,j,k) );
      }
    }
    else {
      w_face = MY_HUGE;
    }
    fz(i,j,k,n) = .5*(epsilon_g(i,j,k-1)+epsilon_g(i,j,k)) * w(i,j,k) * w_face;
  });

  // Compute div(tau) with EB algorithm
  compute_divop(
      BL_TO_FORTRAN_BOX(bx),
      BL_TO_FORTRAN_ANYD((*conv[lev])[*mfi]),
      BL_TO_FORTRAN_ANYD((*vel[lev])[*mfi]),
      BL_TO_FORTRAN_ANYD(fxfab),
      BL_TO_FORTRAN_ANYD(fyfab),
      BL_TO_FORTRAN_ANYD(fzfab),
      BL_TO_FORTRAN_ANYD((*ep_g[lev])[*mfi]),
      BL_TO_FORTRAN_ANYD((*areafrac[0])[*mfi]),
      BL_TO_FORTRAN_ANYD((*areafrac[1])[*mfi]),
      BL_TO_FORTRAN_ANYD((*areafrac[2])[*mfi]),
      BL_TO_FORTRAN_ANYD((*facecent[0])[*mfi]),
      BL_TO_FORTRAN_ANYD((*facecent[1])[*mfi]),
      BL_TO_FORTRAN_ANYD((*facecent[2])[*mfi]),
      BL_TO_FORTRAN_ANYD(flags),
      BL_TO_FORTRAN_ANYD((*volfrac)[*mfi]),
      BL_TO_FORTRAN_ANYD((*bndrycent)[*mfi]),
      domain.loVect(), domain.hiVect(),
      geom[lev].CellSize(), &nghost);

  AMREX_CUDA_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
  {
    const Real coefficient(-1/epsilon_g(i,j,k));
    ugradu(i,j,k,n) *= coefficient;
  });
}



//
// Compute acc using the vel passed in
//
void
mfix::mfix_compute_ugradu_predictor( Vector< std::unique_ptr<MultiFab> >& conv, 
                                     Vector< std::unique_ptr<MultiFab> >& vel,
                                     Real time)
{
    BL_PROFILE("mfix::mfix_compute_ugradu_predictor");

    amrex::Print() << "In predictor at time " << time << std::endl;

    mfix_compute_MAC_velocity_at_faces( time, vel );

    for (int lev=0; lev < nlev; ++lev)
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
       
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vel[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox ();
            
            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
            const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

            if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
            {
                // If tile is completely covered by EB geometry, set slopes
                // value to some very large number so we know if
                // we accidentaly use these covered slopes later in calculations
                conv[lev] -> setVal( 1.2345e300, bx, 0, 3);
            }
            else
            {
                // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
                if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
                {
                    mfix_compute_ugradu(bx, conv, vel, &mfi, domain, lev);
                }
                else
                {
                    mfix_compute_ugradu_eb(bx, conv, vel, &mfi, areafrac, facecent,
                                           volfrac, bndrycent, domain, flags, lev);
                }
            }
        }

    }
}

//
// Compute acc using the vel passed in
//
void
mfix::mfix_compute_ugradu_corrector( Vector< std::unique_ptr<MultiFab> >& conv, 
                                     Vector< std::unique_ptr<MultiFab> >& vel,
                                     Real time)
{
    BL_PROFILE("mfix::mfix_compute_ugradu_corrector");

    amrex::Print() << "In corrector at time " << time << std::endl;

    mfix_compute_MAC_velocity_at_faces( time, vel );

    for (int lev=0; lev < nlev; ++lev)
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
       
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*vel[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox ();
            
            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
            const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

            if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
            {
                // If tile is completely covered by EB geometry, set slopes
                // value to some very large number so we know if
                // we accidentaly use these covered slopes later in calculations
                conv[lev] -> setVal( 1.2345e300, bx, 0, 3);
            }
            else
            {
                // No cut cells in tile grown by nghost -> use non-eb routine
                if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
                {
                    mfix_compute_ugradu(bx, conv, vel, &mfi, domain, lev);
                }
                else
                {
                    mfix_compute_ugradu_eb(bx, conv, vel, &mfi, areafrac, facecent,
                                           volfrac, bndrycent, domain, flags, lev);
                }
            }
        }

    }
}
