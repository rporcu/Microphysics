#include <mfix_divop.hpp>

#include <cmath>
#include <limits>

#define MY_HUGE 1.e200

namespace divop_aux {

AMREX_GPU_HOST_DEVICE
Real 
interp_to_face_centroid_x(const int i,
                          const int j,
                          const int k,
                          Array4<Real> const& var,
                          const int n,
                          Array4<const Real> const& afrac,
                          Array4<const Real> const& cent,
                          EBCellFlag& flag)
{
  Real result(0);

  Real frac_y(0), frac_z(0);

  if (afrac(i,j,k) == 0)
    result = 0;
  else if (afrac(i,j,k) == 1)
    result = var(i,j,k,n);
  else
  {
    if (cent(i,j,k,0) < 0)
    {
      frac_y = -1 * cent(i,j,k,0) * flag.isConnected({AMREX_D_DECL(0,-1,0)});
      if (cent(i,j,k,1) <= 0)
      {
        frac_z = -1 * cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,0,-1)});
        result = (1-frac_z) * (   frac_y * var(i,j-1,k  ,n)  +
                              (1-frac_y) * var(i,j  ,k  ,n)) +
                     frac_z * (   frac_y * var(i,j-1,k-1,n)  +
                              (1-frac_y) * var(i,j  ,k-1,n));
      }
      else
      {
        frac_z = cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,0,1)});
        result = (1-frac_z) * (   frac_y * var(i,j-1,k  ,n)  +
                              (1-frac_y) * var(i,j  ,k  ,n)) +
                     frac_z * (   frac_y * var(i,j-1,k+1,n)  +
                              (1-frac_y) * var(i,j  ,k+1,n));
      }
    }
    else
    {
      frac_y = cent(i,j,k,0) * flag.isConnected({AMREX_D_DECL(0,1,0)});
      if (cent(i,j,k,1) <= 0)
      {
        frac_z = -1 * cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,0,-1)});
        result = (1-frac_z) * (   frac_y * var(i,j+1,k  ,n)  +
                              (1-frac_y) * var(i,j  ,k  ,n)) +
                     frac_z * (   frac_y * var(i,j+1,k-1,n)  +
                              (1-frac_y) * var(i,j  ,k-1,n));
      }
      else
      {
        frac_z = cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,0,1)});
        result = (1-frac_z) * (   frac_y * var(i,j+1,k  ,n)  +
                              (1-frac_y) * var(i,j  ,k  ,n)) +
                     frac_z * (   frac_y * var(i,j+1,k+1,n)  +
                              (1-frac_y) * var(i,j  ,k+1,n));
      }
    }
  }

  return result;
}

AMREX_GPU_HOST_DEVICE
Real 
interp_to_face_centroid_y(const int i,
                          const int j,
                          const int k,
                          Array4<Real> const& var,
                          const int n,
                          Array4<const Real> const& afrac,
                          Array4<const Real> const& cent,
                          EBCellFlag& flag)
{
  Real result(0);

  Real frac_x(0), frac_z(0);

  if (afrac(i,j,k) == 0)
    result = 0;
  else if (afrac(i,j,k) == 1)
    result = var(i,j,k,n);
  else
  {
    if (cent(i,j,k,0) < 0)
    {
      frac_x = -1 * cent(i,j,k,0) * flag.isConnected({AMREX_D_DECL(-1,0,0)});
      if (cent(i,j,k,1) <= 0)
      {
        frac_z = -1 * cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,0,-1)});
        result = (1-frac_z) * (   frac_x * var(i-1,j,k  ,n)  +
                              (1-frac_x) * var(i  ,j,k  ,n)) +
                     frac_z * (   frac_x * var(i-1,j,k-1,n)  +
                              (1-frac_x) * var(i  ,j,k-1,n));
      }
      else
      {
         frac_z =  cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,0,1)});
         result = (1-frac_z) * (   frac_x * var(i-1,j,k  ,n)  +
                               (1-frac_x) * var(i  ,j,k  ,n)) +
                      frac_z * (   frac_x * var(i-1,j,k+1,n)  +
                               (1-frac_x) * var(i  ,j,k+1,n));
      }
    }
    else
    {
      frac_x = cent(i,j,k,0) * flag.isConnected({AMREX_D_DECL(1,0,0)});
      if (cent(i,j,k,1) <= 0)
      {
        frac_z = -1 * cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,0,-1)});
        result = (1-frac_z) * (   frac_x * var(i+1,j,k  ,n)  +
                              (1-frac_x) * var(i  ,j,k  ,n)) +
                     frac_z * (   frac_x * var(i+1,j,k-1,n)  +
                              (1-frac_x) * var(i  ,j,k-1,n));
      }
      else
      {
        frac_z =  cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,0,1)});
        result = (1-frac_z) * (   frac_x * var(i+1,j,k  ,n)  +
                              (1-frac_x) * var(i  ,j,k  ,n)) +
                     frac_z * (   frac_x * var(i+1,j,k+1,n)  +
                              (1-frac_x) * var(i  ,j,k+1,n));
      }
    } 
  }

  return result;
}

AMREX_GPU_HOST_DEVICE
Real 
interp_to_face_centroid_z(const int i,
                          const int j,
                          const int k,
                          Array4<Real> const& var,
                          const int n,
                          Array4<const Real> const& afrac,
                          Array4<const Real> const& cent,
                          EBCellFlag& flag)
{
  Real result(0);

  Real frac_x(0), frac_y(0);

  if (afrac(i,j,k) == 0)
    result = 0;
  else if (afrac(i,j,k) == 1)
    result = var(i,j,k,n);
  else
  {
    if (cent(i,j,k,0) < 0)
    {
      frac_x = -1 * cent(i,j,k,0) * flag.isConnected({AMREX_D_DECL(-1,0,0)});
      if (cent(i,j,k,1) <= 0)
      {
        frac_y = -1 * cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,-1,0)});
        result = (1-frac_y) * (   frac_x * var(i-1,j  ,k,n)  +
                              (1-frac_x) * var(i  ,j  ,k,n)) +
                     frac_y * (   frac_x * var(i-1,j-1,k,n)  +
                              (1-frac_x) * var(i  ,j-1,k,n));
      }
      else
      {
        frac_y =  cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,1,0)});
        result = (1-frac_y) * (   frac_x * var(i-1,j  ,k,n)  +
                              (1-frac_x) * var(i  ,j  ,k,n)) +
                     frac_y * (   frac_x * var(i-1,j+1,k,n)  +
                              (1-frac_x) * var(i  ,j+1,k,n));
      }
    }
    else
    {
      frac_x = cent(i,j,k,0) * flag.isConnected({AMREX_D_DECL(1,0,0)});
      if (cent(i,j,k,1) <= 0)
      {
        frac_y = -1 * cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,-1,0)});
        result = (1-frac_y) * (   frac_x * var(i+1,j  ,k,n)  +
                              (1-frac_x) * var(i  ,j  ,k,n)) +
                     frac_y * (   frac_x * var(i+1,j-1,k,n)  +
                              (1-frac_x) * var(i  ,j-1,k,n));
      }
      else
      {
        frac_y =  cent(i,j,k,1) * flag.isConnected({AMREX_D_DECL(0,1,0)});
        result = (1-frac_y) * (   frac_x * var(i+1,j  ,k,n)  +
                              (1-frac_x) * var(i  ,j  ,k,n)) +
                     frac_y * (   frac_x * var(i+1,j+1,k,n)  +
                              (1-frac_x) * var(i  ,j+1,k,n));
      }
    }
  }

  return result;
}

AMREX_GPU_HOST_DEVICE
void compute_dphidn_3d(Real* dphidn,
                       const int i,
                       const int j,
                       const int k,
                       Array4<Real> const& velocity,
                       Array4<const Real> const& bndrycent,
                       const Real* phib,
                       const Real anrmx,
                       const Real anrmy,
                       const Real anrmz,
                       Array4<const Real> const& vfrac)
{
  Real bct[3] = {bndrycent(i,j,k,0), bndrycent(i,j,k,1), bndrycent(i,j,k,2)};

  Real vf = vfrac(i,j,k);

  //Real dx_eb = amrex_get_dx_eb(vf);
  Real dx_eb = std::max(0.3, (vf*vf - 0.25)/(2.*vf));

  const double coeff = std::max(abs(anrmx), std::max(abs(anrmy),abs(anrmz)));

  Real dg = dx_eb / coeff;

  Real gx = bct[0] - dg*anrmx;
  Real gy = bct[1] - dg*anrmy;
  Real gz = bct[2] - dg*anrmz;

  double sx =  std::copysign(1., anrmx);
  double sy =  std::copysign(1., anrmy);
  double sz =  std::copysign(1., anrmz);
  
  gx *= sx;
  gy *= sy;
  gz *= sz;

  int ii = i - int(sx);
  int jj = j - int(sy);
  int kk = k - int(sz);

  Real gxy = gx*gy;
  Real gxz = gx*gz;
  Real gyz = gy*gz;

  Real gxyz = gx*gy*gz;

  for(unsigned int n(0); n < 3; ++n)
  {
    Real phig(0);

    phig = (1 + gx + gy + gz + gxy + gxz + gyz + gxyz) * velocity(i ,j ,k ,n) -
           (gz + gxz + gyz + gxyz)                     * velocity(i ,j ,kk,n) -
           (gy + gxy + gyz + gxyz)                     * velocity(i ,jj,k ,n) +
           (gyz + gxyz)                                * velocity(i ,jj,kk,n) -
           (gx + gxy + gxz + gxyz)                     * velocity(ii,j ,k ,n) +
           (gxz + gxyz)                                * velocity(ii,j ,kk,n) +
           (gxy + gxyz)                                * velocity(ii,jj,k ,n) -
           (gxyz)                                      * velocity(ii,jj,kk,n);

    dphidn[n] = (phib[n] - phig)/dg;
  }
}

AMREX_GPU_HOST_DEVICE
void compute_diff_wallfluxes(Real* divw,
                             const Real* dx,
                             const int i,
                             const int j,
                             const int k,
                             Array4<Real> const& velocity,
                             const Array4<Real>* mu,
                             Array4<const Real> const& bndrycent,
                             Array4<const EBCellFlag> const& flags,
                             Array4<const Real> const& afrac_x,
                             Array4<const Real> const& afrac_y,
                             Array4<const Real> const& afrac_z,
                             Array4<const Real> const& vfrac,
                             const int* do_explicit_diffusion)
{
  const Real dxinv[3] = {1/dx[0], 1/dx[1], 1/dx[2]};

  Real dapx = afrac_x(i+1,j,k) - afrac_x(i,j,k);
  Real dapy = afrac_y(i,j+1,k) - afrac_y(i,j,k);
  Real dapz = afrac_z(i,j,k+1) - afrac_z(i,j,k);

  Real apnorm = std::sqrt(dapx*dapx + dapy*dapy + dapz*dapz);

  const Real tolerance = std::numeric_limits<Real>::epsilon();

  //if (apnorm == 0)
  if (apnorm <= tolerance)
    amrex::Abort("compute_diff_wallflux: we are in trouble.");

  Real apnorminv = 1/apnorm;
  Real anrmx = -1 * dapx * apnorminv; // unit vector pointing toward the wall
  Real anrmy = -1 * dapy * apnorminv;
  Real anrmz = -1 * dapz * apnorminv;

  // Value on wall -- here we enforce no-slip therefore 0 for all components
  const Real phib[3] = {0, 0, 0};

  Real dveldn[3] = {0, 0, 0};

  compute_dphidn_3d(dveldn, i, j, k, velocity, bndrycent, phib,
                    anrmx, anrmy, anrmz, vfrac);

  // transform them to d/dx, d/dy and d/dz given transverse derivatives are zero
  Real dudx = dveldn[0] * anrmx;
  Real dudy = dveldn[0] * anrmy;
  Real dudz = dveldn[0] * anrmz;

  Real dvdx = dveldn[1] * anrmx;
  Real dvdy = dveldn[1] * anrmy;
  Real dvdz = dveldn[1] * anrmz;

  Real dwdx = dveldn[2] * anrmx;
  Real dwdy = dveldn[2] * anrmy;
  Real dwdz = dveldn[2] * anrmz;

  Real divu = dudx + dvdy + dwdz;
  Real tautmp = -(2./3) * (*mu)(i,j,k) * divu; // This MUST be verified

  Real tauxx = (*mu)(i,j,k) * (dudx + dudx) + tautmp;
  Real tauxy = (*mu)(i,j,k) * (dudy + dvdx);
  Real tauxz = (*mu)(i,j,k) * (dudz + dwdx);

  Real tauyx = (*mu)(i,j,k) * (dvdx + dudy);
  Real tauyy = (*mu)(i,j,k) * (dvdy + dvdy) + tautmp;
  Real tauyz = (*mu)(i,j,k) * (dvdz + dwdy);

  Real tauzx = (*mu)(i,j,k) * (dwdx + dudz);
  Real tauzy = (*mu)(i,j,k) * (dwdy + dvdz);
  Real tauzz = (*mu)(i,j,k) * (dwdz + dwdz) + tautmp;

  if(*do_explicit_diffusion == 0)
  {
    //
    // Subtract diagonal terms of stress tensor, to be obtained through
    // implicit solve instead.
    //
    tauxx -= (*mu)(i,j,k) * dudx;
    tauxy -= (*mu)(i,j,k) * dudy;
    tauxz -= (*mu)(i,j,k) * dudz;

    tauyx -= (*mu)(i,j,k) * dvdx;
    tauyy -= (*mu)(i,j,k) * dvdy;
    tauyz -= (*mu)(i,j,k) * dvdz;

    tauzx -= (*mu)(i,j,k) * dwdx;
    tauzy -= (*mu)(i,j,k) * dwdy;
    tauzz -= (*mu)(i,j,k) * dwdz;
  }

  divw[0] = dxinv[0] * (dapx*tauxx + dapy*tauxy + dapz*tauxz);
  divw[1] = dxinv[1] * (dapx*tauyx + dapy*tauyy + dapz*tauyz);
  divw[2] = dxinv[2] * (dapx*tauzx + dapy*tauzy + dapz*tauzz);
}

void
step2(const Box& grown1_bx,
      const Box& grown2_bx,
      Array4<Real> const& optmp,
      Array4<Real> const& ep_g,
      Array4<Real> const& divc,
      Array4<Real> const& delm,
      Array4<const Real> const& vfrac,
      Array4<Real> const& mask,
      Array4<const EBCellFlag> const& flags)
{
  // TODO isn't it already initialized with zeroes?
  AMREX_HOST_DEVICE_FOR_3D(grown2_bx, i, j, k,
  {
    optmp(i,j,k) = 0;
  });

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  AMREX_HOST_DEVICE_FOR_3D(grown1_bx, i, j, k,
  {
    if(flags(i,j,k).isSingleValued())
    {
      Real divnc = 0;
      Real vtot = 0;

      Real epvfrac = 0;

      // TODO unroll this
      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            // Check if we have to include also cell (i,j,k) itself
            if((ii != 0 or jj != 0 or kk != 0) and 
                (flags(i,j,k).isConnected({AMREX_D_DECL(ii,jj,kk)}) == 1))
            {
              epvfrac = vfrac(i+ii,j+jj,k+kk) * ep_g(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
              vtot += epvfrac;
              divnc += epvfrac * divc(i+ii,j+jj,k+kk);
            }

      divnc /= vtot;
      epvfrac = vfrac(i,j,k) * ep_g(i,j,k);
      optmp(i,j,k) = (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k));
      delm(i,j,k) = -1 * epvfrac * optmp(i,j,k);
    }
    else
      delm(i,j,k) = 0;
  });
}

void
step3(const Box& grown1_bx,
      Array4<Real> const& optmp,
      Array4<Real> const& ep_g,
      Array4<Real> const& delm,
      Array4<const Real> const& vfrac,
      Array4<Real> const& mask,
      Array4<const EBCellFlag> const& flags)
{
  AMREX_HOST_DEVICE_FOR_3D(grown1_bx, i, j, k,
  {
    if(flags(i,j,k).isSingleValued())
    {
      Real wtot = 0;
      
      // TODO unroll this
      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            // Check if we have to include also cell (i,j,k) itself
            if((ii != 0 or jj != 0 or kk != 0) and
                (flags(i,j,k).isConnected({AMREX_D_DECL(ii,jj,kk)}) == 1))
            {
              wtot += ep_g(i+ii,j+jj,k+kk) * vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
            }

      wtot = 1/wtot;
      
      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            // Check if we have to include also cell (i,j,k) itself
            if((ii != 0 or jj != 0 or kk != 0) and
                (flags(i,j,k).isConnected({AMREX_D_DECL(ii,jj,kk)}) == 1))
            {
#ifdef AMREX_USE_CUDA
              Gpu::Atomic::Add(&optmp(i+ii,j+jj,k+kk),
                                delm(i,j,k) * wtot * mask(i+ii,j+jj,k+kk));
#else
              optmp(i+ii,j+jj,k+kk) += 
                delm(i,j,k) * wtot * mask(i+ii,j+jj,k+kk);
#endif
            }
    }
  });
}

} // end namespace amrex::divop_aux

using namespace divop_aux;

void
compute_divop(Box& bx,
              Array4<Real> const& divergence,
              Array4<Real> const& velocity,
              Array4<Real> const& fx,
              Array4<Real> const& fy,
              Array4<Real> const& fz,
              Array4<Real> const& ep_g,
              MFIter* mfi,
              Array<const MultiCutFab*, AMREX_SPACEDIM>& areafrac,
              Array<const MultiCutFab*, AMREX_SPACEDIM>& facecent,
              const EBCellFlagFab& flags_fab,
              const MultiFab* volfrac,
              const MultiCutFab* bndrycent_fab,
              const int cyclic_x,
              const int cyclic_y,
              const int cyclic_z,
              Box& domain,
              const Real* dx,
              const int* nghost,
              const Array4<Real>* mu, 
              const int* do_explicit_diffusion)
{
  // Check number of ghost cells
  AMREX_ASSERT_WITH_MESSAGE(*nghost >= 4, "Compute divop(): ng must be >=4");

  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  const Real i_dx (1/dx[0]), i_dy(1/dx[1]), i_dz(1/dx[2]);

  Array4<const Real> const& areafrac_x = areafrac[0]->array(*mfi);
  Array4<const Real> const& areafrac_y = areafrac[1]->array(*mfi);
  Array4<const Real> const& areafrac_z = areafrac[2]->array(*mfi);

  Array4<const Real> const& facecent_x = facecent[0]->array(*mfi);
  Array4<const Real> const& facecent_y = facecent[1]->array(*mfi);
  Array4<const Real> const& facecent_z = facecent[2]->array(*mfi);

  Array4<const EBCellFlag> const& flags = flags_fab.array();

  Array4<const Real> const& vfrac = volfrac->array(*mfi);
  Array4<const Real> const& bndrycent = bndrycent_fab->array(*mfi);

  bool is_dirichlet(false);

  // Check if we are computing divergence for viscous term by checking if mu was
  // passed in
  if(mu != nullptr)
    is_dirichlet = true;

  const Real tolerance = std::numeric_limits<Real>::epsilon();

  if((std::abs(dx[0] - dx[1]) > tolerance) or 
     (std::abs(dx[0] - dx[2]) > tolerance) or
     (std::abs(dx[1] - dx[2]) > tolerance))
    amrex::Abort("Compute divop(): grid spacing must be uniform");

  const Box& grown1_bx = amrex::grow(bx,1);
  const Box& grown2_bx = amrex::grow(bx,2);

  FArrayBox divc_fbx(grown2_bx);
  FArrayBox optmp_fbx(grown2_bx);
  FArrayBox delm_fbx(grown2_bx);
  FArrayBox mask_fbx(grown2_bx);

  Array4<Real> const& divc = divc_fbx.array();
  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& delm = delm_fbx.array();
  Array4<Real> const& mask = mask_fbx.array();

  //
  // Array "mask" is used to sever the link to ghost cells when the BCs are not
  // periodic
  // It is set to 1 when a cell can be used in computations, 0 otherwise
  //
  AMREX_HOST_DEVICE_FOR_3D(grown2_bx, i, j, k,
  {
    if(((not cyclic_x) and (i < dom_low.x or i > dom_high.x)) or
       ((not cyclic_y) and (j < dom_low.y or j > dom_high.y)) or
       ((not cyclic_z) and (k < dom_low.z or k > dom_high.z)))
      mask(i,j,k) = 0;
    else
      mask(i,j,k) = 1;
  });

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  //
  // We use the EB algorithm to compute the divergence at cell centers
  //
  for(unsigned int n(0); n < 3; ++n)
  {
    //
    // Step 1: compute conservative divergence on stencil (lo-2,hi-2)
    //

    AMREX_HOST_DEVICE_FOR_3D(grown2_bx, i, j, k,
    {
      if(flags(i,j,k).isCovered())
        divc(i,j,k) = MY_HUGE;
      else if(flags(i,j,k).isSingleValued())
      {
        EBCellFlag flag = flags(i,j,k);

        Real fxp = interp_to_face_centroid_x(i+1, j, k, fx, n,
                                             areafrac_x, facecent_x, flag);
        Real fxm = interp_to_face_centroid_x(i  , j, k, fx, n,
                                             areafrac_x, facecent_x, flag);

        Real fyp = interp_to_face_centroid_y(i, j+1, k, fy, n,
                                             areafrac_y, facecent_y, flag);
        Real fym = interp_to_face_centroid_y(i, j  , k, fy, n,
                                             areafrac_y, facecent_y, flag);
        
        Real fzp = interp_to_face_centroid_z(i, j, k+1, fz, n,
                                             areafrac_z, facecent_z, flag);
        Real fzm = interp_to_face_centroid_z(i, j, k  , fz, n,
                                             areafrac_z, facecent_z, flag);

        divc(i,j,k) = ((fxp*areafrac_x(i+1,j,k) - (fxm*areafrac_x(i,j,k))) * i_dx +
                       (fyp*areafrac_y(i,j+1,k) - (fym*areafrac_y(i,j,k))) * i_dy +
                       (fzp*areafrac_z(i,j,k+1) - (fzm*areafrac_z(i,j,k))) * i_dz) / vfrac(i,j,k);

        // Add viscous wall fluxes (compute three components only during the
        // first pass, i.e. for n=0)

        if(is_dirichlet)
        {
          Real divdiff_w[3];

          compute_diff_wallfluxes(divdiff_w, dx, i, j, k, velocity, mu,
                                  bndrycent, flags, areafrac_x, areafrac_y, areafrac_z,
                                  vfrac, do_explicit_diffusion);

          divc(i,j,k) -= divdiff_w[n] / (dx[n]*vfrac(i,j,k));
        }
      }
      else
      {
        divc(i,j,k) = ((fx(i+1,j,k,n) - fx(i,j,k,n)) * i_dx +
                       (fy(i,j+1,k,n) - fy(i,j,k,n)) * i_dy +
                       (fz(i,j,k+1,n) - fz(i,j,k,n)) * i_dz);
      }
    });

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

    //
    // Step 2: compute delta M (mass gain or loss) on (lo-1,lo+1)
    //
    step2(grown1_bx, grown2_bx, optmp, ep_g, divc, delm, vfrac, mask, flags);

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

    //
    // Step 3: redistribute excess/loss of mass
    //
    step3(grown1_bx, optmp, ep_g, delm, vfrac, mask, flags);

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

    //
    // Resume the correct sign, AKA return the negative
    //
    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
    {
      divergence(i,j,k,n) = divc(i,j,k) + optmp(i,j,k);
    });

  }
}
