#include <AMReX_Box.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuUtility.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiCutFab.H>

#include <MFIX_gradient.H>

#include <cmath>

namespace mfix_aux {

namespace gradient_aux {

AMREX_GPU_HOST_DEVICE
inline
void compute_gradient (amrex::Box const& bx,
                       amrex::Array4<amrex::Real> const& grad,
                       amrex::Array4<amrex::Real const> const& field,
                       amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  const auto lo = lbound(bx);
  const auto hi = ubound(bx);

  const amrex::Real dxi = .5 * dxinv[0];
  const amrex::Real dyi = .5 * dxinv[1];
  const amrex::Real dzi = .5 * dxinv[2];

  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
        grad(i,j,k,0) = dxi * (field(i+1,j,k)-field(i-1,j,k));
        grad(i,j,k,1) = dyi * (field(i,j+1,k)-field(i,j-1,k));
        grad(i,j,k,2) = dzi * (field(i,j,k+1)-field(i,j,k-1));
      }
    }
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_compute_gradient (int i, int j, int k,
                          amrex::Array4<amrex::Real> const& grad,
                          amrex::Array4<amrex::Real const> const& field,
                          amrex::Array4<int const> const& ccm,
                          amrex::Array4<amrex::EBCellFlag const> const& flag,
                          amrex::Array4<amrex::Real const> const& vfrc,
                          amrex::Array4<amrex::Real const> const& apx,
                          amrex::Array4<amrex::Real const> const& apy,
                          amrex::Array4<amrex::Real const> const& apz,
                          amrex::Array4<amrex::Real const> const& fcx,
                          amrex::Array4<amrex::Real const> const& fcy,
                          amrex::Array4<amrex::Real const> const& fcz,
                          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
                          bool already_on_centroids)
{
  const amrex::Real dxi = .5 * dxinv[0];
  const amrex::Real dyi = .5 * dxinv[1];
  const amrex::Real dzi = .5 * dxinv[2];

  if (flag(i,j,k).isCovered())
  {
    grad(i,j,k,0) = 0.0;
    grad(i,j,k,1) = 0.0;
    grad(i,j,k,2) = 0.0;
  }
  else if (flag(i,j,k).isRegular())
  {
    grad(i,j,k,0) = dxi * (field(i+1,j,k)-field(i-1,j,k));
    grad(i,j,k,1) = dyi * (field(i,j+1,k)-field(i,j-1,k));
    grad(i,j,k,2) = dzi * (field(i,j,k+1)-field(i,j,k-1));
  }
  else if (already_on_centroids)
  {
    const amrex::Real vfrc_value = vfrc(i,j,k);

    grad(i,j,k,0) = (dxi/vfrc_value) * (apx(i+1,j,k)*field(i+1,j,k)-apx(i-1,j,k)*field(i-1,j,k));
    grad(i,j,k,1) = (dyi/vfrc_value) * (apy(i,j+1,k)*field(i,j+1,k)-apy(i,j-1,k)*field(i,j-1,k));
    grad(i,j,k,2) = (dzi/vfrc_value) * (apz(i,j,k+1)*field(i,j,k+1)-apz(i,j,k-1)*field(i,j,k-1));
  }
  else
  {
    const amrex::Real vfrc_value = vfrc(i,j,k);

    amrex::Real fxm = field(i-1,j,k);
    amrex::Real fxp = field(i+1,j,k);
    amrex::Real fym = field(i,j-1,k);
    amrex::Real fyp = field(i,j+1,k);
    amrex::Real fzm = field(i,j,k-1);
    amrex::Real fzp = field(i,j,k+1);

    const amrex::Real apx_minus_value = apx(i-1,j,k);
    if (apx_minus_value != 0.0 and apx_minus_value != 1.0)
    {
      int jj = j + static_cast<int>(amrex::Math::copysign(1.0, fcx(i-1,j,k,0)));
      int kk = k + static_cast<int>(amrex::Math::copysign(1.0, fcx(i-1,j,k,1)));

      amrex::Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? amrex::Math::abs(fcx(i-1,j,k,0)) : 0.0;
      amrex::Real fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? amrex::Math::abs(fcx(i-1,j,k,1)) : 0.0;

      fxm = (1.0-fracy)*(1.0-fracz)*fxm
          +      fracy *(1.0-fracz)*field(i-1,jj,k )
          +      fracz *(1.0-fracy)*field(i-1,j ,kk)
          +      fracy *     fracz *field(i-1,jj,kk);
    }

    const amrex::Real apx_plus_value = apx(i+1,j,k);
    if (apx_plus_value != 0.0 and apx_plus_value != 1.0)
    {
      int jj = j + static_cast<int>(amrex::Math::copysign(1.0, fcx(i+1,j,k,0)));
      int kk = k + static_cast<int>(amrex::Math::copysign(1.0, fcx(i+1,j,k,1)));

      amrex::Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? amrex::Math::abs(fcx(i+1,j,k,0)) : 0.0;
      amrex::Real fracz = (ccm(i,j,kk) || ccm(i+1,j,kk)) ? amrex::Math::abs(fcx(i+1,j,k,1)) : 0.0;

      fxp = (1.0-fracy)*(1.0-fracz)*fxp
          +      fracy *(1.0-fracz)*field(i+1,jj,k )
          +      fracz *(1.0-fracy)*field(i+1,j ,kk)
          +      fracy *     fracz *field(i+1,jj,kk);
    }

    const amrex::Real apy_minus_value = apy(i,j-1,k);
    if (apy_minus_value != 0.0 and apy_minus_value != 1.0)
    {
      int ii = i + static_cast<int>(amrex::Math::copysign(1.0, fcy(i,j-1,k,0)));
      int kk = k + static_cast<int>(amrex::Math::copysign(1.0, fcy(i,j-1,k,1)));

      amrex::Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? amrex::Math::abs(fcy(i,j-1,k,0)) : 0.0;
      amrex::Real fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? amrex::Math::abs(fcy(i,j-1,k,1)) : 0.0;

      fym = (1.0-fracx)*(1.0-fracz)*fym
          +      fracx *(1.0-fracz)*field(ii,j-1,k )
          +      fracz *(1.0-fracx)*field(i ,j-1,kk)
          +      fracx *     fracz *field(ii,j-1,kk);
    }

    const amrex::Real apy_plus_value = apy(i,j+1,k);
    if (apy_plus_value != 0.0 and apy_plus_value != 1.0)
    {
      int ii = i + static_cast<int>(amrex::Math::copysign(1.0, fcy(i,j+1,k,0)));
      int kk = k + static_cast<int>(amrex::Math::copysign(1.0, fcy(i,j+1,k,1)));

      amrex::Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? amrex::Math::abs(fcy(i,j+1,k,0)) : 0.0;
      amrex::Real fracz = (ccm(i,j,kk) || ccm(i,j+1,kk)) ? amrex::Math::abs(fcy(i,j+1,k,1)) : 0.0;

      fyp = (1.0-fracx)*(1.0-fracz)*fyp
          +      fracx *(1.0-fracz)*field(ii,j+1,k )
          +      fracz *(1.0-fracx)*field(i ,j+1,kk)
          +      fracx *     fracz *field(ii,j+1,kk);
    }

    const amrex::Real apz_minus_value = apz(i,j,k-1);
    if (apz_minus_value != 0.0 and apz_minus_value != 1.0)
    {
      int ii = i + static_cast<int>(amrex::Math::copysign(1.0, fcz(i,j,k-1,0)));
      int jj = j + static_cast<int>(amrex::Math::copysign(1.0, fcz(i,j,k-1,1)));

      amrex::Real fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? amrex::Math::abs(fcz(i,j,k-1,0)) : 0.0;
      amrex::Real fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? amrex::Math::abs(fcz(i,j,k-1,1)) : 0.0;

      fzm = (1.0-fracx)*(1.0-fracy)*fzm
          +      fracx *(1.0-fracy)*field(ii,j ,k-1)
          +      fracy *(1.0-fracx)*field(i ,jj,k-1)
          +      fracx *     fracy *field(ii,jj,k-1);
    }

    const amrex::Real apz_plus_value = apz(i,j,k+1);
    if (apz_plus_value != 0.0 and apz_plus_value != 1.0)
    {
      int ii = i + static_cast<int>(amrex::Math::copysign(1.0, fcz(i,j,k+1,0)));
      int jj = j + static_cast<int>(amrex::Math::copysign(1.0, fcz(i,j,k+1,1)));

      amrex::Real fracx = (ccm(ii,j,k) || ccm(ii,j,k+1)) ? amrex::Math::abs(fcz(i,j,k+1,0)) : 0.0;
      amrex::Real fracy = (ccm(i,jj,k) || ccm(i,jj,k+1)) ? amrex::Math::abs(fcz(i,j,k+1,1)) : 0.0;

      fzp = (1.0-fracx)*(1.0-fracy)*fzp
          +      fracx *(1.0-fracy)*field(ii,j ,k+1)
          +      fracy *(1.0-fracx)*field(i ,jj,k+1)
          +      fracx *     fracy *field(ii,jj,k+1);
    }

    grad(i,j,k,0) = (dxi/vfrc_value) * (apx_plus_value*fxp - apx_minus_value*fxm);
    grad(i,j,k,1) = (dyi/vfrc_value) * (apy_plus_value*fyp - apy_minus_value*fym);
    grad(i,j,k,2) = (dzi/vfrc_value) * (apz_plus_value*fzp - apz_minus_value*fzm);
  }
}

} // end namespace gradient_aux

using namespace amrex;
using namespace gradient_aux;

void computeGradient (MultiFab& grad,
                      const MultiFab& field,
                      const Geometry& geom)
{
  AMREX_ASSERT(grad.nComp() == AMREX_SPACEDIM);

  const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(grad, TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();

    const auto& gradarr = grad.array(mfi);
    const auto& fieldarr = field.const_array(mfi);

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
    {
      compute_gradient(tbx, gradarr, fieldarr, dxinv);
    });
  }
}

void EB_computeGradient (MultiFab& grad,
                         const MultiFab& field,
                         const Geometry& geom,
                         bool already_on_centroids)
{
  AMREX_ASSERT(grad.nComp() == AMREX_SPACEDIM);
  
  if (not grad.hasEBFabFactory())
  {
    computeGradient(grad, field, geom);
  }
  else
  {
    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(grad.Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();
    const auto& vfrac = factory.getVolFrac();
    const auto& area = factory.getAreaFrac();
    const auto& fcent = factory.getFaceCent();

    iMultiFab cc_mask;
    if (not already_on_centroids)
    {
      cc_mask.define(grad.boxArray(), grad.DistributionMap(), 1, 1);
      cc_mask.BuildMask(geom.Domain(), geom.periodicity(), 1, 0, 0, 1);
    }

    const GpuArray<Real,AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    MFItInfo info;

    if (Gpu::notInLaunchRegion())
      info.EnableTiling().SetDynamic(true);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(grad,info); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      const auto& flagfab = flags[mfi];

      Array4<Real> const& gradarr = grad.array(mfi);
      Array4<Real const> const& fieldarr = field.const_array(mfi);

      const auto fabtyp = flagfab.getType(bx);

      if (fabtyp == FabType::covered)
      {
        AMREX_HOST_DEVICE_FOR_4D(bx, grad.nComp(), i, j, k, n,
            { gradarr(i,j,k,n) = 0.0; });
      }
      else if (fabtyp == FabType::regular)
      {
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, b,
            { compute_gradient(b, gradarr, fieldarr, dxinv); });
      }
      else
      {
        Array4<int const> const& ccm = (already_on_centroids) ?
            Array4<int const>{} : cc_mask.const_array(mfi);

        Array4<Real const> const& vol = vfrac.const_array(mfi);

        Array4<Real const> const& apx = area[0]->const_array(mfi);
        Array4<Real const> const& apy = area[1]->const_array(mfi);
        Array4<Real const> const& apz = area[2]->const_array(mfi);

        Array4<Real const> const& fcx = fcent[0]->const_array(mfi);
        Array4<Real const> const& fcy = fcent[1]->const_array(mfi);
        Array4<Real const> const& fcz = fcent[2]->const_array(mfi);

        Array4<EBCellFlag const> const& flagarr = flagfab.const_array();

        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
        {
          eb_compute_gradient(i, j, k, gradarr, fieldarr,
                              ccm, flagarr, vol, apx, apy, apz,
                              fcx, fcy, fcz, dxinv, already_on_centroids);
        });
      }
    }
  }
}

} // end namespace mfix_aux
