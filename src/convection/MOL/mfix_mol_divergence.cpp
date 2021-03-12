#include <MOL.H>

using namespace amrex;

void
mol::compute_convective_rate (Box const& bx, int ncomp, int icomp,
                              Array4<Real> const& dUdt,
                              Array4<Real const> const& fx,
                              Array4<Real const> const& fy,
                              Array4<Real const> const& fz,
                              Geometry& geom)
{
    const auto dxinv = geom.InvCellSizeArray();
    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        dUdt(i,j,k,n+icomp) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
            +                 dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n))
            +                 dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n));
    });
}

void
mol::compute_convective_rate_eb (Box const& bx, int ncomp, int icomp,
                                 Array4<Real> const& dUdt,
                                 Array4<Real const> const& fx,
                                 Array4<Real const> const& fy,
                                 Array4<Real const> const& fz,
                                 Array4<EBCellFlag const> const& flag,
                                 Array4<Real const> const& vfrac,
                                 Array4<Real const> const& apx,
                                 Array4<Real const> const& apy,
                                 Array4<Real const> const& apz,
                                 Geometry& geom)
{
    const auto dxinv = geom.InvCellSizeArray();
    const Box dbox   = geom.growPeriodicDomain(2);
    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (!dbox.contains(IntVect(i,j,k)) or flag(i,j,k).isCovered()) {
            dUdt(i,j,k,n+icomp) = 0.0;
        } else if (flag(i,j,k).isRegular()) {
            dUdt(i,j,k,n+icomp) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
                +                 dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n))
                +                 dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n));
        } else {
            dUdt(i,j,k,n+icomp) = (1.0/vfrac(i,j,k)) *
                ( dxinv[0] * (apx(i,j,k)*fx(i,j,k,n) - apx(i+1,j,k)*fx(i+1,j,k,n))
                + dxinv[1] * (apy(i,j,k)*fy(i,j,k,n) - apy(i,j+1,k)*fy(i,j+1,k,n))
                + dxinv[2] * (apz(i,j,k)*fz(i,j,k,n) - apz(i,j,k+1)*fz(i,j,k+1,n)) );
        }
    });
}
