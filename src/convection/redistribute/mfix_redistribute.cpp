#include <redistribution.H>
#include <AMReX_IArrayBox.H>

using namespace amrex;

void redistribution::redistribute_eb (Box const& bx, int ncomp, int scomp,
                                      Array4<Real      > const& dUdt_out,
                                      Array4<Real const> const& dUdt_in,
                                      Array4<Real const> const& U_in,
                                      Array4<Real> const& scratch,
                                      amrex::Array4<amrex::Real const> const& umac,
                                      amrex::Array4<amrex::Real const> const& vmac,
                                      amrex::Array4<amrex::Real const> const& wmac,
                                      Array4<EBCellFlag const> const& flag,
                                      amrex::Array4<amrex::Real const> const& apx,
                                      amrex::Array4<amrex::Real const> const& apy,
                                      amrex::Array4<amrex::Real const> const& apz,
                                      amrex::Array4<amrex::Real const> const& vfrac,
                                      amrex::Array4<amrex::Real const> const& fcx,
                                      amrex::Array4<amrex::Real const> const& fcy,
                                      amrex::Array4<amrex::Real const> const& fcz,
                                      amrex::Array4<amrex::Real const> const& ccc,
                                      Geometry& lev_geom, Real dt, std::string advection_type)
{
    int redist_type;
    // redist_type = 0;   // no redistribution
    // redist_type = 1;   // flux_redistribute
    // redist_type = 2;   // state_redistribute
    // redist_type = 3;   // merge_redistribute update
    // redist_type = 4;   // merge_redistribute full

    amrex::Abort("This is not used yet.");
#if 0
    amrex::Print() << "redistribution::redistribute_eb \n";

    if (advection_type == "MOL")
        redist_type = 1;   // flux_redistribute
    else if (advection_type == "Godunov")
        redist_type = 3;   // merge_redistribute update

    // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,1),8);
    itracker.setVal(0);

    if (redist_type == 1)
    {
      flux_redistribute_eb (bx, ncomp, scomp, dUdt_out, dUdt_in, scratch, flag, vfrac, lev_geom);

    } else if (redist_type == 2) {
      state_redistribute_eb(bx, ncomp, scomp, dUdt_out, dUdt_in, flag,
                              AMREX_D_DECL(apx, apy, apz), vfrac,
                              AMREX_D_DECL(fcx, fcy, fcz), ccc, lev_geom);

    } else if (redist_type == 3) {
        Array4<int> itr = itracker.array();
        make_itracker(bx,
                      AMREX_D_DECL(apx, apy, apz), vfrac,
                      itr, lev_geom);

        merge_redistribute_update(bx, ncomp, scomp, dUdt_out, dUdt_in,
                              AMREX_D_DECL(apx, apy, apz), vfrac,
                              itr, lev_geom);

#if 0
    } else if (redist_type == 4) {
      merge_redistribute_full(bx, ncomp, scomp, dUdt_out, dUdt_in, U_in,
                              AMREX_D_DECL(umac, vmac, wmac), flag,
                              AMREX_D_DECL(apx, apy, apz), vfrac,
                              AMREX_D_DECL(fcx, fcy, fcz), ccc, lev_geom, dt);

#endif
    } else if (redist_type == 0) {
        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dUdt_out(i,j,k,n+scomp) = dUdt_in(i,j,k,n+scomp);
            }
        );

    } else {
       amrex::Error("Not a legit redist_type");
    }
#endif
}
