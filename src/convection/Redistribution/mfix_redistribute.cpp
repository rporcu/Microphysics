#include <Redistribution.H>
#include <AMReX_EB_utils.H>

using namespace amrex;

void redistribution::redistribute_eb (amrex::Box const& bx, int ncomp, int icomp,
                                      amrex::Array4<amrex::Real      > const& dUdt_out,
                                      amrex::Array4<amrex::Real      > const& dUdt_in,
                                      amrex::Array4<amrex::Real const> const& U_in,
                                      amrex::Array4<amrex::Real      > const& scratch,
                                      amrex::Array4<amrex::Real const> const& ep_g,
                                      amrex::Array4<amrex::Real const> const& umac,
                                      amrex::Array4<amrex::Real const> const& vmac,
                                      amrex::Array4<amrex::Real const> const& wmac,
                                      amrex::Array4<amrex::EBCellFlag const> const& flag,
                                      amrex::Array4<amrex::Real const> const& apx,
                                      amrex::Array4<amrex::Real const> const& apy,
                                      amrex::Array4<amrex::Real const> const& apz,
                                      amrex::Array4<amrex::Real const> const& vfrac,
                                      amrex::Array4<amrex::Real const> const& fcx,
                                      amrex::Array4<amrex::Real const> const& fcy,
                                      amrex::Array4<amrex::Real const> const& fcz,
                                      amrex::Array4<amrex::Real const> const& ccc,
                                      amrex::Geometry& lev_geom,
                                      amrex::Real dt, std::string redistribution_type)
{

    int redist_type;
    // redistribution_type = "NoRedist";          // no redistribution
    // redistribution_type = "FluxRedist";        // flux_redistribute
    // redistribution_type = "MergeRedistUpdate"; // merge_redistribute update
    // redistribution_type = "MergeRedistFull";   // merge_redistribute update
    // redistribution_type = "StateRedistUpdate"; // state_redistribute update
    // redistribution_type = "StateRedistFull";   // merge_redistribute full

#if (AMREX_SPACEDIM == 2)
    // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,1),4);
#else
    // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,1),8);
#endif

    amrex::Print() << "REDISTRIBUTION TYPE " << redistribution_type << std::endl;

    if (redistribution_type == "FluxRedist")
    {

      amrex::apply_flux_redistribution(bx, dUdt_out, dUdt_in, ep_g, icomp, ncomp,
                                       flag, vfrac, lev_geom);

    } else if (redistribution_type == "MergeRedistUpdate" || redistribution_type == "MergeRedistFull") {

        Array4<int> itr = itracker.array();
        make_itracker(bx,
                      AMREX_D_DECL(apx, apy, apz), vfrac,
                      itr, lev_geom);

        if (redistribution_type == "MergeRedistFull")
        {
            amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    dUdt_in(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n);
                }
            );
        }

        merge_redistribute(bx, ncomp, icomp, dUdt_out, dUdt_in,
                           AMREX_D_DECL(apx, apy, apz), vfrac,
                           itr, lev_geom);

        if (redistribution_type == "MergeRedistFull")
        {
           amrex::ParallelFor(bx, ncomp,
           [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
               {
                   dUdt_out(i,j,k,n+icomp) = (dUdt_out(i,j,k,n+icomp) - U_in(i,j,k,n)) / dt;
               }
           );
        }


    } else if (redistribution_type == "StateRedistUpdate" || redistribution_type == "StateRedistFull") {
        if (redistribution_type == "StateRedistFull")
        {
            amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    dUdt_in(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n);
                }
            );
        }

        state_redistribute(bx, ncomp, icomp, dUdt_out, dUdt_in, flag,
                           AMREX_D_DECL(apx, apy, apz), vfrac,
                           AMREX_D_DECL(fcx, fcy, fcz), ccc, lev_geom);

        if (redistribution_type == "StateRedistFull")
        {
            amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    dUdt_out(i,j,k,n+icomp) = (dUdt_out(i,j,k,n+icomp) - U_in(i,j,k,n)) / dt;
                }
            );
        }

    } else if (redistribution_type == "NoRedist") {
        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dUdt_out(i,j,k,n+icomp) = dUdt_in(i,j,k,n);
            }
        );

    } else {
       amrex::Error("Not a legit redist_type");
    }
}

void redistribution::redistribute_initial_data (Box const& bx, int ncomp, int icomp,
                                                Array4<Real      > const& U_inout,
                                                Array4<EBCellFlag const> const& flag,
                                                AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                                             amrex::Array4<amrex::Real const> const& apy,
                                                             amrex::Array4<amrex::Real const> const& apz),
                                                amrex::Array4<amrex::Real const> const& vfrac,
                                                AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                                             amrex::Array4<amrex::Real const> const& fcy,
                                                             amrex::Array4<amrex::Real const> const& fcz),
                                                amrex::Array4<amrex::Real const> const& ccc,
                                                Geometry& lev_geom, std::string redistribution_type)
{
#if 0

  int redist_type;
    // redistribution_type = "MergeRedistFull";   // merge_redistribute update
    // redistribution_type = "StateRedistFull";   // merge_redistribute full

#if (AMREX_SPACEDIM == 2)
    // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,1),4);
#else
    // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,1),8);
#endif

    amrex::Print() << "REDISTRIBUTION TYPE " << redistribution_type << std::endl;

    FArrayBox U_out(grow(bx,2),ncomp);
    Array4<Real> uout_array = U_out.array();
    amrex::ParallelFor(Box(uout_array),ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        uout_array(i,j,k,n) = 0.;
    });

    if (redistribution_type == "MergeRedistFull") {
        Array4<int> itr = itracker.array();
        make_itracker(bx,
                      AMREX_D_DECL(apx, apy, apz), vfrac,
                      itr, lev_geom);

        merge_redistribute(bx, ncomp, icomp, uout_array, U_inout,
                           AMREX_D_DECL(apx, apy, apz), vfrac,
                           itr, lev_geom);

    } else if (redistribution_type == "StateRedistFull") {
        state_redistribute(bx, ncomp, icomp, uout_array, U_inout, flag,
                           AMREX_D_DECL(apx, apy, apz), vfrac,
                           AMREX_D_DECL(fcx, fcy, fcz), ccc, lev_geom);

    } else {
       amrex::Error("Should not be here with this redist type");
    }

    // Return the new (redistributed) data in the same U_in that came in
    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        U_inout(i,j,k,n) = uout_array(i,j,k,n);
    });
#endif
}
