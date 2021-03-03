#include <Redistribution.H>
#include <AMReX_EB_utils.H>

using namespace amrex;

void redistribution::redistribute_eb (amrex::Box const& bx, int ncomp, int icomp,
                                      amrex::Array4<amrex::Real      > const& dUdt_out,
                                      amrex::Array4<amrex::Real      > const& dUdt_in,
                                      amrex::Array4<amrex::Real const> const& U_in,
                                      amrex::Array4<amrex::Real      > const& scratch,
                                      amrex::Array4<amrex::Real const> const& ep_g,
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
    // redistribution_type = "NoRedist";      // no redistribution
    // redistribution_type = "FluxRedist"     // flux_redistribute
    // redistribution_type = "MergeRedist";   // merge redistribute
    // redistribution_type = "StateRedist";   // state redistribute

#if (AMREX_SPACEDIM == 2)
    // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,4),4);
#else
    // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,4),8);
#endif

    Array4<int> itr = itracker.array();

    amrex::ParallelFor(bx,ncomp, [=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        dUdt_out(i,j,k,n+icomp) = 0.;
    });

    if (redistribution_type == "FluxRedist")
    {

      amrex::apply_flux_redistribution(bx, dUdt_out, dUdt_in, ep_g, icomp, ncomp,
                                       flag, vfrac, lev_geom);

    } else if (redistribution_type == "MergeRedist") {

      amrex::ParallelFor(Box(dUdt_in), ncomp, [=]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        dUdt_in(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n);
      });

      make_itracker(bx, AMREX_D_DECL(apx, apy, apz), vfrac, itr, lev_geom, "Merge");

      merge_redistribute(bx, ncomp, icomp, dUdt_out, dUdt_in, vfrac, itr, lev_geom);

      amrex::ParallelFor(bx, ncomp, [=]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        dUdt_out(i,j,k,n+icomp) = (dUdt_out(i,j,k,n+icomp) - U_in(i,j,k,n)) / dt;
      });

    } else if (redistribution_type == "StateRedist") {

        Box domain_per_grown = lev_geom.Domain();
        if (lev_geom.isPeriodic(0)) domain_per_grown.grow(0,1);
        if (lev_geom.isPeriodic(1)) domain_per_grown.grow(1,1);
        if (lev_geom.isPeriodic(2)) domain_per_grown.grow(2,1);

        Box const& bxg1 = grow(bx,1);

        // At any external Dirichlet domain boundaries we need to set dUdt_in to 0
        //    in the cells just outside the domain because those values will be used
        //    in the slope computation in state redistribution.  We assume here that
        //    the ext_dir values of U_in itself have already been set.
        if (!domain_per_grown.contains(bxg1)){
            amrex::ParallelFor(bxg1,ncomp, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (!domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k))))
                    dUdt_in(i,j,k,n) = 0.;
            });
        }

        amrex::ParallelFor(Box(dUdt_in), ncomp, [=]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            dUdt_in(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n);
        });

        make_itracker(bx, apx, apy, apz, vfrac, itr, lev_geom, "State");

        state_redistribute(bx, ncomp, icomp, dUdt_out, dUdt_in, flag, vfrac,
                           fcx, fcy, fcz, ccc, itr, lev_geom);

        amrex::ParallelFor(bx, ncomp, [=]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          dUdt_out(i,j,k,n+icomp) = (dUdt_out(i,j,k,n+icomp) - U_in(i,j,k,n)) / dt;
        });


    } else if (redistribution_type == "NoRedist") {
        amrex::ParallelFor(bx, ncomp, [=]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            dUdt_out(i,j,k,n+icomp) = dUdt_in(i,j,k,n);
        });

    } else {
       amrex::Error("Not a legit redist_type");
    }
}

void redistribution::redistribute_initial_data (Box const& bx, int ncomp, int icomp,
                                                Array4<Real      > const& U_out,
                                                Array4<Real      > const& U_in,
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
    // redistribution_type = "MergeRedist";   // merge redistribute
    // redistribution_type = "StateRedist";   // state redistribute

    // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
    //    use the first component of this for the number of neighbors
    IArrayBox itracker(grow(bx,4),8);

    amrex::ParallelFor(bx,ncomp, [=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        U_out(i,j,k,n+icomp) = 0.;
    });

    Array4<int> itr = itracker.array();

    if (redistribution_type == "MergeRedist") {

        make_itracker(bx, apx, apy, apz, vfrac, itr, lev_geom, "Merge");

        merge_redistribute(bx, ncomp, icomp, U_out, U_in, vfrac, itr, lev_geom);

    } else if (redistribution_type == "StateRedist") {

        make_itracker(bx, apx, apy, apz, vfrac, itr, lev_geom, "State");

        state_redistribute(bx, ncomp, icomp, U_out, U_in, flag, vfrac,
                           fcx, fcy, fcz, ccc, itr, lev_geom);

    } else {
       amrex::Error("Shouldn't be here with this redist type");
    }

#endif
}
