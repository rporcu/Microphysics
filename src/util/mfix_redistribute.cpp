#include <param_mod_F.H>
#include <cmath>
#include <limits>
#include <mfix.H>
#include <mfix_util_F.H>

namespace eb_redistribute {


void
mfix_apply_eb_redistribution ( Box& bx,
                               MultiFab& div_mf,
                               MultiFab& divc_mf,
                               MultiFab& ep_g,
                               MFIter* mfi,
                               const int icomp,
                               const int ncomp,
                               const EBCellFlagFab& flags_fab,
                               const MultiFab* volfrac,
                               Box& domain,
                               const int cyclic_x,
                               const int cyclic_y,
                               const int cyclic_z,
                               const Real* dx)
{

    //
    // Check that grid is uniform
    //
    const Real tolerance = std::numeric_limits<Real>::epsilon();

    if((std::abs(dx[0] - dx[1]) > tolerance) or
       (std::abs(dx[0] - dx[2]) > tolerance) or
       (std::abs(dx[1] - dx[2]) > tolerance))
        amrex::Abort("Compute divop(): grid spacing must be uniform");


    const amrex::Dim3 dom_low  = amrex::lbound(domain);
    const amrex::Dim3 dom_high = amrex::ubound(domain);


    //
    // Get array from arguments
    //
    Array4<Real> const& div  = div_mf.array(*mfi);
    Array4<Real> const& divc = divc_mf.array(*mfi);
    Array4<Real> const& ep   = ep_g.array(*mfi);
    auto const&        flags = flags_fab.array();
    auto const&        vfrac = volfrac->array(*mfi);

    //
    // Working arrays
    //
    const Box& grown1_bx = amrex::grow(bx,1);
    const Box& grown2_bx = amrex::grow(bx,2);

    FArrayBox  delm_fab(grown1_bx,ncomp);
    FArrayBox  optmp_fab(grown2_bx,ncomp);
    FArrayBox  mask_fab(grown2_bx);

    Array4<Real> const& optmp = optmp_fab.array();
    Array4<Real> const& mask  = mask_fab.array();
    Array4<Real> const& delm  = delm_fab.array();

    //
    // Array "mask" is used to sever the link to ghost cells when the BCs are not
    // periodic
    // It is set to 1 when a cell can be used in computations, 0 otherwise
    //
    AMREX_FOR_3D(grown2_bx, i, j, k,
    {
        if(((not cyclic_x) and (i < dom_low.x or i > dom_high.x)) or
           ((not cyclic_y) and (j < dom_low.y or j > dom_high.y)) or
           ((not cyclic_z) and (k < dom_low.z or k > dom_high.z)))
            mask(i,j,k) = 0;
        else
            mask(i,j,k) = 1;
    });

    //
    // Init to zero tmp array
    //
    AMREX_FOR_4D(grown2_bx, ncomp, i, j, k, n,
    {
        optmp(i,j,k,n) = 0;
    });

    //
    // Step 2: compute delta M (mass gain or loss) on (lo-1,lo+1)
    //
    AMREX_FOR_4D(grown1_bx, ncomp, i, j, k, n,
    {
        if(flags(i,j,k).isSingleValued())
        {
            Real divnc(0.0);
            Real vtot(0.0);
            Real epvfrac(0.0);

            for(int ii(-1); ii <= 1; ii++)
                for(int jj(-1); jj <= 1; jj++)
                    for(int kk(-1); kk <= 1; kk++)
                        if( (ii != 0 or jj != 0 or kk != 0) and
                            (flags(i,j,k).isConnected({ii,jj,kk}) == 1))
                        {
                            epvfrac = vfrac(i+ii,j+jj,k+kk) * ep(i+ii,j+jj,k+kk)
                                * mask(i+ii,j+jj,k+kk);
                            vtot   += epvfrac;
                            divnc  += epvfrac * divc(i+ii,j+jj,k+kk,n);
                        }

            divnc /= vtot;

            // We need to multiply divc by mask to make sure optmp is zero for cells
            // outside the domain for non-cyclic BCs
            optmp(i,j,k,n) =  (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k,n) * mask(i,j,k));
            delm(i,j,k,n)  = -(    vfrac(i,j,k)) * optmp(i,j,k,n);
        }
        else
        {
            delm(i,j,k,n) = 0;
        }
    });


    //
    // Step 3: redistribute excess/loss of mass
    //
    AMREX_FOR_4D(grown1_bx, ncomp, i, j, k, n,
    {
        if(flags(i,j,k).isSingleValued())
        {
            Real wtot(0.0);

            for(int ii(-1); ii <= 1; ii++)
                for(int jj(-1); jj <= 1; jj++)
                    for(int kk(-1); kk <= 1; kk++)
                        if((ii != 0 or jj != 0 or kk != 0) and
                           (flags(i,j,k).isConnected({ii,jj,kk}) == 1))
                        {
                            wtot += ep(i+ii,j+jj,k+kk) * vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
                        }


            wtot = 1/wtot;


            // Note -- based on testing conservation of a conservatively advected tracer ( (rho T)_t + div (rho U T) = 0)
            //         we do *not* want the epsilon-weighting in defining delm but we *do* need the epsilon weighting below
            // Caveat -- to test in the presence of epsilon != 1 we must turn off the advection of particles or the
            //        sum will change because of the change in epsilon
            for(int ii(-1); ii <= 1; ii++)
                for(int jj(-1); jj <= 1; jj++)
                    for(int kk(-1); kk <= 1; kk++)
                        if((ii != 0 or jj != 0 or kk != 0) and
                           (flags(i,j,k).isConnected({ii,jj,kk}) == 1))
                        {
#ifdef AMREX_USE_CUDA
                            Gpu::Atomic::Add(&optmp(i+ii,j+jj,k+kk,n),
                                             delm(i,j,k,n) * wtot * mask(i+ii,j+jj,k+kk) * ep(i+ii,j+jj,k+kk));
#else
                            optmp(i+ii,j+jj,k+kk,n) += delm(i,j,k,n) * wtot * mask(i+ii,j+jj,k+kk) * ep(i+ii,j+jj,k+kk);
#endif
                        }
        }
    });

    //
    // Resume the correct sign, AKA return the negative
    //
    AMREX_FOR_4D(bx, ncomp, i, j, k, n,
    {
        div(i,j,k,icomp+n) = divc(i,j,k,n) + optmp(i,j,k,n);
    });

    Gpu::synchronize();
}


} // end namespace eb_redistribute

using namespace eb_redistribute;

void
mfix::mfix_redistribute ( int lev,
                          MultiFab& conv_tmp_in,
                          MultiFab& conv_out,
                          int conv_comp,
                          int ncomp )
{
    Box domain(geom[lev].Domain());

    EB_set_covered(conv_tmp_in, covered_val);
    conv_tmp_in.FillBoundary(geom[lev].periodicity());

    // Get EB geometric info
    const amrex::MultiFab*   volfrac;
    volfrac   = &(ebfactory[lev] -> getVolFrac());

    for (MFIter mfi(conv_out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Tilebox
        Box bx = mfi.tilebox ();

        // this is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  conv_fab = static_cast<EBFArrayBox const&>(conv_out[mfi]);
        const EBCellFlagFab&  flags = conv_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
        {
            // If tile is completely covered by EB geometry, set slopes
            // value to some very large number so we know if
            // we accidentally use these covered slopes later in calculations
            setFabVal(conv_out[mfi], get_my_huge(), bx, conv_comp, ncomp);
        }
        else
        {
            // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
            if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
            {
                conv_out[mfi].copy(conv_tmp_in[mfi],conv_comp,0,ncomp);
            }
            else
            {
                const int cyclic_x = geom[0].isPeriodic(0) ? 1 : 0;
                const int cyclic_y = geom[0].isPeriodic(1) ? 1 : 0;
                const int cyclic_z = geom[0].isPeriodic(2) ? 1 : 0;

                // Compute div(tau) with EB algorithm
                mfix_apply_eb_redistribution(bx, conv_out, conv_tmp_in, *ep_g[lev], &mfi,
                                             conv_comp, ncomp, flags, volfrac, domain,
                                             cyclic_x, cyclic_y, cyclic_z,
                                             geom[lev].CellSize());

            }
        }
    }
}
