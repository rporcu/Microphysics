#include <mfix.H>
#include <mfix_F.H>
#include <mfix_proj_F.H>
#include <mfix_mac_F.H>

//
// Compute the slopes of each velocity component in all three directions
// 
void
mfix::mfix_compute_velocity_slopes (int lev, Real time, MultiFab& Sborder)
{
    BL_PROFILE("mfix::mfix_compute_velocity_slopes");

    Box domain(geom[lev].Domain());
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Sborder,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       // Tilebox
       Box bx = mfi.tilebox ();

       // this is to check efficiently if this tile contains any eb stuff
       const EBFArrayBox&  v_fab = static_cast<EBFArrayBox const&>(Sborder[mfi]);
       const EBCellFlagFab&  flags = v_fab.getEBCellFlagFab();

       if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
       {
          // If tile is completely covered by EB geometry, set slopes
          // value to some very large number so we know if
          // we accidentaly use these covered slopes later in calculations
          xslopes[lev] -> setVal( 1.2345e300, bx, 0, 3);
          yslopes[lev] -> setVal( 1.2345e300, bx, 0, 3);
          zslopes[lev] -> setVal( 1.2345e300, bx, 0, 3);
       }
       else
       {
          // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
          if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
          {
             const auto& vel_fab =      Sborder.array(mfi);
             const auto&  xs_fab = xslopes[lev]->array(mfi);
             const auto&  ys_fab = yslopes[lev]->array(mfi);
             const auto&  zs_fab = zslopes[lev]->array(mfi);

             int ncomp = Sborder.nComp();
             amrex::ParallelFor(bx, ncomp,
                   [=] (int i, int j, int k, int n)
             {
                // X direction
                Real du_xl = 2.0*(vel_fab(i  ,j,k,n) - vel_fab(i-1,j,k,n));
                Real du_xr = 2.0*(vel_fab(i+1,j,k,n) - vel_fab(i  ,j,k,n));
                Real du_xc = 0.5*(vel_fab(i+1,j,k,n) - vel_fab(i-1,j,k,n));

                Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                xs_fab(i,j,k,n) = (du_xc       > 0.0) ? xslope : -xslope;
 
                // Y direction
                Real du_yl = 2.0*(vel_fab(i,j  ,k,n) - vel_fab(i,j-1,k,n));
                Real du_yr = 2.0*(vel_fab(i,j+1,k,n) - vel_fab(i,j  ,k,n));
                Real du_yc = 0.5*(vel_fab(i,j+1,k,n) - vel_fab(i,j-1,k,n));

                Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                ys_fab(i,j,k,n) = (du_yc       > 0.0) ? yslope : -yslope;
 
                // Z direction
                Real du_zl = 2.0*(vel_fab(i,j,k  ,n) - vel_fab(i,j,k-1,n));
                Real du_zr = 2.0*(vel_fab(i,j,k+1,n) - vel_fab(i,j,k  ,n));
                Real du_zc = 0.5*(vel_fab(i,j,k+1,n) - vel_fab(i,j,k-1,n));

                Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                zs_fab(i,j,k,n) = (du_zc       > 0.0) ? zslope : -zslope;
             });

          }
          else 
          {
             const auto& vel_fab  =      Sborder.array(mfi);
             const auto&  xs_fab  = xslopes[lev]->array(mfi);
             const auto&  ys_fab  = yslopes[lev]->array(mfi);
             const auto&  zs_fab  = zslopes[lev]->array(mfi);
             const auto& flag_fab =         flags.array();

             Real my_huge = 1.e100;

             int ncomp = Sborder.nComp();
             amrex::ParallelFor(bx, ncomp,
                   [=] (int i, int j, int k, int n)
             {
                if (flag_fab(i,j,k).isCovered())
                {
                    xs_fab(i,j,k,n) = 0.0;
                    ys_fab(i,j,k,n) = 0.0;
                    zs_fab(i,j,k,n) = 0.0;

                } else {

                    // X direction
                    Real du_xl = (flag_fab(i-1,j,k).isCovered()) ? 0.0 :
                                 2.0*(vel_fab(i  ,j,k,n) - vel_fab(i-1,j,k,n));
                    Real du_xr = (flag_fab(i+1,j,k).isCovered()) ? 0.0 :
                                 2.0*(vel_fab(i+1,j,k,n) - vel_fab(i  ,j,k,n));
                    Real du_xc = 0.5*(vel_fab(i+1,j,k,n) - vel_fab(i-1,j,k,n));

                    Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                    xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                    xs_fab(i,j,k,n) = (du_xc       > 0.0) ? xslope : -xslope;
 
                    // Y direction
                    Real du_yl = (flag_fab(i,j-1,k).isCovered()) ? 0.0 :
                                 2.0*(vel_fab(i,j  ,k,n) - vel_fab(i,j-1,k,n));
                    Real du_yr = (flag_fab(i,j+1,k).isCovered()) ? 0.0 :
                                 2.0*(vel_fab(i,j+1,k,n) - vel_fab(i,j  ,k,n));
                    Real du_yc = 0.5*(vel_fab(i,j+1,k,n) - vel_fab(i,j-1,k,n));
    
                    Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                    yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                    ys_fab(i,j,k,n) = (du_yc       > 0.0) ? yslope : -yslope;
     
                    // Z direction
                    Real du_zl = (flag_fab(i,j,k-1).isCovered()) ? 0.0 :
                                 2.0*(vel_fab(i,j,k  ,n) - vel_fab(i,j,k-1,n));
                    Real du_zr = (flag_fab(i,j,k+1).isCovered()) ? 0.0 :
                                 2.0*(vel_fab(i,j,k+1,n) - vel_fab(i,j,k  ,n));
                    Real du_zc = 0.5*(vel_fab(i,j,k+1,n) - vel_fab(i,j,k-1,n));
    
                    Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                    zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                    zs_fab(i,j,k,n) = (du_zc       > 0.0) ? zslope : -zslope;
                }
             });
          }

          correct_slopes( BL_TO_FORTRAN_BOX(bx),
                          BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                          (*xslopes[lev])[mfi].dataPtr (),
                          (*yslopes[lev])[mfi].dataPtr (),
                          BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
                          BL_TO_FORTRAN_ANYD(flags),
                          domain.loVect (), domain.hiVect (),
                          bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                          bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                          bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                          &nghost );
       }
    }

    xslopes[lev] -> FillBoundary(geom[lev].periodicity());
    yslopes[lev] -> FillBoundary(geom[lev].periodicity());
    zslopes[lev] -> FillBoundary(geom[lev].periodicity());
}

