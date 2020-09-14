#include <mfix.H>
#include <mfix_algorithm.H>
#include <AMReX_Slopes_K.H>
#include <AMReX_EB_slopes_K.H>


//
// Compute the slopes of Sborder (velocity, density, temperature or tracer)
//
void
mfix::mfix_compute_slopes (int lev,
                           Real time,
                           MultiFab& Sborder,
                           Vector< MultiFab* > const& xslopes_in,
                           Vector< MultiFab* > const& yslopes_in,
                           Vector< MultiFab* > const& zslopes_in,
                           int slopes_comp,
                           std::map<std::string, Gpu::ManagedVector<int>>& bc_types)
{
    BL_PROFILE("mfix::mfix_compute_slopes");

    EB_set_covered(Sborder, 0, Sborder.nComp(), 1, covered_val);

    Box domain(geom[lev].Domain());

    int ncomp = Sborder.nComp();

    // We initialize slopes to zero in the grown domain ... this is essential
    //    to handle the up-winding at outflow faces
    xslopes_in[lev]->setVal(0.0, slopes_comp, ncomp, xslopes_in[lev]->nGrow());
    yslopes_in[lev]->setVal(0.0, slopes_comp, ncomp, yslopes_in[lev]->nGrow());
    zslopes_in[lev]->setVal(0.0, slopes_comp, ncomp, zslopes_in[lev]->nGrow());

    // ... then set them to this large number in the interior in order to be sure
    //     that no "bad" values go unnoticed
    xslopes_in[lev]->setVal(1.2345e300, slopes_comp, ncomp, 0);
    yslopes_in[lev]->setVal(1.2345e300, slopes_comp, ncomp, 0);
    zslopes_in[lev]->setVal(1.2345e300, slopes_comp, ncomp, 0);

    const auto cellcent = &(ebfactory[lev] -> getCentroid());
    const auto fcx = ebfactory[lev] -> getFaceCent()[0];
    const auto fcy = ebfactory[lev] -> getFaceCent()[1];
    const auto fcz = ebfactory[lev] -> getFaceCent()[2];

    const int* bct_Dirichlet = (bc_types["Dirichlet"]).data();
    const int bct_size = (bc_types["Dirichlet"]).size();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Sborder,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       // Tilebox
       Box bx = mfi.tilebox();

       // This is to check efficiently if this tile contains any eb stuff
       const EBFArrayBox& Sborder_fab = static_cast<EBFArrayBox const&>(Sborder[mfi]);
       const EBCellFlagFab& flags = Sborder_fab.getEBCellFlagFab();

       if (flags.getType(amrex::grow(bx,0)) != FabType::covered)
       {
           const auto& state_fab = Sborder.array(mfi);
           const auto& xs_fab = xslopes_in[lev]->array(mfi);
           const auto& ys_fab = yslopes_in[lev]->array(mfi);
           const auto& zs_fab = zslopes_in[lev]->array(mfi);

           const auto& ilo_ifab = bc_ilo[lev]->array();
           const auto& ihi_ifab = bc_ihi[lev]->array();
           const auto& jlo_ifab = bc_jlo[lev]->array();
           const auto& jhi_ifab = bc_jhi[lev]->array();
           const auto& klo_ifab = bc_klo[lev]->array();
           const auto& khi_ifab = bc_khi[lev]->array();

           const int domain_ilo = domain.smallEnd(0);
           const int domain_ihi = domain.bigEnd(0);
           const int domain_jlo = domain.smallEnd(1);
           const int domain_jhi = domain.bigEnd(1);
           const int domain_klo = domain.smallEnd(2);
           const int domain_khi = domain.bigEnd(2);

           // No cut cells in tile + 1-cell width halo -> use non-eb routine
           if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
           {
             amrex::ParallelFor(bx, ncomp,
               [state_fab,xs_fab,ys_fab,zs_fab,slopes_comp,bct_Dirichlet,bct_size,
                ilo_ifab,ihi_ifab,jlo_ifab,jhi_ifab,klo_ifab,khi_ifab,
                domain_ilo,domain_ihi,domain_jlo,domain_jhi,domain_klo,domain_khi,bx]
               AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
               {
                   bool ed_ilo = (domain_ilo >= bx.smallEnd(0) and domain_ilo <= bx.bigEnd(0)) and
                                 (i == domain_ilo) and 
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(ilo_ifab(i-1,j,k,0))));

                   bool ed_ihi = (domain_ihi >= bx.smallEnd(0) and domain_ihi <= bx.bigEnd(0)) and
                                 (i == domain_ihi) and 
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(ihi_ifab(i+1,j,k,0))));

                   bool ed_jlo = (domain_jlo >= bx.smallEnd(1) and domain_jlo <= bx.bigEnd(1)) and
                                 (j == domain_jlo) and
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(jlo_ifab(i,j-1,k,0))));

                   bool ed_jhi = (domain_jhi >= bx.smallEnd(1) and domain_jhi <= bx.bigEnd(1)) and
                                 (j == domain_jhi) and 
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(jhi_ifab(i,j+1,k,0))));

                   bool ed_klo = (domain_klo >= bx.smallEnd(2) and domain_klo <= bx.bigEnd(2)) and
                                 (k == domain_klo) and 
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(klo_ifab(i,j,k-1,0))));

                   bool ed_khi = (domain_khi >= bx.smallEnd(2) and domain_khi <= bx.bigEnd(2)) and
                                 (k == domain_khi) and
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(khi_ifab(i,j,k+1,0))));

                   // X direction
                   xs_fab(i,j,k,slopes_comp+n) = amrex_calc_xslope_extdir(i,j,k,n,2,state_fab,
                                                   ed_ilo,ed_ihi,domain_ilo,domain_ihi);

                   // Y direction
                   ys_fab(i,j,k,slopes_comp+n) = amrex_calc_yslope_extdir(i,j,k,n,2,state_fab,
                                                   ed_jlo,ed_jhi,domain_jlo,domain_jhi);

                   // Z direction
                   zs_fab(i,j,k,slopes_comp+n) = amrex_calc_zslope_extdir(i,j,k,n,2,state_fab,
                                                   ed_klo,ed_khi,domain_klo,domain_khi);
               });
           }
           else
           {
               const auto& flag_fab = flags.array();
               const auto& ccent_fab = cellcent->array(mfi);
               const auto& fcx_fab = fcx->array(mfi);
               const auto& fcy_fab = fcy->array(mfi);
               const auto& fcz_fab = fcz->array(mfi);

               amrex::ParallelFor(bx, ncomp,
               [state_fab,xs_fab,ys_fab,zs_fab,slopes_comp,flag_fab,ccent_fab,bct_Dirichlet,bct_size,
                fcx_fab,fcy_fab,fcz_fab,
                ilo_ifab,ihi_ifab,jlo_ifab,jhi_ifab,klo_ifab,khi_ifab,
                domain_ilo,domain_ihi,domain_jlo,domain_jhi,domain_klo,domain_khi,bx]
               AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
               {
                   bool ed_ilo = (domain_ilo >= bx.smallEnd(0) and domain_ilo <= bx.bigEnd(0)) and
                                 (i == domain_ilo) and 
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(ilo_ifab(i-1,j,k,0))));

                   bool ed_ihi = (domain_ihi >= bx.smallEnd(0) and domain_ihi <= bx.bigEnd(0)) and
                                 (i == domain_ihi) and 
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(ihi_ifab(i+1,j,k,0))));

                   bool ed_jlo = (domain_jlo >= bx.smallEnd(1) and domain_jlo <= bx.bigEnd(1)) and
                                 (j == domain_jlo) and 
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(jlo_ifab(i,j-1,k,0))));

                   bool ed_jhi = (domain_jhi >= bx.smallEnd(1) and domain_jhi <= bx.bigEnd(1)) and
                                 (j == domain_jhi) and 
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(jhi_ifab(i,j+1,k,0))));

                   bool ed_klo = (domain_klo >= bx.smallEnd(2) and domain_klo <= bx.bigEnd(2)) and
                                 (k == domain_klo) and 
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(klo_ifab(i,j,k-1,0))));

                   bool ed_khi = (domain_khi >= bx.smallEnd(2) and domain_khi <= bx.bigEnd(2)) and
                                 (k == domain_khi) and 
                                 (aux::any_of(&bct_Dirichlet[0], &bct_Dirichlet[bct_size],
                                    aux::is_equal<int>(khi_ifab(i,j,k+1,0))));

                   if (flag_fab(i,j,k).isCovered())
                   {
                       xs_fab(i,j,k,slopes_comp+n) = 0.0;
                       ys_fab(i,j,k,slopes_comp+n) = 0.0;
                       zs_fab(i,j,k,slopes_comp+n) = 0.0;
                   }
                   else
                   {
                     // X direction
                     if( flag_fab(i  ,j,k).isSingleValued() or
                        (flag_fab(i-1,j,k).isSingleValued() or not flag_fab(i,j,k).isConnected(-1,0,0)) or
                        (flag_fab(i+1,j,k).isSingleValued() or not flag_fab(i,j,k).isConnected( 1,0,0))) {

                       auto eb_slopes = amrex_calc_slopes_extdir_eb(i,j,k,n,
                                          state_fab,ccent_fab,fcx_fab,fcy_fab,fcz_fab,flag_fab,
                                          ed_ilo,ed_jlo,ed_klo,ed_ihi,ed_jhi,ed_khi,
                                          domain_ilo,domain_jlo,domain_klo,
                                          domain_ihi,domain_jhi,domain_khi);

                       xs_fab(i,j,k,slopes_comp+n) = eb_slopes[0];

                     } else {
                       
                       xs_fab(i,j,k,slopes_comp+n) = amrex_calc_xslope_extdir(i,j,k,n,2,state_fab,
                                                        ed_ilo,ed_ihi,domain_ilo,domain_ihi);
                     }

                     // Y direction
                     if(flag_fab(i,j  ,k).isSingleValued() or
                       (flag_fab(i,j-1,k).isSingleValued() or not flag_fab(i,j,k).isConnected(0,-1,0)) or
                       (flag_fab(i,j+1,k).isSingleValued() or not flag_fab(i,j,k).isConnected(0, 1,0))) {

                       auto eb_slopes = amrex_calc_slopes_extdir_eb(i,j,k,n,
                                          state_fab,ccent_fab,fcx_fab,fcy_fab,fcz_fab,flag_fab,
                                          ed_ilo,ed_jlo,ed_klo,ed_ihi,ed_jhi,ed_khi,
                                          domain_ilo,domain_jlo,domain_klo,
                                          domain_ihi,domain_jhi,domain_khi);

                       ys_fab(i,j,k,slopes_comp+n) = eb_slopes[1];
                     } else {

                       ys_fab(i,j,k,slopes_comp+n) = amrex_calc_yslope_extdir(i,j,k,n,2,state_fab,
                                                        ed_jlo,ed_jhi,domain_jlo,domain_jhi);

                     }

                     // Z direction
                     if(flag_fab(i,j,k  ).isSingleValued() or
                       (flag_fab(i,j,k-1).isSingleValued() or not flag_fab(i,j,k).isConnected(0,0,-1)) or
                       (flag_fab(i,j,k+1).isSingleValued() or not flag_fab(i,j,k).isConnected(0,0, 1))) {

                       auto eb_slopes = amrex_calc_slopes_extdir_eb(i,j,k,n,
                                          state_fab,ccent_fab,fcx_fab,fcy_fab,fcz_fab,flag_fab,
                                          ed_ilo,ed_jlo,ed_klo,ed_ihi,ed_jhi,ed_khi,
                                          domain_ilo,domain_jlo,domain_klo,
                                          domain_ihi,domain_jhi,domain_khi);

                       zs_fab(i,j,k,slopes_comp+n) = eb_slopes[2];

                     } else {

                       zs_fab(i,j,k,slopes_comp+n) = amrex_calc_zslope_extdir(i,j,k,n,2,state_fab,
                                                        ed_klo,ed_khi,domain_klo,domain_khi);

                     }
                   }
               });
           } // end of cut cell region
       } // not covered
    } // MFIter

    xslopes_in[lev]->FillBoundary(geom[lev].periodicity());
    yslopes_in[lev]->FillBoundary(geom[lev].periodicity());
    zslopes_in[lev]->FillBoundary(geom[lev].periodicity());
}
