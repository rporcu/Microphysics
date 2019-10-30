#include <mfix_F.H>
#include <mfix_eb_F.H>
#include <mfix.H>
#include "mfix_util_F.H"
#include <param_mod_F.H>
#include <bc_mod_F.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>


void mfix::mfix_calc_volume_fraction(Real& sum_vol)
{
    BL_PROFILE("mfix::mfix_calc_volume_fraction()");

    // Start the timers ...
    const Real strttime = ParallelDescriptor::second();


    if (solve_dem)
    {
       // This re-calculates the volume fraction within the domain
       // but does not change the values outside the domain

      int fortran_volume_comp = 5;

      MultiFab* mf_pointer[nlev];
      // Vector< std::unique_ptr<amrex::MultiFab>> mf_pointer(nlev);

      int ncomp = 1;

      if (nlev > 2)
        amrex::Abort("For right now mfix::mfix_calc_volume_fraction can only handle up to 2 levels");

      for (int lev = 0; lev < nlev; lev++) {

        bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                             (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

        if (lev == 0 && OnSameGrids) {

          // If we are already working with the internal mf defined on the
          // particle_box_array, then we just work with this.
          mf_pointer[lev] = ep_g[lev].get();

        } else if (lev == 0 && !OnSameGrids)  {
          // If ep_g is not defined on the particle_box_array, then we need
          // to make a temporary here and copy into ep_g at the end.
          mf_pointer[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                         pc->ParticleDistributionMap(lev),
                                         ncomp, ep_g[lev]->nGrow());

        } else {
          // If lev > 0 we make a temporary at the coarse resolution
          BoxArray ba_crse(amrex::coarsen(pc->ParticleBoxArray(lev),this->m_gdb->refRatio(0)));
          mf_pointer[lev] = new MultiFab(ba_crse, pc->ParticleDistributionMap(lev),
                                         ncomp, 1);
        }

        // We must have ghost cells for each FAB so that a particle in one grid can spread
        // its effect to an adjacent grid by first putting the value into ghost cells of its
        // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
        // to another grid's valid region.
        if (mf_pointer[lev]->nGrow() < 1)
          amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

        mf_pointer[lev]->setVal(0.0,0,1,mf_pointer[lev]->nGrow());
      }

      const Geometry& gm  = Geom(0);
      const FabArray<EBCellFlagFab>* flags;
      const MultiFab* volfrac;

      for (int lev = 0; lev < nlev; lev++) {

        // Use level 0 to define the EB factory
        if (lev == 0) {
          flags   = &(ebfactory[lev]->getMultiEBCellFlagFab());
          volfrac = &(ebfactory[lev]->getVolFrac());

        } else {

          Vector<int> ngrow = {1,1,1};
          std::unique_ptr<EBFArrayBoxFactory> crse_factory;

          crse_factory = makeEBFabFactory(gm, mf_pointer[lev]->boxArray(),
                                          mf_pointer[lev]->DistributionMap(),
                                          ngrow, EBSupport::volume);

          flags   = &(crse_factory->getMultiEBCellFlagFab());
          volfrac = &(crse_factory->getVolFrac());
        }

        // This call deposits the particle volume onto the grid in a PIC-like manner
        pc->TrilinearDepositionScalar(lev, *mf_pointer[lev], volfrac, flags, fortran_volume_comp);


        // Move any field deposited outside the domain back into the domain
        // when BC is pressure inlet and mass inflow.
        Box domain(Geom(lev).Domain());

        const int minf = bc_list.get_minf();
        const int pinf = bc_list.get_pinf();

        for (MFIter mfi(*mf_pointer[lev]); mfi.isValid(); ++mfi) {

          const Box& sbx = (*mf_pointer[lev])[mfi].box();

          // Extract the lower and upper boundaries of Box sbx and Domain
          const IntVect sbx_lo(sbx.loVect()), sbx_hi(sbx.hiVect());
          const amrex::Dim3& dom_lo = amrex::lbound(domain);
          const amrex::Dim3& dom_hi = amrex::ubound(domain);

          // Create a 2D Box collapsing sbx on x-direction
          IntVect sbx_yz_hi(sbx.hiVect());
          sbx_yz_hi[0] = sbx_lo[0];
          const Box sbx_yz(sbx_lo, sbx_yz_hi);

          // Create a 2D Box collapsing sbx on y-direction
          IntVect sbx_xz_hi(sbx.hiVect());
          sbx_xz_hi[1] = sbx_lo[1];
          const Box sbx_xz(sbx_lo, sbx_xz_hi);

          // Create a 2D Box collapsing sbx on z-direction
          IntVect sbx_xy_hi(sbx.hiVect());
          sbx_xy_hi[2] = sbx_lo[2];
          const Box sbx_xy(sbx_lo, sbx_xy_hi);

          Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
          Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();
          Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
          Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();
          Array4<int> const& bc_klo_type = bc_klo[lev]->array();
          Array4<int> const& bc_khi_type = bc_khi[lev]->array();

          Array4<Real> const& vol = mf_pointer[lev]->array(mfi);

          if(sbx_lo[0] < dom_lo.x)
            {
              const int ilo = dom_lo.x;
              AMREX_FOR_3D(sbx_yz, i, j, k,
             {
               if(bc_ilo_type(dom_lo.x-1,j,k,0) == pinf or
                  bc_ilo_type(dom_lo.x-1,j,k,0) == minf)
               {
                 vol(ilo,j,k) += vol(ilo-1,j,k);
                 vol(ilo-1,j,k) = 0;
               }
             });
           }

           if(sbx_hi[0] > dom_hi.x)
           {
             const int ihi = dom_hi.x;
             AMREX_FOR_3D(sbx_yz, i, j, k,
             {
               if(bc_ihi_type(dom_hi.x+1,j,k,0) == pinf or
                  bc_ihi_type(dom_hi.x+1,j,k,0) == minf)
               {
                 vol(ihi,j,k) += vol(ihi+1,j,k);
                 vol(ihi+1,j,k) = 0;
               }
             });
           }

           if(sbx_lo[1] < dom_lo.y)
           {
             const int jlo = dom_lo.y;
             AMREX_FOR_3D(sbx_xz, i, j, k,
             {
               if(bc_jlo_type(i,dom_lo.y-1,k,0) == pinf or
                  bc_jlo_type(i,dom_lo.y-1,k,0) == minf)
               {
                 vol(i,jlo,k) += vol(i,jlo-1,k);
                 vol(i,jlo-1,k) = 0;
               }
             });
           }

           if(sbx_hi[1] > dom_hi.y)
           {
             const int jhi = dom_hi.y;
             AMREX_FOR_3D(sbx_xz, i, j, k,
             {
               if(bc_jhi_type(i,dom_hi.y+1,k,0) == pinf or
                  bc_jhi_type(i,dom_hi.y+1,k,0) == minf)
               {
                 vol(i,jhi,k) += vol(i,jhi+1,k);
                 vol(i,jhi+1,k) = 0;
               }
             });
           }

           if(sbx_lo[2] < dom_lo.z)
           {
             const int klo = dom_lo.z;
             AMREX_FOR_3D(sbx_xy, i, j, k,
             {
               if(bc_klo_type(i,j,dom_lo.z-1,0) == pinf or
                  bc_klo_type(i,j,dom_lo.z-1,0) == minf)
               {
                 vol(i,j,klo) += vol(i,j,klo-1);
                 vol(i,j,klo-1) = 0;
               }
             });
           }

           if(sbx_hi[2] > dom_hi.z)
           {
             const int khi = dom_hi.z;
             AMREX_FOR_3D(sbx_xy, i, j, k,
             {
               if(bc_khi_type(i,j,dom_hi.z+1,0) == pinf or
                  bc_khi_type(i,j,dom_hi.z+1,0) == minf)
               {
                 vol(i,j,khi) += vol(i,j,khi+1);
                 vol(i,j,khi+1) = 0;
               }
             });
           }

           // NOTE: here we do not need host-device synchronization since it is
           // already included in the MFIter destructor
         }
       }



    int  src_nghost = 1;
    int dest_nghost = 0;
    for (int lev = 1; lev < nlev; lev++)
        mf_pointer[0]->copy(*mf_pointer[lev],0,0,ncomp,src_nghost,dest_nghost,gm.periodicity(),FabArrayBase::ADD);

    mf_pointer[0]->SumBoundary(gm.periodicity());

    if (nlev > 1)
    {
        IntVect ref_ratio(this->m_gdb->refRatio(0));

        // Now interpolate from the coarse grid to define the fine grid ep-g
        Interpolater* mapper = &cell_cons_interp;
        int lo_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
        int hi_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
        Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

        BndryFuncArray bfunc(phifill);

        Real time = 0.0;
        for (int lev = 1; lev < nlev; lev++)
        {
            PhysBCFunct<BndryFuncArray> cphysbc(Geom(lev-1), bcs, bfunc);
            PhysBCFunct<BndryFuncArray> fphysbc(Geom(lev  ), bcs, bfunc);
            ep_g[lev]->setVal(0.0);
            amrex::InterpFromCoarseLevel(*ep_g[lev], time, *mf_pointer[lev-1],
                                         0, 0, 1, Geom(lev-1), Geom(lev),
                                         cphysbc, 0, fphysbc, 0,
                                         ref_ratio, mapper,
                                         bcs, 0);
        }
    }

    // If ep_g is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into ep_g. I believe that we don't
    // need any information in ghost cells so we don't copy those.

    if (mf_pointer[0] != ep_g[0].get())
       ep_g[0]->copy(*mf_pointer[0],0,0,ncomp);

    for (int lev = 0; lev < nlev; lev++)
       if (mf_pointer[lev] != ep_g[lev].get())
          delete mf_pointer[lev];

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICDeposition time: " << stoptime << '\n';
    }
      /// remainder goes here






       // At this point, we have the particle volume on the fluid grid (ep_s).
       // We will diffuse it first, then convert it to ep_g.
       mfix_diffuse_eps (ep_g);

       // for (MFIter mfi(ep_g); mfi.isValid(); ++mfi) {
       //   const Box& sbx = ep_g[mfi].box(); // UNUSED_VARIABLE

       //   const Real max_pack = 0.42;
       //   const Real max_pack = 0.21; // UNUSED_VARIABLE

       //   Array4<Real> const& ep_g = ep_g.array(mfi); // UNUSED_VARIABLE

       //   // These lines are commented because this code represents a functionality which
       //   // may be added in future developments
       //   AMREX_FOR_3D(sbx, i, j, k, {
       //       ep_g(i,j,k) = std::max(max_pack, ep_g(i,j,k));
       //     });
       // }


       for (int lev = 0; lev < nlev; lev++)
         {
           // Now define this mf = (1 - particle_vol)
           ep_g[lev]->mult(-1.0,ep_g[lev]->nGrow());
           ep_g[lev]->plus( 1.0,ep_g[lev]->nGrow());

           // We set ep_g to 1 rather than 0 in covered cells so that when we divide by ep_g
           //    following the projection we don't have to protect against divide by 0.
           EB_set_covered(*ep_g[lev],1.0);

         }

       // HACK -- we really should average down (ep_g * volfrac) not ep_g.
       for (int lev = nlev - 1; lev > 0; lev --)
         {
           amrex::EB_average_down(* ep_g[lev], * ep_g[lev - 1],
                                  0, 1, m_gdb->refRatio(lev - 1));
         }
    }
    else
    {
       for (int lev = 0; lev < nlev; lev++)
          ep_g[lev]->setVal(1.);
    }

    // This sets the values outside walls or periodic boundaries
    for (int lev = 0; lev < nlev; lev++)
        ep_g[lev]->FillBoundary(geom[lev].periodicity());

    // Sum up all the values of ep_g[lev], weighted by each cell's EB volfrac
    // Note ep_g = 1 - particle_volume / this_cell_volume where
    //    this_cell_volume = (volfrac * dx * dy * dz)
    // When we define the sum we add up (ep_g * volfrac) so that the total sum
    //    does not depend on whether a particle is in a full or cut cell.
    int lev = 0; int comp = 0;

#ifdef AMREX_USE_CUDA
     bool notInLaunchRegionStatus = Gpu::notInLaunchRegion();

     if(notInLaunchRegionStatus)
       Gpu::setLaunchRegion(true);
#endif

    sum_vol = volWgtSum(lev,*ep_g[lev],comp);

#ifdef AMREX_USE_CUDA
    Gpu::setLaunchRegion(notInLaunchRegionStatus);
#endif
}
