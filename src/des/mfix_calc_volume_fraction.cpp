#include <mfix.H>
#include <MFIX_FilCC.H>

#include <AMReX_FillPatchUtil.H>
#include <MFIXParticleContainer.H>
#include <MFIX_MFHelpers.H>

#include <MFIX_DEM_Parms.H>
#include <DiffusionOp.H>

void mfix::mfix_calc_volume_fraction (Real& sum_vol)
{
  BL_PROFILE("mfix::mfix_calc_volume_fraction()");

  // Start the timers ...
  const Real strttime = ParallelDescriptor::second();

  if (DEM::solve)
  {
    // This re-calculates the volume fraction within the domain
    // but does not change the values outside the domain

    MultiFab* mf_pointer[nlev];

    if (nlev > 2)
      amrex::Abort("For right now mfix::mfix_calc_volume_fraction can only handle up to 2 levels");

    for (int lev = 0; lev < nlev; lev++) {

      bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) and
                           (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

      if (lev == 0 and OnSameGrids) {

        // If we are already working with the internal mf defined on the
        // particle_box_array, then we just work with this.
        mf_pointer[lev] = m_leveldata[lev]->ep_g;

      } else if (lev == 0 and (not OnSameGrids))  {
        // If ep_g is not defined on the particle_box_array, then we need
        // to make a temporary here and copy into ep_g at the end.
        mf_pointer[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                       pc->ParticleDistributionMap(lev),
                                       m_leveldata[lev]->ep_g->nComp(),
                                       m_leveldata[lev]->ep_g->nGrow());

      } else {
        // If lev > 0 we make a temporary at the coarse resolution
        BoxArray ba_crse(amrex::coarsen(pc->ParticleBoxArray(lev),this->m_gdb->refRatio(0)));
        mf_pointer[lev] = new MultiFab(ba_crse, pc->ParticleDistributionMap(lev),
                                       m_leveldata[lev]->ep_g->nComp(), 1);
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

        flags   = &(particle_ebfactory[lev]->getMultiEBCellFlagFab());
        volfrac = &(particle_ebfactory[lev]->getVolFrac());

      } else {

        Vector<int> ngrow = {1,1,1};
        EBFArrayBoxFactory* crse_factory;

        crse_factory = (makeEBFabFactory(gm, mf_pointer[lev]->boxArray(),
                                        mf_pointer[lev]->DistributionMap(),
                                        ngrow, EBSupport::volume)).release();

        flags   = &(crse_factory->getMultiEBCellFlagFab());
        volfrac = &(crse_factory->getVolFrac());

        delete crse_factory;
      }

      // Deposit particle volume to the grid
      pc->ScalarDeposition(lev, *mf_pointer[lev], volfrac, flags);
    }

    {
      // The deposition occurred on level 0, thus the next few operations
      // only need to be carried out on level 0.
      int lev(0);

      // Move any volume deposited outside the domain back into the domain
      // when BC is either a pressure inlet or mass inflow.
      mfix_deposition_bcs(lev, *mf_pointer[lev]);

      // Sum grid boundaries to capture any material that was deposited into
      // your grid from an adjacent gird.
      mf_pointer[lev]->SumBoundary(gm.periodicity());

      // Fill the boundaries so we calculate the correct average
      // solids volume fraction for periodic boundaries.
      mf_pointer[lev]->FillBoundary(gm.periodicity());

      // Create a copy of the solids volume fraction to use
      // in the redistribution.
      MultiFab* eps_tmp;
      eps_tmp = (MFHelpers::createFrom(*mf_pointer[lev])).release();

      // Clear the grid boundaries to prepare for redistribution which
      // could put material in your ghost cells. Do this AFTER
      // we created the copy so the copy has the correct boundary info.
      mf_pointer[lev]->setBndry(0.0);


      // Move excessive solids volume from small cells to neighboring cells. A copy
      // of the deposition field is made so that when an average is calc
      mfix_redistribute_deposition(lev, *eps_tmp, *mf_pointer[lev], volfrac, flags,
                                   mfix::m_max_solids_volume_fraction);

      // Sum grid boundaries to capture any material that was deposited into
      // your grid from an adjacent gird.
      mf_pointer[lev]->SumBoundary(gm.periodicity());
      mf_pointer[lev]->FillBoundary(gm.periodicity());

      // we no longer need the copy.
      delete eps_tmp;

    }


    int  src_nghost = 1;
    int dest_nghost = 0;
    for (int lev = 1; lev < nlev; lev++)
      mf_pointer[0]->copy(*mf_pointer[lev],0,0, m_leveldata[lev]->ep_g->nComp(),
                          src_nghost, dest_nghost, gm.periodicity(), FabArrayBase::ADD);

    if (nlev > 1)
    {
        // IntVect ref_ratio(this->m_gdb->refRatio(0));

        // Now interpolate from the coarse grid to define the fine grid ep-g
        Interpolater* mapper = &cell_cons_interp;
        int lo_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
        int hi_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
        Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

        BndryFuncArray bfunc(mfix_aux::filcc);

        Real time = 0.0;
        for (int lev = 1; lev < nlev; lev++)
        {
            PhysBCFunct<BndryFuncArray> cphysbc(Geom(lev-1), bcs, bfunc);
            PhysBCFunct<BndryFuncArray> fphysbc(Geom(lev  ), bcs, bfunc);
            m_leveldata[lev]->ep_g->setVal(0.0);
            amrex::InterpFromCoarseLevel(*(m_leveldata[lev]->ep_g), time,
                                         *mf_pointer[lev-1],
                                         0, 0, 1, Geom(lev-1), Geom(lev),
                                         cphysbc, 0, fphysbc, 0,
                                         ref_ratio[0], mapper, bcs, 0);
        }
    }

    // If ep_g is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into ep_g. I believe that we don't
    // need any information in ghost cells so we don't copy those.

    if (mf_pointer[0] != m_leveldata[0]->ep_g)
      m_leveldata[0]->ep_g->copy(*mf_pointer[0], 0, 0, m_leveldata[0]->ep_g->nComp());

    for (int lev = 0; lev < nlev; lev++)
       if (mf_pointer[lev] != m_leveldata[lev]->ep_g)
          delete mf_pointer[lev];

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICDeposition time: " << stoptime << '\n';
    }

    // At this point, we have the particle volume on the fluid grid (ep_s).
    // We will diffuse it first, then convert it to ep_g.
    if(mfix::m_deposition_diffusion_coeff > 0.)
      diffusion_op->diffuse_volfrac(get_ep_g(), mfix::m_deposition_diffusion_coeff);

    for (int lev = 0; lev < nlev; lev++) {
      // Now define this mf = (1 - particle_vol)
      m_leveldata[lev]->ep_g->mult(-1.0, m_leveldata[lev]->ep_g->nGrow());
      m_leveldata[lev]->ep_g->plus( 1.0, m_leveldata[lev]->ep_g->nGrow());

      // We set ep_g to 1 rather than 0 in covered cells so that when we divide by ep_g
      //    following the projection we don't have to protect against divide by 0.
      EB_set_covered(*(m_leveldata[lev]->ep_g),1.0);
    }

    // HACK -- we really should average down (ep_g * volfrac) not ep_g.
    for (int lev = nlev - 1; lev > 0; lev--) {
      amrex::EB_average_down(*(m_leveldata[lev]->ep_g),
                             *(m_leveldata[lev-1]->ep_g),
                             0, 1, m_gdb->refRatio(lev-1));
    }
  } else {
    for (int lev = 0; lev < nlev; lev++)
      m_leveldata[lev]->ep_g->setVal(1.);
  }

  // This sets the values outside walls or periodic boundaries
  for (int lev = 0; lev < nlev; lev++)
    m_leveldata[lev]->ep_g->FillBoundary(geom[lev].periodicity());

  mfix_set_epg_bcs(get_ep_g());

  // Sum up all the values of ep_g[lev], weighted by each cell's EB volfrac
  // Note ep_g = 1 - particle_volume / this_cell_volume where
  //    this_cell_volume = (volfrac * dx * dy * dz)
  // When we define the sum we add up (ep_g * volfrac) so that the total sum
  //    does not depend on whether a particle is in a full or cut cell.
  int lev = 0;
  int comp = 0;

  sum_vol = volWgtSum(lev, *(m_leveldata[lev]->ep_g), comp);
}
