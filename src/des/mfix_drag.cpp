#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_K.H>
#include <MFIX_FilCC.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>
#include <MFIX_MFHelpers.H>
#include <MFIX_DEM_Parms.H>

void
mfix::mfix_calc_drag_fluid (Real time)
{
  const Real strttime = ParallelDescriptor::second();

  mfix_calc_particle_beta(time);

  // ******************************************************************************
  // Now use the beta of individual particles to create the drag terms on the fluid
  // ******************************************************************************
  for (int lev = 0; lev < nlev; lev++)
    m_leveldata[lev]->drag->setVal(0);

  if (nlev > 2)
    amrex::Abort("For right now"
        " MFIXParticleContainer::TrilinearDepositionFluidDragForce can only"
        " handle up to 2 levels");

  Vector< MultiFab* > drag_ptr(nlev, nullptr);

  for (int lev = 0; lev < nlev; lev++) {

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    if (lev == 0 and OnSameGrids) {

      // If we are already working with the internal mf defined on the
      // particle_box_array, then we just work with this.
      drag_ptr[lev] = m_leveldata[lev]->drag;

    } else if (lev == 0 and (not OnSameGrids)) {

      // If beta_mf is not defined on the particle_box_array, then we need
      // to make a temporary here and copy into beta_mf at the end.
      drag_ptr[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                   pc->ParticleDistributionMap(lev),
                                   m_leveldata[lev]->drag->nComp(),
                                   m_leveldata[lev]->drag->nGrow());

    } else {
      // If lev > 0 we make a temporary at the coarse resolution
      BoxArray ba_crse(amrex::coarsen(pc->ParticleBoxArray(lev),this->m_gdb->refRatio(0)));
      drag_ptr[lev] = new MultiFab(ba_crse, pc->ParticleDistributionMap(lev),
                                   m_leveldata[lev]->drag->nComp(), 1);
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (drag_ptr[lev]->nGrow() < 1)
      amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

    drag_ptr[lev]->setVal(0.0, 0, 4, drag_ptr[lev]->nGrow());
  }

  const Geometry& gm = Geom(0);
  const FabArray<EBCellFlagFab>* flags = nullptr;
  const MultiFab* volfrac = nullptr;

  Vector< MultiFab* > tmp_eps(nlev);

  for (int lev = 0; lev < nlev; lev++) {

    tmp_eps[lev] = (MFHelpers::createFrom(*drag_ptr[lev], 0.0)).release();

    // Use level 0 to define the EB factory. If we are not on level 0
    // then create a copy of the coarse factory to use.

    if (lev == 0) {
      flags   = &(particle_ebfactory[lev]->getMultiEBCellFlagFab());
      volfrac = &(particle_ebfactory[lev]->getVolFrac());

    } else {

      Vector<int> ngrow = {1,1,1};
      EBFArrayBoxFactory* crse_factory;

      crse_factory = (makeEBFabFactory(gm, drag_ptr[lev]->boxArray(),
                                      drag_ptr[lev]->DistributionMap(),
                                      ngrow, EBSupport::volume)).release();

      flags   = &(crse_factory->getMultiEBCellFlagFab());
      volfrac = &(crse_factory->getVolFrac());

      delete crse_factory;
    }

    // Deposit the drag force to the grid (beta and beta*particle_vel)
    pc->FluidDragForceDeposition(lev, *tmp_eps[lev], *drag_ptr[lev], volfrac, flags);
  }

  {

    // The deposition occurred on level 0, thus the next few operations
    // only need to be carried out on level 0.
    int lev(0);

    // Move any volume deposited outside the domain back into the domain
    // when BC is either a pressure inlet or mass inflow.
    mfix_deposition_bcs(lev, *drag_ptr[lev]);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    drag_ptr[lev]->SumBoundary(gm.periodicity());
    drag_ptr[lev]->setBndry(0.0);

    // Sum grid boundaries then fill with correct ghost values.
    tmp_eps[lev]->SumBoundary(gm.periodicity());
    tmp_eps[lev]->FillBoundary(gm.periodicity());

    // Move excessive solids volume from small cells to neighboring cells.
    // Note that we don't change tmp_eps but use the redistribution of
    // particle volume to determine how to redistribute the drag forces.
    mfix_redistribute_deposition(lev, *tmp_eps[lev], *drag_ptr[lev], volfrac, flags,
                                 mfix::m_max_solids_volume_fraction);

    // Sum the boundaries again to recapture any solids moved across
    // grid boundaries during the redistribute
    drag_ptr[lev]->SumBoundary(gm.periodicity());
    drag_ptr[lev]->FillBoundary(gm.periodicity());
  }

  // This might not need to exist on all levels. Maybe only level 0.
  for (int lev(0); lev < nlev; ++lev)
    delete tmp_eps[lev];

  int  src_nghost = 1;
  int dest_nghost = 0;
  for (int lev = 1; lev < nlev; lev++) {
    drag_ptr[0]->copy(*drag_ptr[lev],0,0,drag_ptr[0]->nComp(),
                      src_nghost,dest_nghost,gm.periodicity(),FabArrayBase::ADD);
  }

  if (nlev > 1) {

    // IntVect ref_ratio(this->m_gdb->refRatio(0));

    // Now interpolate from the coarse grid to define the fine grid ep-g
    Interpolater* mapper = &cell_cons_interp;
    int lo_bc[3] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
    int hi_bc[3] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
    Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

    BndryFuncArray bfunc(mfix_aux::filcc);

    for (int lev = 1; lev < nlev; lev++) {

      PhysBCFunct<BndryFuncArray> cphysbc(Geom(lev-1), bcs, bfunc);
      PhysBCFunct<BndryFuncArray> fphysbc(Geom(lev  ), bcs, bfunc);
      m_leveldata[lev]->drag->setVal(0);
      amrex::InterpFromCoarseLevel(*m_leveldata[lev]->drag, time, *drag_ptr[lev-1],
                                   0, 0, 1, Geom(lev-1), Geom(lev),
                                   cphysbc, 0, fphysbc, 0,
                                   ref_ratio[0], mapper,
                                   bcs, 0);
    }
  }

  // If mf_to_be_filled is not defined on the particle_box_array, then we need
  // to copy here from drag_ptr into mf_to_be_filled. I believe that we don't
  // need any information in ghost cells so we don't copy those.

  if (drag_ptr[0] != m_leveldata[0]->drag) {
    m_leveldata[0]->drag->copy(*drag_ptr[0], 0, 0, m_leveldata[0]->drag->nComp());
  }

  for (int lev = 0; lev < nlev; lev++) {
    if (drag_ptr[lev] != m_leveldata[lev]->drag)
      delete drag_ptr[lev];
  }

  if (m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "MFIXParticleContainer::TrilinearDepositionFluidDragForce time: " << stoptime << '\n';
  }

  if(mfix::m_deposition_diffusion_coeff > 0.) {
    // Apply mean field diffusion to drag force
    diffusion_op->diffuse_drag(get_drag(), mfix::m_deposition_diffusion_coeff);
  }

  // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
  for (int lev = 0; lev < nlev; lev++)
    m_leveldata[lev]->drag->FillBoundary(geom[lev].periodicity());
}

void
mfix::mfix_calc_drag_particle (Real time)
{
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  BL_PROFILE("mfix::mfix_calc_drag_particle()");

  // Extrapolate velocity Dirichlet bc's to ghost cells
  int extrap_dir_bcs = 1;

  mfix_set_velocity_bcs(time, get_vel_g(), extrap_dir_bcs);

  for (int lev = 0; lev < nlev; lev++)
  {
    Box domain(geom[lev].Domain());
    MultiFab gp_tmp;

    gp_tmp.define(grids[lev],dmap[lev],3,1,MFInfo(),*ebfactory[lev]);

    MultiFab::Copy(gp_tmp, *m_leveldata[lev]->gp, 0, 0, 3, 1);
    gp_tmp.FillBoundary(geom[lev].periodicity());

    //
    // NOTE -- it is essential that we call set_gradp_bcs after calling FillBoundary
    //         because the set_gradp_bcs call hopefully sets the ghost cells exterior
    //         to the domain from ghost cells interior to the domain
    //

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(gp_tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      set_gradp_bcs(bx, lev, gp_tmp[mfi], domain);
    }

    gp_tmp.FillBoundary(geom[lev].periodicity());

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    // Pointer to Multifab for interpolation
    MultiFab* interp_ptr;

    // This is just a sanity check to make sure we're not using covered values
    // We can remove these lines once we're confident in the algorithm
    EB_set_covered(*m_leveldata[0]->vel_g, 0, 3, 1, covered_val);
    EB_set_covered( gp_tmp  , 0, 3, 1, covered_val);

    const int interp_ng = 1;    // Only one layer needed for interpolation
    const int interp_comp = 6;  // Four components (3 vel_g + 3 gp)

    if (OnSameGrids)
    {
      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(grids[lev], dmap[lev], interp_comp, interp_ng, MFInfo(), *ebfactory[lev]);

      // Copy fluid velocity
      interp_ptr->copy(*m_leveldata[lev]->vel_g, 0, 0,
                        m_leveldata[lev]->vel_g->nComp(),
                        m_leveldata[lev]->vel_g->nGrow(),
                        interp_ng);

      // Copy volume fraction
      interp_ptr->copy(gp_tmp,  0, 3,
                       gp_tmp.nComp(),
                       gp_tmp.nGrow(),
                       interp_ng);

      interp_ptr->FillBoundary(geom[lev].periodicity());

    }
    else
    {
      const BoxArray&            pba = pc->ParticleBoxArray(lev);
      const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

      EBFArrayBoxFactory ebfactory_loc(*eb_levels[lev], geom[lev], pba, pdm,
                                       {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                                       EBSupport::full);


      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_comp, interp_ng, MFInfo(), ebfactory_loc);

      // Copy fluid velocity
      interp_ptr->copy(*m_leveldata[lev]->vel_g, 0, 0,
                        m_leveldata[lev]->vel_g->nComp(),
                        m_leveldata[lev]->vel_g->nGrow(),
                        interp_ng);

      // Copy pressure gradient
      interp_ptr->copy(gp_tmp,  0, 3, gp_tmp.nComp(), gp_tmp.nGrow(), interp_ng);

      interp_ptr->FillBoundary(geom[lev].periodicity());
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      const auto dx_array  = geom[lev].CellSizeArray();
      const auto dxi_array = geom[lev].InvCellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const amrex::RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
      const amrex::RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
      const amrex::RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(interp_ptr->Factory());

      const auto cellcent = &(factory.getCentroid());
      const auto bndrycent = &(factory.getBndryCent());
      const auto areafrac = factory.getAreaFrac();

      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        const int np = particles.size();

        Box bx = pti.tilebox ();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  interp_fab = static_cast<EBFArrayBox const&>((*interp_ptr)[pti]);
        const EBCellFlagFab&  flags = interp_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered)
        {
          const auto& interp_array = interp_ptr->array(pti);

          const auto& flags_array = flags.array();

          // We need this until we remove static attribute from mfix::gp0;
          const RealVect gp0_dev(gp0);

          const amrex::Real pmult = DEM::solve ? 1.0 : 0.0;

          if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
          {
            amrex::ParallelFor(np,
              [pstruct,interp_array,interp_comp,gp0_dev,plo,dxi,pmult]
              AMREX_GPU_DEVICE (int pid) noexcept
              {
                // Local array storing interpolated values
                amrex::Real interp_loc[interp_comp];

                MFIXParticleContainer::ParticleType& particle = pstruct[pid];

                trilinear_interp(particle.pos(), &interp_loc[0],
                                 interp_array, plo, dxi, interp_comp);

                Real pbeta = particle.rdata(realData::dragcoeff);

                // Particle drag calculation.  We multiply the particle velocity
                // by "pmult" so that DEM uses the slip velocity. For PIC we
                // only want the fluid velocity as it uses a pseudo implicit
                // slip velocity for parcels.
                particle.rdata(realData::dragx) =
                  pbeta * ( interp_loc[0] - pmult*particle.rdata(realData::velx) ) -
                  (interp_loc[3] + gp0_dev[0]) * particle.rdata(realData::volume);

                particle.rdata(realData::dragy) =
                  pbeta * ( interp_loc[1] - pmult*particle.rdata(realData::vely) ) -
                  (interp_loc[4] + gp0_dev[1]) * particle.rdata(realData::volume);

                particle.rdata(realData::dragz) =
                  pbeta * ( interp_loc[2] - pmult*particle.rdata(realData::velz) ) -
                  (interp_loc[5] + gp0_dev[2]) * particle.rdata(realData::volume);
              });
          }
          else // FAB not all regular
          {

            // Cell centroids
            const auto& ccent_fab = cellcent->array(pti);
            // Centroid of EB
            const auto& bcent_fab = bndrycent->array(pti);
            // Area fractions
            const auto& apx_fab = areafrac[0]->array(pti);
            const auto& apy_fab = areafrac[1]->array(pti);
            const auto& apz_fab = areafrac[2]->array(pti);

            amrex::ParallelFor(np,
              [pstruct,interp_array,interp_comp,flags_array,gp0_dev, pmult,
              plo,dx,dxi,ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab]
              AMREX_GPU_DEVICE (int pid) noexcept
              {
                // Local array storing interpolated values
                amrex::Real interp_loc[interp_comp];

                MFIXParticleContainer::ParticleType& particle = pstruct[pid];
                Real pbeta = particle.rdata(realData::dragcoeff);

                // Cell containing particle centroid
                const int ip = floor((particle.pos(0) - plo[0])*dxi[0]);
                const int jp = floor((particle.pos(1) - plo[1])*dxi[1]);
                const int kp = floor((particle.pos(2) - plo[2])*dxi[2]);

                // The particle is in a covered cell.
                if (flags_array(ip,jp,kp).isCovered())
                {
                  particle.rdata(realData::dragx) = 0.0;
                  particle.rdata(realData::dragy) = 0.0;
                  particle.rdata(realData::dragz) = 0.0;

                // Cut or regular cell and none of the cells in the stencil is covered
                // (Note we can't assume regular cell has no covered cells in the stencil
                //      because of the diagonal case)
                } else {

                  // Upper cell in trilinear stencil
                  const int i = std::floor((particle.pos(0) - plo[0])*dxi[0] + 0.5);
                  const int j = std::floor((particle.pos(1) - plo[1])*dxi[1] + 0.5);
                  const int k = std::floor((particle.pos(2) - plo[2])*dxi[2] + 0.5);

                  // All cells in the stencil are regular. Use
                  // traditional trilinear interpolation
                  if (flags_array(i-1,j-1,k-1).isRegular() and
                      flags_array(i  ,j-1,k-1).isRegular() and
                      flags_array(i-1,j  ,k-1).isRegular() and
                      flags_array(i  ,j  ,k-1).isRegular() and
                      flags_array(i-1,j-1,k  ).isRegular() and
                      flags_array(i  ,j-1,k  ).isRegular() and
                      flags_array(i-1,j  ,k  ).isRegular() and
                      flags_array(i  ,j  ,k  ).isRegular()) {

                    trilinear_interp(particle.pos(), &interp_loc[0],
                                     interp_array, plo, dxi, interp_comp);

                  // At least one of the cells in the stencil is cut or covered
                  } else {

                  const int scomp = 3;
                  fe_interp(particle.pos(), ip, jp, kp, dx, dxi,
                            flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                            interp_array, &interp_loc[0], interp_comp, scomp);


                  } // Cut cell

                  particle.rdata(realData::dragx) =
                    pbeta * ( interp_loc[0] - pmult*particle.rdata(realData::velx) ) -
                    (interp_loc[3] + gp0_dev[0]) * particle.rdata(realData::volume);

                  particle.rdata(realData::dragy) =
                    pbeta * ( interp_loc[1] - pmult*particle.rdata(realData::vely) ) -
                    (interp_loc[4] + gp0_dev[1]) * particle.rdata(realData::volume);

                  particle.rdata(realData::dragz) =
                    pbeta * ( interp_loc[2] - pmult*particle.rdata(realData::velz) ) -
                    (interp_loc[5] + gp0_dev[2]) * particle.rdata(realData::volume);

                } // Not covered
              }); // particle loop
          } // if box not all regular
        } // FAB not covered
      } // pti
    } // omp region

    delete interp_ptr;

  } // lev

  // Reset velocity Dirichlet bc's to face values
  extrap_dir_bcs = 0;
  mfix_set_velocity_bcs(time, get_vel_g(), extrap_dir_bcs);
}
