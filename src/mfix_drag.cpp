#include <mfix_F.H>
#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_drag_K.H>
#include "mfix_util_F.H"

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>

void
mfix::mfix_calc_drag_fluid(Real time)
{

  const Real strttime = ParallelDescriptor::second();

  mfix_calc_particle_beta(time);

  // ******************************************************************************
  // Now use the beta of individual particles to create the drag terms on the fluid
  // ******************************************************************************
  for (int lev = 0; lev < nlev; lev++)

    drag[lev] ->setVal(0.0L);

  int fortran_beta_comp = 15;
  int fortran_vel_comp  =  9;


  if (nlev > 2)
    amrex::Abort("For right now MFIXParticleContainer::TrilinearDepositionFluidDragForce can only handle up to 2 levels");

  MultiFab*  drag_ptr[nlev];

  for (int lev = 0; lev < nlev; lev++) {

    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    if (lev == 0 && OnSameGrids) {

      // If we are already working with the internal mf defined on the
      // particle_box_array, then we just work with this.
      drag_ptr[lev] = drag[lev].get();

    } else if (lev == 0 && !OnSameGrids) {

      // If beta_mf is not defined on the particle_box_array, then we need
      // to make a temporary here and copy into beta_mf at the end.
      drag_ptr[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                   pc->ParticleDistributionMap(lev),
                                   drag[lev]->nComp(),
                                   drag[lev]->nGrow());

    } else {
      // If lev > 0 we make a temporary at the coarse resolution
      BoxArray ba_crse(amrex::coarsen(pc->ParticleBoxArray(lev),this->m_gdb->refRatio(0)));
      drag_ptr[lev] = new MultiFab(ba_crse, pc->ParticleDistributionMap(lev),drag[lev]->nComp(),1);
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (drag_ptr[lev]->nGrow() < 1)
      amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

    drag_ptr[lev]->setVal(0.0,0,4,drag_ptr[lev]->nGrow());
  }

  const Geometry& gm  = Geom(0);
  const FabArray<EBCellFlagFab>* flags;
  const MultiFab* volfrac;

  for (int lev = 0; lev < nlev; lev++) {

    // Use level 0 to define the EB factory. If we are not on level 0
    // then create a copy of the coarse factory to use.

    if (lev == 0) {
      flags   = &(ebfactory[lev]->getMultiEBCellFlagFab());
      volfrac = &(ebfactory[lev]->getVolFrac());

    } else {

      Vector<int> ngrow = {1,1,1};
      std::unique_ptr<EBFArrayBoxFactory> crse_factory;

      crse_factory = makeEBFabFactory(gm, drag_ptr[lev]->boxArray(),
                                      drag_ptr[lev]->DistributionMap(),
                                      ngrow, EBSupport::volume);

      flags   = &(crse_factory->getMultiEBCellFlagFab());
      volfrac = &(crse_factory->getVolFrac());
    }


    if (m_deposition_scheme == DepositionScheme::trilinear) {

      pc -> TrilinearDepositionFluidDragForce(lev, *drag_ptr[lev], volfrac, flags,
                                              fortran_beta_comp, fortran_vel_comp);

    } else {

      amrex::Abort("Don't know this deposition_scheme!");

    }

  }


  int  src_nghost = 1;
  int dest_nghost = 0;
  for (int lev = 1; lev < nlev; lev++) {
    drag_ptr[0]->copy(*drag_ptr[lev],0,0,drag_ptr[0]->nComp(),src_nghost,dest_nghost,gm.periodicity(),FabArrayBase::ADD);
  }

  drag_ptr[0]->SumBoundary(gm.periodicity());

  if (nlev > 1) {

    IntVect ref_ratio(this->m_gdb->refRatio(0));

    // Now interpolate from the coarse grid to define the fine grid ep-g
    Interpolater* mapper = &cell_cons_interp;
    int lo_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
    int hi_bc[] = {BCType::foextrap, BCType::foextrap, BCType::foextrap};
    Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

    BndryFuncArray bfunc(phifill);

    Real time = 0.0;
    for (int lev = 1; lev < nlev; lev++) {

      PhysBCFunct<BndryFuncArray> cphysbc(Geom(lev-1), bcs, bfunc);
      PhysBCFunct<BndryFuncArray> fphysbc(Geom(lev  ), bcs, bfunc);
      drag[lev]->setVal(0.0);
      amrex::InterpFromCoarseLevel(*drag[lev], time, *drag_ptr[lev-1],
                                   0, 0, 1, Geom(lev-1), Geom(lev),
                                   cphysbc, 0, fphysbc, 0,
                                   ref_ratio, mapper,
                                   bcs, 0);
    }
  }

  // If mf_to_be_filled is not defined on the particle_box_array, then we need
  // to copy here from drag_ptr into mf_to_be_filled. I believe that we don't
  // need any information in ghost cells so we don't copy those.

  if (drag_ptr[0] != drag[0].get()) {
    drag[0]->copy(*drag_ptr[0],0,0,drag[0]->nComp());
  }

  for (int lev = 0; lev < nlev; lev++) {
    if (drag_ptr[lev] != drag[lev].get())
      delete drag_ptr[lev];
  }

  if (m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "MFIXParticleContainer::TrilinearDepositionFluidDragForce time: " << stoptime << '\n';
  }

#if 0
  mfix_diffuse_drag(drag);
#endif

  // Impose periodic bc's at domain boundaries and fine-fine copies in the interior
  for (int lev = 0; lev < nlev; lev++)
    drag[lev] -> FillBoundary(geom[lev].periodicity());
}






void
mfix::mfix_calc_drag_particle(Real time)
{
  BL_PROFILE("mfix::mfix_calc_drag_particle()");

  // Extrapolate velocity Dirichlet bc's to ghost cells
  int extrap_dir_bcs = 1;
  mfix_set_velocity_bcs(time, vel_g, extrap_dir_bcs);

  for (int lev = 0; lev < nlev; lev++)
  {
    Box      domain(geom[lev].Domain());
    MultiFab gp_tmp;

    gp_tmp.define(grids[lev],dmap[lev],3,1,MFInfo(),*ebfactory[lev]);

    MultiFab::Copy(gp_tmp, *gp[lev], 0, 0, 3, 1);
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

    // Pointers to Multifabs of interest
    MultiFab*  gp_ptr;
    MultiFab* vel_ptr;

    // Temporararies for dual grid case
    std::unique_ptr<MultiFab>   gp_pba;
    std::unique_ptr<MultiFab>  vel_pba;

    // This is just a sanity check to make sure we're not using covered values
    // We can remove these lines once we're confident in the algoirthm
    EB_set_covered(*vel_g[0], 0, 3, 1, covered_val);
    EB_set_covered( gp_tmp  , 0, 3, 1, covered_val);

    if (OnSameGrids)
    {
      gp_ptr  = &gp_tmp;
      vel_ptr = vel_g[lev].get();
    }
    else
    {
      BoxArray            pba = pc->ParticleBoxArray(lev);
      DistributionMapping pdm = pc->ParticleDistributionMap(lev);

      int ng = gp_tmp.nGrow();
      gp_pba.reset(new MultiFab(pba,pdm,gp_tmp.nComp(),ng));
      gp_pba->copy(gp_tmp,0,0,gp_tmp.nComp(),ng,ng);
      gp_pba->FillBoundary(geom[lev].periodicity());

      EBFArrayBoxFactory ebfactory_loc( * eb_levels[lev], geom[lev], pba, pdm,
                                        {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                         m_eb_full_grow_cells}, EBSupport::basic);

      ng = vel_g[lev]->nGrow();
      vel_pba.reset(new MultiFab(pba,pdm,vel_g[lev]->nComp(),ng,MFInfo(), ebfactory_loc));
      vel_pba->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),ng,ng);
      vel_pba->FillBoundary(geom[lev].periodicity());

      gp_ptr  = gp_pba.get();
      vel_ptr = vel_pba.get();
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      const auto dxi = geom[lev].InvCellSizeArray();
      //const auto dx  = geom[lev].CellSizeArray(); // SET_BUT_NOT_USED
      const auto plo = geom[lev].ProbLoArray();

      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        const int np = particles.size();

        Box bx = pti.tilebox ();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_ptr)[pti]);
        const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered)
        {
          const auto& vel_array = vel_ptr->array(pti);
          const auto&  gp_array =  gp_ptr->array(pti);

          const auto& flags_array = flags.array();

          const amrex::Real gp0_x = gp0[0];
          const amrex::Real gp0_y = gp0[1];
          const amrex::Real gp0_z = gp0[2];

          if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
          {
            AMREX_FOR_1D (np, ip,
            {
              amrex::Real velfp[3];
              amrex::Real gradp[3];

              MFIXParticleContainer::ParticleType& particle = pstruct[ip];

              trilinear_interp(particle, &velfp[0], vel_array, plo, dxi);
              trilinear_interp(particle, &gradp[0],  gp_array, plo, dxi);

              Real pbeta = particle.rdata(realData::dragx);

              // Particle drag calculation
              particle.rdata(realData::dragx) =
                pbeta * ( velfp[0] - particle.rdata(realData::velx) ) -
                (gradp[0] + gp0_x) * particle.rdata(realData::volume);

              particle.rdata(realData::dragy) =
                pbeta * ( velfp[1] - particle.rdata(realData::vely) ) -
                (gradp[1] + gp0_y) * particle.rdata(realData::volume);

              particle.rdata(realData::dragz) =
                pbeta * ( velfp[2] - particle.rdata(realData::velz) ) -
                (gradp[2] + gp0_z) * particle.rdata(realData::volume);
            });

            Gpu::synchronize();
          }
          else // FAB not all regular
          {
            // phi is always on the particles grid -- we do not use the refined level here
            const MultiFab & phi = * level_sets[lev];
            const auto& phi_array = phi.array(pti);

            AMREX_FOR_1D (np, ip,
            {
              amrex::Real velfp[3];
              amrex::Real gradp[3];

              MFIXParticleContainer::ParticleType& particle = pstruct[ip];
              Real pbeta = particle.rdata(realData::dragx);

              // This identifies which cell the particle is in
              int iloc = floor((particle.pos(0) - plo[0])*dxi[0]);
              int jloc = floor((particle.pos(1) - plo[1])*dxi[1]);
              int kloc = floor((particle.pos(2) - plo[2])*dxi[2]);

              // Pick upper cell in the stencil
              Real lx = (particle.pos(0) - plo[0])*dxi[0] + 0.5;
              Real ly = (particle.pos(1) - plo[1])*dxi[1] + 0.5;
              Real lz = (particle.pos(2) - plo[2])*dxi[2] + 0.5;

              int i = std::floor(lx);
              int j = std::floor(ly);
              int k = std::floor(lz);

              // Covered cell
              if (flags_array(iloc,jloc,kloc).isCovered())
              {
                particle.rdata(realData::dragx) = 0.0;
                particle.rdata(realData::dragy) = 0.0;
                particle.rdata(realData::dragz) = 0.0;
              }
              else
              {
                // Cut or regular cell and none of the cells in the stencil is covered
                // (Note we can't assume regular cell has no covered cells in the stencil
                //      because of the diagonal case)
                if (not flags_array(i-1,j-1,k-1).isCovered() and
                    not flags_array(i  ,j-1,k-1).isCovered() and
                    not flags_array(i-1,j  ,k-1).isCovered() and
                    not flags_array(i  ,j  ,k-1).isCovered() and
                    not flags_array(i-1,j-1,k  ).isCovered() and
                    not flags_array(i  ,j-1,k  ).isCovered() and
                    not flags_array(i-1,j  ,k  ).isCovered() and
                    not flags_array(i  ,j  ,k  ).isCovered())
                {
                  trilinear_interp(particle, &velfp[0], vel_array, plo, dxi);
                  trilinear_interp(particle, &gradp[0],  gp_array, plo, dxi);
                }
                // At least one of the cells in the stencil is covered
                else
                {
                  Real anrm[3];

                  // Compute distance of the particle from the wall.
                  // (This is the same function we call when computing the
                  // particle-wall collisions)
                  int ls_refinement = 1;
                  //Real dist = interp_level_set(particle, ls_refinement,
                  //phi_array, plo, dxi); // UNUSED_VARIABLE

                  // Compute the normal to the wall in this cell -- it doesn't
                  // matter whether we compute it "at the particle location" or
                  // "at the centroid location" because it interpolates from the
                  // same values of phi.
                  level_set_normal(particle, ls_refinement, &anrm[0], phi_array, plo, dxi);

                  // Particle position must be in [-.5:.5] is relative to cell
                  // center and scaled by dx
                  Real gx = particle.pos(0)*dxi[0] - (iloc + 0.5);
                  Real gy = particle.pos(1)*dxi[1] - (jloc + 0.5);
                  Real gz = particle.pos(2)*dxi[2] - (kloc + 0.5);

                  int ii(0); int jj(0); int kk(0);

                  if (anrm[0] < 0)
                  {
                    ii = iloc - 1;
                  }
                  else
                  {
                    ii = iloc + 1;
                    gx = -gx;
                  }

                  if (anrm[1] < 0)
                  {
                    jj = jloc - 1;
                  }
                  else
                  {
                    jj = jloc + 1;
                    gy = -gy;
                  }

                  if (anrm[2] < 0)
                  {
                    kk = kloc - 1;
                  }
                  else
                  {
                    kk = kloc + 1;
                    gz = -gz;
                  }

                  Real gxy = gx*gy;
                  Real gxz = gx*gz;
                  Real gyz = gy*gz;
                  Real gxyz = gx*gy*gz;

                  for (int n = 0; n < 3; n++)
                  {
                    velfp[n] =
                        (1.0+gx+gy+gz+gxy+gxz+gyz+gxyz) * vel_array(iloc,jloc,kloc,n)
                      + (-gz - gxz - gyz - gxyz)        * vel_array(iloc,jloc,kk  ,n)
                      + (-gy - gxy - gyz - gxyz)        * vel_array(iloc,jj  ,kloc,n)
                      + (gyz + gxyz)                    * vel_array(iloc,jj  ,kk  ,n)
                      + (-gx - gxy - gxz - gxyz)        * vel_array(ii  ,jloc,kloc,n)
                      + (gxz + gxyz)                    * vel_array(ii  ,jloc,kk  ,n)
                      + (gxy + gxyz)                    * vel_array(ii  ,jj  ,kloc,n)
                      + (-gxyz)                         * vel_array(ii  ,jj  ,kk  ,n);
                    gradp[n] =
                        (1.0+gx+gy+gz+gxy+gxz+gyz+gxyz) *  gp_array(iloc,jloc,kloc,n)
                      + (-gz - gxz - gyz - gxyz)        *  gp_array(iloc,jloc,kk  ,n)
                      + (-gy - gxy - gyz - gxyz)        *  gp_array(iloc,jj  ,kloc,n)
                      + (gyz + gxyz)                    *  gp_array(iloc,jj  ,kk  ,n)
                      + (-gx - gxy - gxz - gxyz)        *  gp_array(ii  ,jloc,kloc,n)
                      + (gxz + gxyz)                    *  gp_array(ii  ,jloc,kk  ,n)
                      + (gxy + gxyz)                    *  gp_array(ii  ,jj  ,kloc,n)
                      + (-gxyz)                         *  gp_array(ii  ,jj  ,kk  ,n);

                    // Keep the interpolated velocity between the cell value
                    // and the wall value (0)
                    if ( (velfp[n] > 0.0 and velfp[n] > vel_array(iloc,jloc,kloc,n)) or
                         (velfp[n] < 0.0 and velfp[n] < vel_array(iloc,jloc,kloc,n)) )
                      velfp[n] = vel_array(iloc,jloc,kloc,n);

                    if ( (gradp[n] > 0.0 and gradp[n] > gp_array(iloc,jloc,kloc,n)) or
                         (gradp[n] < 0.0 and gradp[n] < gp_array(iloc,jloc,kloc,n)) )
                      gradp[n] = gp_array(iloc,jloc,kloc,n);
                  }
                } // Cut cell

                particle.rdata(realData::dragx) =
                  pbeta * ( velfp[0] - particle.rdata(realData::velx) ) -
                  (gradp[0] + gp0_x) * particle.rdata(realData::volume);

                particle.rdata(realData::dragy) =
                  pbeta * ( velfp[1] - particle.rdata(realData::vely) ) -
                  (gradp[1] + gp0_y) * particle.rdata(realData::volume);

                particle.rdata(realData::dragz) =
                  pbeta * ( velfp[2] - particle.rdata(realData::velz) ) -
                  (gradp[2] + gp0_z) * particle.rdata(realData::volume);

              } // Not covered
            }); // particle loop

            Gpu::synchronize();
          } // if box not all regular
        } // FAB not covered
      } // pti
    } // omp region
  } // lev

  // Reset velocity Dirichlet bc's to face values
  extrap_dir_bcs = 0;
  mfix_set_velocity_bcs(time, vel_g, extrap_dir_bcs);

}
