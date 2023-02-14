#include <mfix_des_K.H>

#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_bc.H>
#include <mfix_solvers.H>
#include <mfix_monitors.H>
#include <mfix_calc_cell.H>

using namespace amrex;

MFIXParticleContainer::MFIXParticleContainer (AmrCore* amr_core,
                                              MFIXInitialConditions& initial_conditions,
                                              MFIXBoundaryConditions& boundary_conditions,
                                              MFIXSolidsPhase& solids_in,
                                              MFIXDEM& dem,
                                              MFIXPIC& pic,
                                              MFIXFluidPhase& fluid_in,
                                              MFIXReactions& reactions_in)
    : NeighborParticleContainer<0,0,SoArealData::count,SoAintData::count>(amr_core->GetParGDB(), 1)
    , m_runtimedata_idxs(solids_in.nspecies()*solids_in.solve_species(),
                        reactions_in.nreactions()*reactions_in.solve())
    , nlev (amr_core->maxLevel() + 1)
    , m_initial_conditions(initial_conditions)
    , m_boundary_conditions(boundary_conditions)
    , fluid(fluid_in)
    , solids(solids_in)
    , m_dem(dem)
    , m_pic(pic)
    , reactions(reactions_in)
{
    define();
}


MFIXParticleContainer::MFIXParticleContainer (const Geometry& geom,
                                              const DistributionMapping& dmap,
                                              const BoxArray& ba,
                                              const int nlevel,
                                              MFIXInitialConditions& initial_conditions,
                                              MFIXBoundaryConditions& boundary_conditions,
                                              MFIXSolidsPhase& solids_in,
                                              MFIXDEM& dem,
                                              MFIXPIC& pic,
                                              MFIXFluidPhase& fluid_in,
                                              MFIXReactions& reactions_in)
    : NeighborParticleContainer<0, 0, SoArealData::count,SoAintData::count>(geom, dmap, ba, 1)
    , m_runtimedata_idxs(solids_in.nspecies()*solids_in.solve_species(),
                        reactions_in.nreactions()*reactions_in.solve())
    , nlev(nlevel)
    , m_initial_conditions(initial_conditions)
    , m_boundary_conditions(boundary_conditions)
    , fluid(fluid_in)
    , solids(solids_in)
    , m_dem(dem)
    , m_pic(pic)
    , reactions(reactions_in)
{
    define();
}


void
MFIXParticleContainer::define ()
{
    ReadStaticParameters();
    ReadParameters();

    this->SetVerbose(0);

    // turn off certain components for ghost particle communication
    setRealCommComp(0, true);   // posx
    setRealCommComp(1, true);   // posy
    setRealCommComp(2, true);   // posz
    setRealCommComp(3, true);   // radius
    setRealCommComp(4, false);  // volume
    setRealCommComp(5, false);  // mass
    setRealCommComp(6, false);  // density
    setRealCommComp(7, false);  // oneOverI
    setRealCommComp(8, true);   // velx
    setRealCommComp(9, true);   // vely
    setRealCommComp(10, true);  // velz
    setRealCommComp(11, true);  // omegax
    setRealCommComp(12, true);  // omegay
    setRealCommComp(13, true);  // omegaz
    setRealCommComp(14, false); // statwt
    setRealCommComp(15, false); // dragcoeff
    setRealCommComp(16, false); // dragx
    setRealCommComp(17, false); // dragy
    setRealCommComp(18, false); // dragz
    setRealCommComp(19, false); // cp_s
    setRealCommComp(20, true);  // temperature
    setRealCommComp(21, true);  // convection

#if defined(AMREX_DEBUG) || defined(AMREX_USE_ASSERTION)
    setIntCommComp(0, true);  // id
    setIntCommComp(1, true);  // cpu
#else
    setIntCommComp(0, false); // id
    setIntCommComp(1, false); // cpu
#endif
    setIntCommComp(2, true);  // phase
    setIntCommComp(3, true);  // state
#if MFIX_POLYDISPERSE
    setIntCommComp(4, true);  // ptype
#endif

    // Add solids nspecies components
    for (int n(0); n < m_runtimedata_idxs.count; ++n) {
      AddRealComp(true); // Turn on comm for redistribute on ghosting
      setRealCommComp(21+n, false); // turn off for ghosting
    }
}


void
MFIXParticleContainer::AllocData ()
{
    reserveData();
    resizeData();
}


void
MFIXParticleContainer::PrintParticleCounts ()
{
  const int lev = 0;
  amrex::AllPrintToFile("load_balance") << "Particles on each box: \n";
  long local_count = 0;
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      long np = pti.numParticles();
      local_count += np;
      amrex::AllPrintToFile("load_balance") << "Box:" << pti.index() << ", count: " << np << std::endl;
    }
  amrex::AllPrintToFile("load_balance") << "Total for this process: " << local_count << std::endl << std::endl;
}


void
MFIXParticleContainer::ReadStaticParameters ()
{
    static bool initialized = false;

    if (!initialized)
        initialized = true;
}


void
MFIXParticleContainer::ReadParameters ()
{
    ParmParse pp("solids.newton_solver");

    pp.query("absolute_tol", newton_abstol);
    pp.query("relative_tol", newton_reltol);
    pp.query("max_iterations", newton_maxiter);
}


void
MFIXParticleContainer::EvolveParticles (int lev,
                                        int nstep,
                                        Real dt,
                                        Real time,
                                        RealVect& gravity,
                                        EBFArrayBoxFactory* ebfactory,
                                        EBFArrayBoxFactory* particle_ebfactory,
                                        const MultiFab* ls_phi,
                                        const int ls_refinement,
                                        MultiFab* cost,
                                        std::string& knapsack_weight_type,
                                        int& nsubsteps,
                                        int compute_mass_balance)
{
    BL_PROFILE_REGION_START("mfix_dem::EvolveParticles()");
    BL_PROFILE("mfix_dem::EvolveParticles()");

    amrex::Print() << "Evolving particles on level: " << lev
                   << " ... with fluid dt " << dt << std::endl;


    /****************************************************************************
     * Geometry                                                                 *
     ***************************************************************************/

    const Real* dx = Geom(lev).CellSize();

    /****************************************************************************
     * Init substeps                                                            *
     ***************************************************************************/

    Real subdt;
    // des_init_time_loop(&dt, &nsubsteps, &subdt);
    if ( dt >= m_dem.dtsolid() )
    {
       nsubsteps = static_cast<int>(amrex::Math::ceil(dt / static_cast<amrex::Real>(m_dem.dtsolid())));
       subdt     = dt / nsubsteps;
    } else {
       nsubsteps = 1;
       subdt     = dt;
    }

    /****************************************************************************
     * Get particle EB geometric info
     ***************************************************************************/
    const FabArray<EBCellFlagFab>* flags = &(particle_ebfactory->getMultiEBCellFlagFab());

    /****************************************************************************
     * Init temporary storage:                                                  *
     *   -> particle-particle, and particle-wall forces                         *
     *   -> particle-particle, and particle-wall torques                        *
     ***************************************************************************/
    std::map<PairIndex, Gpu::DeviceVector<Real>> tow;
    std::map<PairIndex, Gpu::DeviceVector<Real>> fc, pfor, wfor;

    std::map<PairIndex, Gpu::DeviceVector<Real>> cond;

    std::map<PairIndex, bool> tile_has_walls;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const Box& bx = pti.tilebox();
        PairIndex index(pti.index(), pti.LocalTileIndex());
        tow[index]  = Gpu::DeviceVector<Real>();
        fc[index]   = Gpu::DeviceVector<Real>();
        pfor[index] = Gpu::DeviceVector<Real>();
        wfor[index] = Gpu::DeviceVector<Real>();
        cond[index] = Gpu::DeviceVector<Real>();

        // Determine if this particle tile actually has any walls
        bool has_wall = false;

        if ((particle_ebfactory != NULL)
            && ((*flags)[pti].getType(amrex::grow(bx,1)) == FabType::singlevalued))
        {
            has_wall = true;
        }
        else
        {
            // We need this test for the case of an inflow boundary:
            // inflow does not appear in the EBFactory but
            // the particles see it as a wall

            // Create the nodal refined box based on the current particle tile
            Box refined_box(amrex::convert(amrex::refine(bx,ls_refinement), IntVect{1,1,1}));

            // Set tol to 1/2 dx
            Real tol = amrex::min(dx[0], amrex::min(dx[1], dx[2])) / 2;

#ifdef AMREX_USE_GPU
            Real ls_min_over_box = ((*ls_phi)[pti]).min<RunOn::Gpu>(refined_box,0);
#else
            Real ls_min_over_box = ((*ls_phi)[pti]).min<RunOn::Cpu>(refined_box,0);
#endif

            if (ls_min_over_box < tol) has_wall = true;
        }

        tile_has_walls[index] = has_wall;
    }

    /****************************************************************************
     * Iterate over sub-steps                                                   *
     ***************************************************************************/

    loc_maxvel  = RealVect(0., 0., 0.);  // Tracks max (absolute) velocity
    loc_maxpfor = RealVect(0., 0., 0.);  // Tracks max particle-particle force
    loc_maxwfor = RealVect(0., 0., 0.);  // Tracks max particle-wall force
    int n = 0; // Counts sub-steps


    // Particle inflow
    if (ebfactory != NULL)
      mfix_pc_inflow(lev, dt, time, solids.solve_enthalpy(), ebfactory);


    // sort particles by cell, this can significantly improve the locality
    if (sort_int > 0 && nstep % sort_int == 0) {
      SortParticlesByBin(m_sorting_bin);
      Print() << "   Sort particles at step " << nstep+1 << "\n";
    }

    const int is_IOProc = int(ParallelDescriptor::IOProcessor());

    const Real abstol = newton_abstol;
    const Real reltol = newton_reltol;
    const int maxiter = newton_maxiter;

    while (n < nsubsteps)
    {
        // Redistribute particles ever so often BUT always update the neighbour
        // list (Note that this fills the neighbour list after every
        // redistribute operation)
        if (solids.update_momentum()) {
          if (n % 25 == 0) {
              clearNeighbors();
              Redistribute(0, 0, 0, 0); // Do not remove negatives
              fillNeighbors();
#ifdef AMREX_USE_GPU
              if (reduceGhostParticles) {
                selectActualNeighbors(MFIXCheckPair(m_dem.neighborhood()));
                updateNeighbors(true);
              }
#endif
              if( m_dem.pneig_flag() ) {
#if MFIX_POLYDISPERSE
                  buildNeighborList(MFIXCheckPolyPair(SoAintData::ptype,m_dem.nptypes(),m_dem.pneighdata()),
                                    SoAintData::ptype, m_dem.prefratdata(), m_dem.nptypes(), false);
#else
                  amrex::Abort("MFIX not built with POLYDISPERSE support");
#endif
              } else {
                  buildNeighborList(MFIXCheckPair(m_dem.neighborhood()), false);
              }
          } else {
              updateNeighbors();
          }
        }

        /********************************************************************
         * Particles routines                                               *
         *******************************************************************/
        BL_PROFILE_VAR("particles_computation", particles_computation);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            // Timer used for load-balancing
            Real wt = ParallelDescriptor::second();

            //const Box& bx = pti.tilebox(); // UNUSED_VARIABLE
            PairIndex index(pti.index(), pti.LocalTileIndex());

            auto& plev  = GetParticles(lev);
            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            ParticleType* pstruct = aos().dataPtr();

            const auto ntp = aos.size();
            const int  nrp = GetParticles(lev)[index].numRealParticles();

            // For multi-grid neighbor list search, we must
            // loop over the ghost particles to find all coll pairs.
            const int  nlp = (m_dem.pneig_flag()) ? ntp : nrp;

            auto& soa = pti.GetStructOfArrays();
            GpuArray<ParticleReal*, SoArealData::count> p_realarray = soa.realarray();
            GpuArray<int*, SoAintData::count> p_intarray  = soa.intarray();

            // Particle-particle (and particle-wall) forces and torques. We need
            // these to be zero every time we start a new batch (i.e tile and
            // substep) of particles.
            tow[index].clear();
            fc[index].clear();
            cond[index].clear();
            tow[index].resize(3*nrp, 0.0);
            fc[index].resize(3*nrp, 0.0);
            cond[index].resize(nrp, 0.0);

            Real* fc_ptr   = fc[index].dataPtr();
            Real* tow_ptr  = tow[index].dataPtr();
            Real* cond_ptr = cond[index].dataPtr();

            if (solids.update_momentum()) {

              // For debugging: keep track of particle-particle (pfor) and
              // particle-wall (wfor) forces
              pfor[index].clear();
              wfor[index].clear();
              pfor[index].resize(3*nrp, 0.0);
              wfor[index].resize(3*nrp, 0.0);

              // BL_PROFILE_VAR("calc_particle_collisions()", calc_particle_collisions);

              const Geometry& geom = this->Geom(lev);
              const GpuArray<Real,AMREX_SPACEDIM> dxi = geom.InvCellSizeArray();
              const GpuArray<Real,AMREX_SPACEDIM> plo = geom.ProbLoArray();
              const Array4<const Real>& phiarr = ls_phi->array(pti);

              const int walls_in_tile = tile_has_walls[index];

              auto nbor_data = m_neighbor_list[lev][index].data();

              constexpr Real small_number = 1.0e-15;

              const auto& bc_parms = m_boundary_conditions.parameters();
              const auto& solids_parms = solids.parameters();
              const auto& dem_parms = m_dem.parameters();
              const int solve_enthalpy = solids.solve_enthalpy();

              ParticleInteractions particle_interactions(pstruct, p_realarray,
                  p_intarray, bc_parms, solids_parms, dem_parms);

              // now we loop over the neighbor list and compute the forces
              amrex::ParallelFor(nlp,
                  [nrp,pstruct,p_realarray,p_intarray,fc_ptr,tow_ptr,cond_ptr,nbor_data,
                   subdt,walls_in_tile,ls_refinement,phiarr,plo,dxi,
                   solve_enthalpy,particle_interactions]
                AMREX_GPU_DEVICE (int i) noexcept
              {
                RealVect total_force_i(0.);
                RealVect total_tow_force_i(0.);
                Real Q_dot_i(0.);

                //**********************************************************
                // Particle-wall collisions
                //**********************************************************
                if (walls_in_tile && (i < nrp)) {

                  particle_interactions.particle_wall(i, solve_enthalpy,
                      ls_refinement, phiarr, plo, dxi, small_number,
                      subdt, total_force_i, total_tow_force_i, Q_dot_i);

                } // tile has walls

                //**********************************************************
                // Particle-particle collisions
                //**********************************************************
                int has_collisions(0);

                const auto neighbs = nbor_data.getNeighbors(i);

                for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit)
                {
                  const int j = mit.index();

                  RealVect total_force_j(0.);
                  RealVect total_tow_force_j(0.);
                  Real Q_dot_j(0.);

                  particle_interactions.particle_particle(i, mit, solve_enthalpy,
                      small_number, subdt, nrp, has_collisions,
                      total_force_i, total_tow_force_i, Q_dot_i,
                      total_force_j, total_tow_force_j, Q_dot_j);

                  if(j < nrp) {

                    HostDevice::Atomic::Add(&fc_ptr[j      ], total_force_j[0]);
                    HostDevice::Atomic::Add(&fc_ptr[j+nrp  ], total_force_j[1]);
                    HostDevice::Atomic::Add(&fc_ptr[j+2*nrp], total_force_j[2]);

                    HostDevice::Atomic::Add(&tow_ptr[j      ], total_tow_force_j[0]);
                    HostDevice::Atomic::Add(&tow_ptr[j+nrp  ], total_tow_force_j[1]);
                    HostDevice::Atomic::Add(&tow_ptr[j+2*nrp], total_tow_force_j[2]);

                    HostDevice::Atomic::Add(&cond_ptr[j], Q_dot_j);
                  }

                } // end neighbor loop

                if(i < nrp) {

                  HostDevice::Atomic::Add(&fc_ptr[i      ], total_force_i[0]);
                  HostDevice::Atomic::Add(&fc_ptr[i+nrp  ], total_force_i[1]);
                  HostDevice::Atomic::Add(&fc_ptr[i+2*nrp], total_force_i[2]);

                  HostDevice::Atomic::Add(&tow_ptr[i      ], total_tow_force_i[0]);
                  HostDevice::Atomic::Add(&tow_ptr[i+nrp  ], total_tow_force_i[1]);
                  HostDevice::Atomic::Add(&tow_ptr[i+2*nrp], total_tow_force_i[2]);

                  HostDevice::Atomic::Add(&cond_ptr[i], Q_dot_i);
                }

                if(i < nrp) {
                  if ((p_intarray[SoAintData::state][i] == 10) && (!has_collisions))
                      p_intarray[SoAintData::state][i] = 1;
                }
              });

              Gpu::Device::synchronize();

              // BL_PROFILE_VAR_STOP(calc_particle_collisions);

              // BL_PROFILE_VAR("des::solve_particle_velocity_and_position()", des_time_march);
            }

            /********************************************************************
             * Move particles based on collision forces and torques             *
             *******************************************************************/

            const GpuArray<Real,AMREX_SPACEDIM> p_lo = Geom(lev).ProbLoArray();
            const GpuArray<Real,AMREX_SPACEDIM> p_hi = Geom(lev).ProbHiArray();

            //Access to added variables
            auto ptile_data = ptile.getParticleTileData();

            const int nspecies_s = solids.nspecies();

            const auto runtimedata_idxs = m_runtimedata_idxs;

            const int update_mass = solids.update_mass() && solids.solve_species() && reactions.solve();
            const int update_momentum = solids.update_momentum();
            const int solve_enthalpy = solids.solve_enthalpy();
            const int solve_reactions = reactions.solve();

            const Real enthalpy_source = solids.enthalpy_source();

            const int solid_is_a_mixture = solids.isMixture();

            const auto& bc_parms = m_boundary_conditions.parameters();
            const auto& solids_parms = solids.parameters();

            MFIXSolvers::Newton nonlinear_solver(abstol, reltol, maxiter, is_IOProc);

            ParticleUpdates particle_updates(pstruct, p_realarray, p_intarray,
                ptile_data, bc_parms, solids_parms);

            amrex::ParallelFor(nrp,
               [nrp,p_realarray,subdt,particle_updates,ptile_data,nspecies_s,
                runtimedata_idxs,update_mass,fc_ptr,cond_ptr,gravity,tow_ptr,
                p_hi,p_lo,enthalpy_source,update_momentum,solve_reactions,
                time,solid_is_a_mixture,solve_enthalpy,nonlinear_solver]
              AMREX_GPU_DEVICE (int i) noexcept
            {
              GpuArray<Real, MFIXSpecies::NMAX> X_sn;
              X_sn.fill(0.);

              // Get current particle's species mass fractions
              for (int n_s(0); n_s < nspecies_s; ++n_s) {
                const int idx = runtimedata_idxs.species_mass_fractions;
                X_sn[n_s] = ptile_data.m_runtime_rdata[idx+n_s][i];
              }

              // Get current particle's mass
              const Real p_mass_old = p_realarray[SoArealData::mass][i];

              // Flag to stop computing particle's quantities if mass_new < 0,
              // i.e. the particle disappears because of chemical reactions
              int proceed = 1;

              //***************************************************************
              // First step: update particles' mass and density
              //***************************************************************
              if (update_mass) {

                const int solve_dem = 1;

                particle_updates.update_mass(i, nspecies_s, X_sn, subdt,
                    runtimedata_idxs, proceed, solve_dem);

              }

              if (proceed) {

                const Real p_mass_new = p_realarray[SoArealData::mass][i];
                const Real mass_coeff = update_mass ? p_mass_old/p_mass_new : 1.;

                //***************************************************************
                // Second step: update particles' positions and velocities
                //***************************************************************
                if (update_momentum) {

                  particle_updates.update_momentum(i, solve_reactions, mass_coeff,
                      fc_ptr, tow_ptr, nrp, subdt, runtimedata_idxs, gravity,
                      p_lo, p_hi);

                }

                //***************************************************************
                // Third step: update particles' temperature
                //***************************************************************
                if (solve_enthalpy) {

                  particle_updates.update_enthalpy(i, nspecies_s, X_sn,
                      solve_reactions, solid_is_a_mixture, mass_coeff,
                      cond_ptr, subdt, runtimedata_idxs, enthalpy_source,
                      nonlinear_solver);
                }
              }
            });

            Gpu::synchronize();

            usr2_des(nrp, ptile);

            /********************************************************************
             * Update runtime cost (used in load-balancing)                     *
             *******************************************************************/

            if (cost)
            {
                // Runtime cost is either (weighted by tile box size):
                //   * time spent
                //   * number of particles
                const Box& tbx = pti.tilebox();
                if (knapsack_weight_type == "RunTimeCosts")
                {
                    wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
                }
                else if (knapsack_weight_type == "NumParticles")
                {
                    wt = nrp / tbx.d_numPts();
                }
                (*cost)[pti].plus<RunOn::Device>(wt, tbx);
            }
        }
        BL_PROFILE_VAR_STOP(particles_computation);

        // Update substep count
        n += 1;

    } // end of loop over substeps


    if ( compute_mass_balance ) {
      ComputeMassProduction(lev, dt);
      ComputeMassOutflow(lev);
    }

    // Redistribute particles at the end of all substeps (note that the particle
    // neighbour list needs to be reset when redistributing).
    if (solids.update_momentum()) {
      clearNeighbors();
      Redistribute(0, 0, 0, 1);
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int nrp   = NumberOfParticles(pti);
        void* particles = pti.GetArrayOfStructs().data();

        usr3_des(nrp,particles);
    }

    amrex::Print() << "done. \n";

    BL_PROFILE_REGION_STOP("mfix_dem::EvolveParticles()");
}


void MFIXParticleContainer::
ComputeAverageDensities (const int lev,
                         const Real time,
                         const std::string&  basename,
                         const Vector<Real>& avg_ro_p,
                         const Vector<Real>& avg_region_x_w,
                         const Vector<Real>& avg_region_x_e,
                         const Vector<Real>& avg_region_y_s,
                         const Vector<Real>& avg_region_y_n,
                         const Vector<Real>& avg_region_z_b,
                         const Vector<Real>& avg_region_z_t )
{
  // Count number of calls -- Used to determine when to create file from scratch
  static int ncalls = 0;
  ++ncalls;

  int  nregions = avg_region_x_w.size();

  if(avg_ro_p.size() > 0)
  {
    //
    // Check the regions are defined correctly
    //
    if ( ( avg_region_x_e.size() != nregions ) ||
         ( avg_region_y_s.size() != nregions ) ||
         ( avg_region_y_n.size() != nregions ) ||
         ( avg_region_z_b.size() != nregions ) ||
         ( avg_region_z_t.size() != nregions ) )
    {
      amrex::Print() << "ComputeAverageTemperatures: some regions are not properly"
        " defined: skipping.";
      return;
    }

    std::vector<Real> region_ro_p(nregions, 0.0);
    std::vector<long> region_np(nregions, 0);

    for ( int nr = 0; nr < nregions; ++nr )
    {
      //amrex::Print() << "size of avg_ro_p " << avg_ro_p[nr] << "\n";

      // This region isn't needed for particle data.
      if( avg_ro_p[nr] == 0) continue;

      // Create Real box for this region
      RealBox avg_region({avg_region_x_w[nr], avg_region_y_s[nr], avg_region_z_b[nr]},
                         {avg_region_x_e[nr], avg_region_y_n[nr], avg_region_z_t[nr]});

      // Jump to next iteration if this averaging region is not valid
      if (!avg_region.ok())
      {
        amrex::Print() << "ComputeAverageDensities: region " << nr
                       << " is invalid: skipping\n";
        continue;
      }

      // Reduce sum operation for np, Tp
      ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
      ReduceData<long, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        Box bx = pti.tilebox();
        RealBox tile_region(bx, Geom(lev).CellSize(), Geom(lev).ProbLo());

        if (tile_region.intersects(avg_region))
        {
          const int np         = NumberOfParticles(pti);
          const AoS &particles = pti.GetArrayOfStructs();
          const ParticleType* pstruct = particles().dataPtr();

          auto& soa = pti.GetStructOfArrays();
          auto p_realarray = soa.realarray();

          reduce_op.eval(np, reduce_data, [pstruct,p_realarray,avg_region]
            AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
          {
            const ParticleType p = pstruct[p_id];

            long l_np = static_cast<long>(0);
            Real l_ro_p = 0.0;

            if (avg_region.contains(p.pos()))
            {
              l_np = static_cast<long>(1);
              l_ro_p = p_realarray[SoArealData::density][p_id];
            }

            return {l_np, l_ro_p};
          });
        }
      }

      ReduceTuple host_tuple = reduce_data.value();
      region_np[nr]   = amrex::get<0>(host_tuple);
      region_ro_p[nr] = amrex::get<1>(host_tuple);
    }

    // Compute parallel reductions
    ParallelDescriptor::ReduceRealSum(region_ro_p.data(), nregions);
    ParallelDescriptor::ReduceLongSum(region_np.data(), nregions);

    // Only the IO processor takes care of the output
    if (ParallelDescriptor::IOProcessor())
    {
      for ( int nr = 0; nr < nregions; ++nr )
      {
        // Skip this region.
        if( avg_ro_p[nr] == 0 ) continue;

        //
        // Compute averages (NaN if NP=0 )
        //
        if (region_np[nr] == 0) {
          region_ro_p[nr] = 0.0;
        }
        else {
          region_ro_p[nr] /= region_np[nr];
        }

        //
        // Print to file
        //
        std::ofstream  ofs;
        std::string    fname;

        fname = basename + "_ro_p_" + std::to_string(nr) + ".dat";

        // Open file
        if ( ncalls == 1 )
        {
          // Create output files only the first time this function is called
          // Use ios:trunc to delete previous content
          ofs.open ( fname.c_str(), std::ios::out | std::ios::trunc );
        }
        else
        {
          // If this is not the first time we write to this file
          // we append to it
          ofs.open ( fname.c_str(), std::ios::out | std::ios::app );
        }

        // Check if file is good
        if ( !ofs.good() )
          amrex::FileOpenFailed ( fname );

        // Print header if first access
        if ( ncalls == 1 )
          ofs << "#  Time   NP  ro" << std::endl;

        ofs << time << " "
            << region_np[nr] << " "
            << region_ro_p[nr] << std::endl;

        ofs.close();
      }
    }
  }
}

void MFIXParticleContainer::
ComputeAverageVelocities (const int lev,
                          const Real time,
                          const std::string&  basename,
                          const Vector<Real>& avg_vel_p,
                          const Vector<Real>& avg_region_x_w,
                          const Vector<Real>& avg_region_x_e,
                          const Vector<Real>& avg_region_y_s,
                          const Vector<Real>& avg_region_y_n,
                          const Vector<Real>& avg_region_z_b,
                          const Vector<Real>& avg_region_z_t )
{
  // Count number of calls -- Used to determine when to create file from scratch
  static int ncalls = 0;
  ++ncalls;

  int  nregions = avg_region_x_w.size();

  if(avg_vel_p.size() > 0)
  {
    //
    // Check the regions are defined correctly
    //
    if ( ( avg_region_x_e.size() != nregions ) ||
         ( avg_region_y_s.size() != nregions ) ||
         ( avg_region_y_n.size() != nregions ) ||
         ( avg_region_z_b.size() != nregions ) ||
         ( avg_region_z_t.size() != nregions ) )
    {
      amrex::Print() << "ComputeAverageVelocities: some regions are not properly"
        " defined: skipping.";
      return;
    }

    std::vector<Real> region_velx(nregions, 0.0);
    std::vector<Real> region_vely(nregions, 0.0);
    std::vector<Real> region_velz(nregions, 0.0);
    std::vector<Real> region_k_en(nregions, 0.0);
    std::vector<long> region_np(nregions, 0);

    for ( int nr = 0; nr < nregions; ++nr )
    {
      //amrex::Print() << "size of avg_vel_p " << avg_vel_p[nr] << "\n";

      // This region isn't needed for particle data.
      if( avg_vel_p[nr] == 0) continue;

      // Create Real box for this region
      RealBox avg_region({avg_region_x_w[nr], avg_region_y_s[nr], avg_region_z_b[nr]},
                         {avg_region_x_e[nr], avg_region_y_n[nr], avg_region_z_t[nr]});

      // Jump to next iteration if this averaging region is not valid
      if ( !avg_region.ok() )
      {
        amrex::Print() << "ComputeAverageVelocities: region " << nr
                       << " is invalid: skipping\n";
        continue;
      }

      // Reduce sum operation for np, velx, vely, velz, kinetic energy
      ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
      ReduceData<long, Real, Real, Real, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        Box bx = pti.tilebox();
        RealBox tile_region(bx, Geom(lev).CellSize(), Geom(lev).ProbLo());

        if (tile_region.intersects(avg_region))
        {
          AoS& aos = pti.GetArrayOfStructs();
          ParticleType* pstruct = aos().dataPtr();

          const int np         = NumberOfParticles(pti);

          SoA& soa = pti.GetStructOfArrays();
          auto p_realarray = soa.realarray();

          reduce_op.eval(np, reduce_data, [pstruct,p_realarray,avg_region]
            AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
          {
            const ParticleType p = pstruct[p_id];

            long l_np   = static_cast<long>(0);
            Real l_velx = 0._rt;
            Real l_vely = 0._rt;
            Real l_velz = 0._rt;
            Real l_k_en = 0._rt;

            if (avg_region.contains(p.pos()))
            {
              const Real mass = p_realarray[SoArealData::mass][p_id];

              l_np = static_cast<long>(1);
              l_velx = p_realarray[SoArealData::velx][p_id];
              l_vely = p_realarray[SoArealData::vely][p_id];
              l_velz = p_realarray[SoArealData::velz][p_id];
              l_k_en = 0.5*mass*(l_velx*l_velx + l_vely*l_vely + l_velz*l_velz);
            }

            return {l_np, l_velx, l_vely, l_velz, l_k_en};
          });
        }
      }

      ReduceTuple host_tuple = reduce_data.value();
      region_np[nr]   = amrex::get<0>(host_tuple);
      region_velx[nr] = amrex::get<1>(host_tuple);
      region_vely[nr] = amrex::get<2>(host_tuple);
      region_velz[nr] = amrex::get<3>(host_tuple);
      region_k_en[nr] = amrex::get<4>(host_tuple);
    }

    // Compute parallel reductions
    ParallelDescriptor::ReduceRealSum(region_velx.data(), nregions);
    ParallelDescriptor::ReduceRealSum(region_vely.data(), nregions);
    ParallelDescriptor::ReduceRealSum(region_velz.data(), nregions);
    ParallelDescriptor::ReduceRealSum(region_k_en.data(), nregions);
    ParallelDescriptor::ReduceLongSum(region_np.data(),   nregions);

    // Only the IO processor takes care of the output
    if (ParallelDescriptor::IOProcessor())
    {
      for ( int nr = 0; nr < nregions; ++nr )
      {
        // Skip this region.
        if( avg_vel_p[nr] == 0 ) continue;

        //
        // Compute averages (NaN if NP=0 )
        //
        if (region_np[nr] == 0) {
          region_velx[nr] = 0.0;
          region_vely[nr] = 0.0;
          region_velz[nr] = 0.0;
          region_k_en[nr] = 0.0;
        }
        else {
          region_velx[nr] /= region_np[nr];
          region_vely[nr] /= region_np[nr];
          region_velz[nr] /= region_np[nr];
          region_k_en[nr] /= region_np[nr];
        }

        //
        // Print to file
        //
        std::ofstream  ofs;
        std::string    fname;

        fname = basename + "_vel_p_" + std::to_string(nr) + ".dat";

        // Open file
        if ( ncalls == 1 )
        {
          // Create output files only the first time this function is called
          // Use ios:trunc to delete previous content
          ofs.open ( fname.c_str(), std::ios::out | std::ios::trunc );
        }
        else
        {
          // If this is not the first time we write to this file
          // we append to it
          ofs.open ( fname.c_str(), std::ios::out | std::ios::app );
        }

        // Check if file is good
        if ( !ofs.good() )
          amrex::FileOpenFailed ( fname );

        // Print header if first access
        if ( ncalls == 1 )
          ofs << "#  Time   NP  U  V  W  KE" << std::endl;

        ofs << time << " "
            << region_np[nr] << " "
            << region_velx[nr] << " "
            << region_vely[nr] << " "
            << region_velz[nr] << " "
            << region_k_en[nr] << std::endl;

        ofs.close();
      }
    }
  }
}

void MFIXParticleContainer::
ComputeAverageTemperatures (const int lev,
                            const Real time,
                            const std::string&  basename,
                            const Vector<Real>& avg_T_p,
                            const Vector<Real>& avg_region_x_w,
                            const Vector<Real>& avg_region_x_e,
                            const Vector<Real>& avg_region_y_s,
                            const Vector<Real>& avg_region_y_n,
                            const Vector<Real>& avg_region_z_b,
                            const Vector<Real>& avg_region_z_t )
{
  // Count number of calls -- Used to determine when to create file from scratch
  static int ncalls = 0;
  ++ncalls;

  int  nregions = avg_region_x_w.size();

  if(avg_T_p.size() > 0)
  {
    //
    // Check the regions are defined correctly
    //
    if ( ( avg_region_x_e.size() != nregions ) ||
         ( avg_region_y_s.size() != nregions ) ||
         ( avg_region_y_n.size() != nregions ) ||
         ( avg_region_z_b.size() != nregions ) ||
         ( avg_region_z_t.size() != nregions ) )
    {
      amrex::Print() << "ComputeAverageTemperatures: some regions are not properly"
        " defined: skipping.";
      return;
    }

    std::vector<Real> region_Tp(nregions, 0.0);
    std::vector<long> region_np(nregions, 0);

    for ( int nr = 0; nr < nregions; ++nr )
    {
      //amrex::Print() << "size of avg_T_p " << avg_T_p[nr] << "\n";

      // This region isn't needed for particle data.
      if( avg_T_p[nr] == 0) continue;

      // Create Real box for this region
      RealBox avg_region({avg_region_x_w[nr], avg_region_y_s[nr], avg_region_z_b[nr]},
                         {avg_region_x_e[nr], avg_region_y_n[nr], avg_region_z_t[nr]});

      // Jump to next iteration if this averaging region is not valid
      if (!avg_region.ok())
      {
        amrex::Print() << "ComputeAverageTemperatures: region " << nr
                       << " is invalid: skipping\n";
        continue;
      }

      // Reduce sum operation for np, Tp
      ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
      ReduceData<long, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        Box bx = pti.tilebox();
        RealBox tile_region(bx, Geom(lev).CellSize(), Geom(lev).ProbLo());

        if (tile_region.intersects(avg_region))
        {
          const int np         = NumberOfParticles(pti);
          const AoS &particles = pti.GetArrayOfStructs();
          const ParticleType* pstruct = particles().dataPtr();

          auto& soa = pti.GetStructOfArrays();
          auto p_realarray = soa.realarray();

          reduce_op.eval(np, reduce_data, [pstruct,p_realarray,avg_region]
            AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
          {
            const ParticleType p = pstruct[p_id];

            long l_np = static_cast<long>(0);
            Real l_Tp = 0.0;

            if (avg_region.contains(p.pos()))
            {
              l_np = static_cast<long>(1);
              l_Tp = p_realarray[SoArealData::temperature][p_id];
            }

            return {l_np, l_Tp};
          });
        }
      }

      ReduceTuple host_tuple = reduce_data.value();
      region_np[nr] = amrex::get<0>(host_tuple);
      region_Tp[nr] = amrex::get<1>(host_tuple);
    }

    // Compute parallel reductions
    ParallelDescriptor::ReduceLongSum(region_np.data(), nregions);
    ParallelDescriptor::ReduceRealSum(region_Tp.data(), nregions);

    // Only the IO processor takes care of the output
    if (ParallelDescriptor::IOProcessor())
    {
      for ( int nr = 0; nr < nregions; ++nr )
      {
        // Skip this region.
        if( avg_T_p[nr] == 0 ) continue;

        //
        // Compute averages (NaN if NP=0 )
        //
        if (region_np[nr] == 0) {
          region_Tp[nr] = 0.0;
        }
        else {
          region_Tp[nr] /= region_np[nr];
        }

        //
        // Print to file
        //
        std::ofstream  ofs;
        std::string    fname;

        fname = basename + "_T_p_" + std::to_string(nr) + ".dat";

        // Open file
        if ( ncalls == 1 )
        {
          // Create output files only the first time this function is called
          // Use ios:trunc to delete previous content
          ofs.open ( fname.c_str(), std::ios::out | std::ios::trunc );
        }
        else
        {
          // If this is not the first time we write to this file
          // we append to it
          ofs.open ( fname.c_str(), std::ios::out | std::ios::app );
        }

        // Check if file is good
        if ( !ofs.good() )
          amrex::FileOpenFailed ( fname );

        // Print header if first access
        if ( ncalls == 1 )
          ofs << "#  Time   NP  T" << std::endl;

        ofs << time << " "
            << region_np[nr] << " "
            << region_Tp[nr] << std::endl;

        ofs.close();
      }
    }
  }
}



namespace {
  typedef std::pair<int, int> BidNp;

  struct PairCompare {
    bool inverse = false;

    PairCompare(const bool a_inverse=false) :inverse(a_inverse) {};
    bool operator() (const BidNp& lhs, const BidNp& rhs)
    {
      return inverse ? lhs.second > rhs.second : lhs.second < rhs.second;
    };
  };

  typedef std::priority_queue<BidNp, Vector<BidNp>, PairCompare> BidNpHeap;
}


void MFIXParticleContainer::partitionParticleGrids(int lev,
                                                   const BoxArray& fba,
                                                   const DistributionMapping& fdmap,
                                                   Real overload_toler,
                                                   Real underload_toler)
{
  // fluid grid info
  const Vector<int>& fpmap    = fdmap.ProcessorMap();
  BoxList            fbl      = fba.boxList();
  Vector<Box>&       fbl_vec  = fbl.data();
  Vector<Box>        fbl_vec0 = fbl_vec;

  // If the mapping between particle and fluid grids hasn't been set,
  // set an 1to1 mapping.
  if (m_pboxid_to_fboxid.size() == 0) {
    m_pboxid_to_fboxid.resize(fbl.size());
    std::iota(m_pboxid_to_fboxid.begin(), m_pboxid_to_fboxid.end(), 0);
  }

  // count particles in fluid grid
  Vector<int> pcount_fbox(fba.size(), 0);
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {
    int fboxid = m_pboxid_to_fboxid[pti.index()];
    pcount_fbox[fboxid] += pti.numParticles();
  }
  ParallelDescriptor::ReduceIntSum(pcount_fbox.dataPtr(), pcount_fbox.size());

  // count particles by rank
  Vector<int> pcount_rank(ParallelDescriptor::NProcs(), 0);
  for (auto i=0; i<pcount_fbox.size(); ++i) {
    pcount_rank[fpmap[i]] += pcount_fbox[i];
  }

  // count total # particles
  if (m_total_numparticle <= 0) {
    m_total_numparticle = 0;
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      m_total_numparticle += static_cast<long>(pti.numRealParticles());
    ParallelDescriptor::ReduceLongSum(m_total_numparticle);
  }

  // find the indices of the overload and underload fluid boxes
  Real avg_np     = static_cast<Real>(m_total_numparticle)
                  / ParallelDescriptor::NProcs();
  int  o_toler_np = static_cast<int>(avg_np * overload_toler);
  int  u_toler_np = static_cast<int>(avg_np * underload_toler);
  Vector<int> overload_fboxid, underload_ranks;
  BoxList     overload_fbl;
  // find overload fluid boxes
  for (auto i=0; i < pcount_fbox.size(); ++i) {
    if (pcount_fbox[i] > o_toler_np) {
      overload_fboxid.push_back(i);
      overload_fbl.push_back(fbl_vec[i]);
    }
  }
  // find underload ranks
  for (auto i=0; i<pcount_rank.size(); ++i)
    if (pcount_rank[i] < u_toler_np)  underload_ranks.push_back(i);
  // debug
  Print() << "avg np: "              << avg_np
          << " overload tolerance "  << overload_toler
          << " underload tolerance " << underload_toler << "\n";

  // re-initialize the map from particle box to fluid box
  m_pboxid_to_fboxid.resize(fbl.size());
  std::iota(m_pboxid_to_fboxid.begin(), m_pboxid_to_fboxid.end(), 0);

  // construct queues for the greedy algorithm
  BidNpHeap o_q;
  for (int i=0; i<overload_fboxid.size(); ++i)
    o_q.push(std::make_pair(i, pcount_fbox[overload_fboxid[i]]));
  BidNpHeap u_q(PairCompare(true));
  for (int rank: underload_ranks)
    u_q.push(std::make_pair(rank, pcount_rank[rank]));

  Vector<int> new_ppmap(fpmap);

  // 1D greedy algorithm
  if (!greedy_3d) {
    // count the particles in bins of the overload boxes
    //currently, each bin is a thin layer of the box
    Vector<int> pcount_bin;           // particle count of bins of overload fbox
    Vector<int> poffset_bin;          // offset of bins of overload fbox
    IntVect     binsize(1024);
    binsize.setVal(greedy_dir, 1);    // thickness of the layer
                                      // layer covers the other 2 direction
    countParticle(lev, overload_fbl, binsize, pcount_bin, poffset_bin);
    ParallelDescriptor::ReduceIntSum(pcount_bin.dataPtr(), pcount_bin.size());

    Vector<int> left_nbin;            // bins left for an overload box
    for (int i=0; i<poffset_bin.size()-1; ++i)
      left_nbin.push_back(poffset_bin[i+1] - poffset_bin[i]);

    // use greedy algorithm to setup the mapping between overload boxes
    // and their cutoff chunks
    while (o_q.size() > 0 && u_q.size() > 0) {
      const BidNp& o_pair = o_q.top();
      const BidNp& u_pair = u_q.top();
      if (o_pair.second < o_toler_np || u_pair.second > avg_np)  break;

      int   room    = amrex::max(static_cast<int>(o_toler_np - u_pair.second), 0);
      int   o_boxid = overload_fboxid[o_pair.first];
      //int   u_boxid = underload_fboxid[u_pair.first];

      // find # bins to chop off
      int   chop_np = 0, chop_nbin = 0;
      for (int ibin = poffset_bin[o_pair.first] + left_nbin[o_pair.first] - 1;
          ibin >= poffset_bin[o_pair.first]; --ibin) {
        int chop_np_tmp = chop_np + pcount_bin[ibin];
        if (  (chop_np_tmp < room && o_pair.second - chop_np_tmp > u_toler_np)
            || chop_nbin < greedy_min_grid_size) {
          chop_np = chop_np_tmp;
          ++chop_nbin;
        }
        else {
          break;
        }
      }// end for ibin

      if (chop_np > room || (o_pair.second - chop_np) < u_toler_np)
        if (chop_np > room)
          break;

      // chop the most overload box
      int chop_pos = fbl_vec[o_boxid].length(greedy_dir)
                   - chop_nbin * binsize[greedy_dir]
                   + fbl_vec[o_boxid].smallEnd(greedy_dir);
      fbl_vec.push_back(fbl_vec[o_boxid].chop(greedy_dir, chop_pos));

      // update mapping for the new particle grid
      m_pboxid_to_fboxid.push_back(o_boxid);
      new_ppmap.push_back(u_pair.first);
      left_nbin[o_pair.first] -= chop_nbin;

      // update np and heap
      BidNp new_o_pair(o_pair.first, o_pair.second - chop_np);
      o_q.pop();
      if (new_o_pair.second > o_toler_np)  o_q.push(new_o_pair);

      BidNp new_u_pair(u_pair.first, u_pair.second + chop_np);
      u_q.pop();
      if (new_u_pair.second < u_toler_np)  u_q.push(new_u_pair);
    }// end while
  }
  // 3D greedy algorithms
  else {
    // particle counts based on particle grids
    iMultiFab np_mf_p(ParticleBoxArray(lev), ParticleDistributionMap(lev), 1, 0);
    np_mf_p.setVal(0);
    // particle counts based on fluid grids
    iMultiFab np_mf_f(The_Pinned_Arena());
    np_mf_f.define(fba, fdmap, 1, 0);
    // count the particles and copy across grids
    countParticle(lev, np_mf_p);
    np_mf_f.ParallelCopy(np_mf_p);

    std::unordered_map<int, Vector<int>> np_cutoff;
    std::unordered_map<int, Vector<int>> id_cutoff;

    // greedy pass to see which box send how many particles to which rank
    while (o_q.size() > 0 && u_q.size() > 0) {
      const BidNp& o_pair = o_q.top();
      const BidNp& u_pair = u_q.top();
      if (o_pair.second < o_toler_np || u_pair.second > avg_np)  break;

      int o_boxid = overload_fboxid[o_pair.first];
      int chop_np = min(o_pair.second - avg_np, avg_np - u_pair.second);

      // Append a placeholder for the cutoff, shape will be set later
      fbl_vec.push_back(Box(IntVect{std::numeric_limits<int>::min()},
                            IntVect{std::numeric_limits<int>::max()}));
      new_ppmap.push_back(u_pair.first);
      m_pboxid_to_fboxid.push_back(o_boxid);


      // the rank that owns the overload box stores the info
      if (fpmap[o_boxid] == ParallelDescriptor::MyProc()) {
        id_cutoff[o_boxid].push_back(fbl_vec.size()-1);
        np_cutoff[o_boxid].push_back(chop_np);
      }

      BidNp new_o_pair(o_pair.first, o_pair.second - chop_np);
      o_q.pop();
      if (new_o_pair.second > o_toler_np)  o_q.push(new_o_pair);

      BidNp new_u_pair(u_pair.first, u_pair.second + chop_np);
      u_q.pop();
      if (new_u_pair.second < u_toler_np)  u_q.push(new_u_pair);
    }// end while

    // setup the shape of cutoffs
    for (auto& kv: id_cutoff) {
      // get # particles in current overload fluid box
      const IArrayBox&         np_fab = np_mf_f[kv.first];
      Array4<int const> const& np_arr = np_fab.const_array();
      const int*               np_ptr = np_arr.dataPtr();

      Box&    fbox0  = fbl_vec0[kv.first];
      Box     remain = fbox0;
      IntVect lo     = fbox0.smallEnd();
      IntVect stride   {1, fbox0.length(0), fbox0.length(0) * fbox0.length(1)};

      // chop the cutoff chunk to fit in the corresponding underload ranks
      // the chop-off workload has be computed by the greedy pass before
      for (int icutoff=0; icutoff<kv.second.size(); ++icutoff) {
        //AllPrintToFile("debug") << "current remain " << remain << "\n";
        int np_target    = np_cutoff[kv.first][icutoff];
        int np_target_lo = static_cast<int>(np_target * underload_toler);
        int np_target_hi = static_cast<int>(np_target *  overload_toler);
        int min_chop_dir = -1;

        Vector<int> np_cutoff_bydir (AMREX_SPACEDIM);
        Vector<int> np_surface_bydir(AMREX_SPACEDIM);
        Vector<int> chop_pos_bydir  (AMREX_SPACEDIM);

        for (int chop_dir=0; chop_dir<AMREX_SPACEDIM; ++chop_dir) {
          if (remain.length(chop_dir) < 2*greedy_min_grid_size) {
            np_surface_bydir[chop_dir] = std::numeric_limits<int>::max();
            np_cutoff_bydir [chop_dir] = std::numeric_limits<int>::max();
            chop_pos_bydir  [chop_dir] = -1;
            continue;
          }

          int dir0 = (chop_dir + 1) % AMREX_SPACEDIM;
          int dir1 = (chop_dir + 2) % AMREX_SPACEDIM;

          int chop_lo  = remain.smallEnd(chop_dir) + greedy_min_grid_size;
          int chop_hi  = remain.bigEnd  (chop_dir);
          int chop_pos = chop_hi;

          int np_tmp       = 0;
          int np_diff_prev = np_target + 1;

          // find the chop position with the closest np to np_target
          while (chop_pos >= chop_lo) {
            int np_slice = 0;
            for (int idir0=remain.smallEnd(dir0); idir0<=remain.bigEnd(dir0); ++idir0) {
              for (int idir1=remain.smallEnd(dir1); idir1<=remain.bigEnd(dir1); ++idir1) {
              // the index is based on the offset to the lower bound of the original box
                int iCell = stride[dir0]     * (idir0    - lo[dir0])
                          + stride[dir1]     * (idir1    - lo[dir1])
                          + stride[chop_dir] * (chop_pos - lo[chop_dir]);
                np_slice += *(np_ptr + iCell);
              }
            }

            int np_diff = std::abs(np_target - np_tmp - np_slice);
            if (np_diff > np_diff_prev && chop_hi - chop_pos + 1 >= greedy_min_grid_size) {
              ++chop_pos;
              break;
            }
            else {
              np_diff_prev = np_diff;
              np_tmp      += np_slice;
              --chop_pos;
            }
          }// end while
          chop_pos = chop_pos == chop_lo-1 ? chop_lo : chop_pos;

          np_cutoff_bydir[chop_dir] = np_tmp;
          chop_pos_bydir [chop_dir] = chop_pos;

          // count surface particles
          Box remain_copy(remain);
          Box cutoff     (remain_copy.chop(chop_dir, chop_pos));
          int np_surface = 0;
          // z faces
          for (int j=cutoff.smallEnd(1); j<=cutoff.bigEnd(1); ++j) {
            for (int i=cutoff.smallEnd(0); i<=cutoff.bigEnd(0); ++i) {
              np_surface += np_arr(i, j, cutoff.smallEnd(2), 0);
              np_surface += np_arr(i, j, cutoff.bigEnd  (2), 0);
            }
          }
          // y faces
          for (int k=cutoff.smallEnd(2); k<=cutoff.bigEnd(2); ++k) {
            for (int i=cutoff.smallEnd(0); i<=cutoff.bigEnd(0); ++i) {
              np_surface += np_arr(i, cutoff.smallEnd(1), k, 0);
              np_surface += np_arr(i, cutoff.bigEnd  (1), k, 0);
            }
          }
          // x faces
          for (int k=cutoff.smallEnd(2); k<=cutoff.bigEnd(2); ++k) {
            for (int j=cutoff.smallEnd(1); j<=cutoff.bigEnd(1); ++j) {
              np_surface += np_arr(cutoff.smallEnd(0), j, k, 0);
              np_surface += np_arr(cutoff.bigEnd  (0), j, k, 0);
            }
          }
          np_surface_bydir[chop_dir] = np_surface;
        }// end for chop_dir

        Vector<int> within_toler_dirs;
        for (int chop_dir=0; chop_dir<AMREX_SPACEDIM; ++chop_dir) {
          if (  np_cutoff_bydir[chop_dir] >= np_target_lo
             && np_cutoff_bydir[chop_dir] <= np_target_hi) {
            within_toler_dirs.push_back(chop_dir);
          }
        }

        // if none of the dir fits, choose the one closest to fit
        if (within_toler_dirs.empty()) {
          int min_diff = std::numeric_limits<int>::max();
          for (int chop_dir=0; chop_dir<AMREX_SPACEDIM; ++chop_dir) {
            int diff = std::abs(np_cutoff_bydir[chop_dir] - np_target);
            if (diff < min_diff) {
              min_diff     = diff;
              min_chop_dir = chop_dir;
            }
          }
        }
        // choose the one with least surface particles
        else {
          int min_np_surface = std::numeric_limits<int>::max();
          for (int dir: within_toler_dirs) {
            if (np_surface_bydir[dir] < min_np_surface) {
              min_np_surface = np_surface_bydir[dir];
              min_chop_dir   = dir;
            }
          }
        }

        // chop the cutoff chunk
        fbl_vec[id_cutoff[kv.first][icutoff]] = \
          remain.chop(min_chop_dir, chop_pos_bydir[min_chop_dir]);
      }// end for icutoff

      // the last underload rank gets whatever left
      fbl_vec[kv.first] = remain;
    }// end for kv

    // Only the overloaded ranks have set up some of the cutoff boxes
    // use an allreduce to inform all the ranks
    size_t nbox = fbl_vec.size();
    Vector<int> ubound_buf(3*nbox);
    Vector<int> lbound_buf(3*nbox);
    // pack the buffer
    for (size_t i=0; i<nbox; ++i) {
      ubound_buf[3*i]   = fbl_vec[i].bigEnd(0);
      ubound_buf[3*i+1] = fbl_vec[i].bigEnd(1);
      ubound_buf[3*i+2] = fbl_vec[i].bigEnd(2);
      lbound_buf[3*i]   = fbl_vec[i].smallEnd(0);
      lbound_buf[3*i+1] = fbl_vec[i].smallEnd(1);
      lbound_buf[3*i+2] = fbl_vec[i].smallEnd(2);
    }

    ParallelDescriptor::ReduceIntMin(ubound_buf.dataPtr(), 3*nbox);
    ParallelDescriptor::ReduceIntMax(lbound_buf.dataPtr(), 3*nbox);

    for (size_t i=0; i<3*nbox; ++i) {
      if (  ubound_buf[i] == std::numeric_limits<int>::max()
         || lbound_buf[i] == std::numeric_limits<int>::min()) {
        Print() << "Greedy 3D cannot balance the current workload, keep "
                << "current particle box array and distribution map.\n";
        return;
      }
    }

    // unpack the buffer
    for (size_t i=0; i<nbox; ++i) {
      fbl_vec[i] =
        Box(IntVect{lbound_buf[3*i], lbound_buf[3*i+1], lbound_buf[3*i+2]},
            IntVect{ubound_buf[3*i], ubound_buf[3*i+1], ubound_buf[3*i+2]});
    }
  }// end else

  // ba and dmap to particle container
  SetParticleBoxArray(lev, BoxArray(fbl));
  SetParticleDistributionMap(lev, DistributionMapping(new_ppmap));
}


void MFIXParticleContainer::setParticleFluidGridMap(const Vector<int>& boxmap)
{
  m_pboxid_to_fboxid.clear();
  m_pboxid_to_fboxid.insert(m_pboxid_to_fboxid.end(), boxmap.begin(), boxmap.end());
}


Real MFIXParticleContainer::particleImbalance()
{
  // # paricles on this process
  long local_count = 0;
  for (MFIXParIter pti(*this, 0); pti.isValid(); ++pti)
      local_count += static_cast<long>(pti.numParticles());
  // count total # particles if not counted
  if (m_total_numparticle <= 0) {
    m_total_numparticle = local_count;
    ParallelDescriptor::ReduceLongSum(m_total_numparticle);
  }
  // max # particles per process
  ParallelDescriptor::ReduceLongMax(local_count,
                                   ParallelDescriptor::IOProcessorNumber());

  return ( static_cast<Real>(m_total_numparticle)
         / ParallelDescriptor::NProcs()
         / (static_cast<Real>(local_count) + 1e-10));
}


void MFIXParticleContainer::countParticle(int lev,
                                          const BoxList& bl,
                                          const IntVect& binsize,
                                          Vector<int>&   pcounts,
                                          Vector<int>&   poffsets)
{
  const auto* boxes = bl.data().data();
  const int   nbox  = bl.size();

  pcounts.clear();
  poffsets.resize(1, 0);

  // find the total number of bins and offsets for each box's bins
  int total_nbin = 0;
  for (int i=0; i<nbox; ++i) {
    IntVect boxsize = boxes[i].size();
    AMREX_ASSERT_WITH_MESSAGE(AMREX_D_TERM(
         (boxsize[0] < binsize[0] || boxsize[0] % binsize[0] == 0),
      && (boxsize[1] < binsize[1] || boxsize[1] % binsize[1] == 0),
      && (boxsize[2] < binsize[2] || boxsize[2] % binsize[2] == 0)),
      "ERROR: For now, the greedy balance requires the box size to be less than"
      " or divisible over the bin size");

    total_nbin += AMREX_D_TERM(  (boxsize[0] + binsize[0] - 1) / binsize[0],
                              * ((boxsize[1] + binsize[1] - 1) / binsize[1]),
                              * ((boxsize[2] + binsize[2] - 1) / binsize[2]));
    poffsets.push_back(total_nbin);
  }
  pcounts.resize(total_nbin, 0);

  // particle tiles and geometry of this level
  const auto& geom   = Geom(lev);
  const auto  domain = geom.Domain();
  const auto  dx_inv = geom.InvCellSizeArray();
  const auto  prob_lo = geom.ProbLoArray();

  amrex::Gpu::DeviceVector<Box> d_boxes(nbox);
  amrex::Gpu::DeviceVector<int> d_poffsets(poffsets.size());
  amrex::Gpu::DeviceVector<int> d_pcounts(pcounts.size(), 0);

  Gpu::copy(Gpu::hostToDevice, bl.data().begin(), bl.data().end(), d_boxes.begin());
  Gpu::copy(Gpu::hostToDevice, poffsets.begin(), poffsets.end(), d_poffsets.begin());

  Gpu::synchronize();

  const auto* p_d_boxes = d_boxes.data();
  const auto* p_d_poffsets = d_poffsets.data();

  auto* p_d_pcounts  = d_pcounts.dataPtr();

  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {
    const auto& aos     = pti.GetArrayOfStructs();
    const auto* pstruct = aos().dataPtr();
    int         np      = pti.numParticles();

    ParallelFor(np, [pstruct, prob_lo, dx_inv, domain, nbox, p_d_boxes,
                     binsize, p_d_pcounts, p_d_poffsets]
      AMREX_GPU_DEVICE (int i) noexcept
      {
        IntVect cell_ijk = getParticleCell(pstruct[i], prob_lo, dx_inv, domain);
        Box     box_tmp;
        for (int ibox=0; ibox<nbox; ++ibox) {
          if (p_d_boxes[ibox].contains(cell_ijk)) {
            int ibin = getTileIndex(cell_ijk, p_d_boxes[ibox], true, binsize, box_tmp);
            ibin += p_d_poffsets[ibox];
            Gpu::Atomic::AddNoRet(p_d_pcounts + ibin, 1);
          }
        }
      });// end parallel for
  }// end for pariter

  Gpu::synchronize();
  Gpu::copy(Gpu::deviceToHost, d_pcounts.begin(), d_pcounts.end(), pcounts.begin());
}


void MFIXParticleContainer::countParticle(int lev, iMultiFab& np_mf)
{
  // particle tiles and geometry of this level
  const auto& geom    = Geom(lev);
  const auto  domain  = geom.Domain();
  const auto  dx_inv  = geom.InvCellSizeArray();
  const auto  prob_lo = geom.ProbLoArray();

  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {
    const auto& aos     = pti.GetArrayOfStructs();
    const auto* pstruct = aos().dataPtr();
    int         np      = pti.numParticles();

    Array4<int> const& np_a4 = np_mf.array(pti);

    ParallelFor(np, [pstruct, np_a4, prob_lo, dx_inv, domain]
      AMREX_GPU_DEVICE (int i) noexcept
      {
        IntVect     ijk = getParticleCell(pstruct[i], prob_lo, dx_inv, domain);
        Gpu::Atomic::AddNoRet(&np_a4(ijk[0], ijk[1], ijk[2], 0), 1);
      });// end parallel for
  }// end for pariter

  Gpu::synchronize();
}


void MFIXParticleContainer::verifyParticleCount()
{
  // # paricles on this process
  long local_count = 0;
  for (MFIXParIter pti(*this, 0); pti.isValid(); ++pti)
      local_count += static_cast<long>(pti.numParticles());
  // count total # particles if not counted
  if (m_total_numparticle <= 0) {
    m_total_numparticle = local_count;
    ParallelDescriptor::ReduceLongSum(m_total_numparticle);
    Print() << "total # particles updated to " << m_total_numparticle;
  }
  else {
    long total_numparticle = local_count;
    ParallelDescriptor::ReduceLongSum(total_numparticle);
    if (total_numparticle != m_total_numparticle)
      Print() << "total # particles does not match"
              << " old " << m_total_numparticle
              << " new " << total_numparticle << "\n";
  }
}


void MFIXParticleContainer::printGhostParticleCount()
{
  const int lev = 0;

  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
  {
    PairIndex  index(pti.index(), pti.LocalTileIndex());
    auto  ptile = GetParticles(lev)[index];
    auto& aos   = ptile.GetArrayOfStructs();
    ParticleType* pstruct = aos().dataPtr();
    const int  nrp = ptile.numRealParticles();
    const int  ngp = ptile.numTotalParticles() - nrp;

    // create vector to mark if a ghost is a neighbor of some real particle.
    Gpu::DeviceVector<int> nbrids(ngp, 0);
    auto* pnbrids  = nbrids.dataPtr();
    Print(2) << "nrp " << nrp << " ngp " << ngp << "\n";

    const Box  box = pti.validbox();
    const auto& geom     = Geom(lev);
    const auto  domain   = geom.Domain();
    const auto  dx_inv   = geom.InvCellSizeArray();
    const auto  prob_lo  = geom.ProbLoArray();

    auto  nbr_data = m_neighbor_list[lev][index].data();

    amrex::ParallelFor(nrp, [pstruct, prob_lo, dx_inv, domain, nbr_data,
                            pnbrids, box, nrp]
                      AMREX_GPU_DEVICE (int i) noexcept {
      const auto nbrs = nbr_data.getNeighbors(i);
      for (auto mit = nbrs.begin(); mit != nbrs.end(); ++mit) {
        IntVect cell_ijk = getParticleCell(*mit, prob_lo, dx_inv, domain);
        int     pid      = mit.index();
        if (!box.contains(cell_ijk)) {
          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pid >= nrp, "ghost particle has real index");
          HostDevice::Atomic::Add(pnbrids + (pid - nrp), 1);
        }
      }
    });
    Gpu::Device::synchronize();

    int  nnbr  = 0;
    int* pnnbr = &nnbr;
    amrex::ParallelFor(ngp, [pnnbr, pnbrids] AMREX_GPU_DEVICE (int i) noexcept {
      if (pnbrids[i] > 0)  HostDevice::Atomic::Add(pnnbr, 1);
    });
    Gpu::Device::synchronize();

    AllPrintToFile("nbr") << "pbox " << pti.index() << " "
      << nnbr << " neighbors "
      << ngp << " ghost particles\n";
  }// end for pti
}
