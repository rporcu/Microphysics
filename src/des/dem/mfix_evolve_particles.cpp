#include <mfix_pc.H>
#include <mfix_des_K.H>
#include <mfix_pc_interactions_K.H>
#include <mfix_pc_updates_K.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_bc.H>
#include <mfix_solvers.H>
#include <mfix_monitors.H>
#include <mfix_calc_cell.H>

using namespace amrex;
using namespace Solvers;

void MFIXParticleContainer::EvolveParticles (int lev,
                                             int /*nstep*/,
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
     * Get EB wall temps/rbox in device array
     ***************************************************************************/
    int bc_tw_count(0);
    for (int bcv(0); bcv < m_boundary_conditions.bc().size(); ++bcv) {
      if (m_boundary_conditions.bc(bcv).type == BCList::eb) {
        bc_tw_count++;
      }
    }

    Gpu::HostVector<RealBox> h_bc_rbv(bc_tw_count);
    Gpu::HostVector<Real>    h_bc_twv(bc_tw_count);
    if (bc_tw_count > 0) {
      int lc0(0);
      for (int bcv(0); bcv < m_boundary_conditions.bc().size(); ++bcv) {
        if (m_boundary_conditions.bc(bcv).type == BCList::eb) {
          h_bc_rbv[lc0] = *(m_boundary_conditions.bc(bcv).region);
          h_bc_twv[lc0] =   m_boundary_conditions.bc(bcv).eb.temperature;
          lc0++;
        }
      }
    }

    Gpu::AsyncArray<RealBox> d_bc_rbv(h_bc_rbv.data(), h_bc_rbv.size());
    Gpu::AsyncArray<Real>    d_bc_twv(h_bc_twv.data(), h_bc_twv.size());

    RealBox* p_bc_rbv = (bc_tw_count > 0) ?  d_bc_rbv.data() : nullptr;
    Real*    p_bc_twv = (bc_tw_count > 0) ?  d_bc_twv.data() : nullptr;

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
    if (ebfactory != NULL) {
      mfix_pc_inflow(lev, 1, 0, dt, solids.solve_enthalpy(), ebfactory);
    }

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
            auto p_realarray = soa.realarray();
            auto p_intarray  = soa.intarray();

            // Number of particles including neighbor particles
            int ntot = nrp;

            // Particle-particle (and particle-wall) forces and torques. We need
            // these to be zero every time we start a new batch (i.e tile and
            // substep) of particles.
            tow[index].clear();
            fc[index].clear();
            cond[index].clear();
            tow[index].resize(3*ntot, 0.0);
            fc[index].resize(3*ntot, 0.0);
            cond[index].resize(ntot, 0.0);

            Real* fc_ptr   = fc[index].dataPtr();
            Real* tow_ptr  = tow[index].dataPtr();
            Real* cond_ptr = cond[index].dataPtr();

            if (solids.update_momentum()) {

              // For debugging: keep track of particle-particle (pfor) and
              // particle-wall (wfor) forces
              pfor[index].clear();
              wfor[index].clear();
              pfor[index].resize(3*ntot, 0.0);
              wfor[index].resize(3*ntot, 0.0);

              // BL_PROFILE_VAR("calc_particle_collisions()", calc_particle_collisions);

              auto& geom = this->Geom(lev);
              const auto dxi = geom.InvCellSizeArray();
              const auto plo = geom.ProbLoArray();
              const auto& phiarr = ls_phi->array(pti);

              const int walls_in_tile = tile_has_walls[index];

              auto nbor_data = m_neighbor_list[lev][index].data();

              const auto& solids_parms = solids.parameters();
              const int solve_enthalpy = solids.solve_enthalpy();

              // now we loop over the neighbor list and compute the forces
              amrex::ParallelFor(nlp,
                  [nrp,pstruct,p_realarray,p_intarray,fc_ptr,tow_ptr,cond_ptr,
                   nbor_data,subdt,ntot,walls_in_tile,ls_refinement,phiarr,plo,
                   dxi,solids_parms,solve_enthalpy,bc_tw_count,p_bc_rbv,p_bc_twv,
                   mew=m_dem.mew(),mew_w=m_dem.mew_w(),kn=m_dem.kn(),
                   kn_w=m_dem.kn_w(),etan=m_dem.etan(),etan_w=m_dem.etan_w(),
                   k_g=m_dem.k_g_dem()]
                AMREX_GPU_DEVICE (int i) noexcept
                {
                    auto particle = pstruct[i];

                    RealVect pos1(particle.pos());

                    RealVect total_force(0.);
                    RealVect total_tow_force(0.);

                    const int istate(p_intarray[SoAintData::state][i]);

                    //**********************************************************
                    // Particle-wall collisions
                    //**********************************************************
                    if (walls_in_tile && (i < nrp)) {

                      particle_walls(p_realarray, p_intarray, i, solve_enthalpy,
                          ls_refinement, phiarr, plo, dxi, subdt, pos1,
                          solids_parms, bc_tw_count, p_bc_rbv, p_bc_twv, k_g,
                          kn_w, etan_w, mew_w, total_force, total_tow_force,
                          cond_ptr, istate);

                    } // tile has walls

                    //**********************************************************
                    // Particle-particle collisions
                    //**********************************************************
                    int has_collisions(0);

                    const auto neighbs = nbor_data.getNeighbors(i);

                    particle_particles(particle, neighbs, p_realarray,
                        p_intarray, i, solve_enthalpy, subdt, pos1, nrp,
                        solids_parms, k_g, kn, etan, mew, total_force,
                        total_tow_force, fc_ptr, tow_ptr, cond_ptr,
                        has_collisions, ntot, istate);

              });

              Gpu::Device::synchronize();

              // BL_PROFILE_VAR_STOP(calc_particle_collisions);

              // BL_PROFILE_VAR("des::solve_particle_velocity_and_position()", des_time_march);

            }

            /********************************************************************
             * Move particles based on collision forces and torques             *
             *******************************************************************/

            const auto p_lo = Geom(lev).ProbLoArray();
            const auto p_hi = Geom(lev).ProbHiArray();

            GpuArray<int,6> lo_hi_bc{m_boundary_conditions.domain_bc(0),
                                     m_boundary_conditions.domain_bc(1),
                                     m_boundary_conditions.domain_bc(2),
                                     m_boundary_conditions.domain_bc(3),
                                     m_boundary_conditions.domain_bc(4),
                                     m_boundary_conditions.domain_bc(5)};

            //Access to added variables
            auto ptile_data = ptile.getParticleTileData();

            const int nspecies_s = solids.nspecies();

            const int idx_X_sn = m_runtimeRealData.X_sn;
            const int idx_mass_txfr = m_runtimeRealData.mass_txfr;
            const int idx_vel_txfr = m_runtimeRealData.vel_txfr;
            const int idx_h_txfr = m_runtimeRealData.h_txfr;

            const int update_mass = solids.update_mass() && solids.solve_species() && reactions.solve();
            const int update_momentum = solids.update_momentum();
            const int solve_enthalpy = solids.solve_enthalpy();
            const int solve_reactions = reactions.solve();

            const Real enthalpy_source = solids.enthalpy_source();

            const int solid_is_a_mixture = solids.isMixture();

            const auto& solids_parms = solids.parameters();

            amrex::ParallelFor(nrp, [pstruct,p_realarray,p_intarray,subdt,
                ptile_data,nspecies_s,idx_X_sn,idx_mass_txfr,idx_vel_txfr,
                idx_h_txfr,update_mass,fc_ptr,cond_ptr,ntot,gravity,tow_ptr,
                p_hi,p_lo,lo_hi_bc,enthalpy_source,update_momentum,time,
                solid_is_a_mixture,solids_parms,solve_enthalpy,solve_reactions,
                abstol,reltol,maxiter]
              AMREX_GPU_DEVICE (int i) noexcept
            {
              ParticleType& p = pstruct[i];

              // Get current particle's mass
              const Real p_mass_old = p_realarray[SoArealData::mass][i];
              Real p_mass_new(p_mass_old);

              // Get current particle's oneOverI
              const Real p_oneOverI_old = p_realarray[SoArealData::oneOverI][i];
              Real p_oneOverI_new(p_oneOverI_old);

              Real p_enthalpy_old(0);
              if (solve_enthalpy) {
                const Real Tp = p_realarray[SoArealData::temperature][i];
                if (solid_is_a_mixture) {
                  for (int n_s(0); n_s < nspecies_s; ++n_s) {
                    p_enthalpy_old += ptile_data.m_runtime_rdata[idx_X_sn+n_s][i]*solids_parms.calc_h_sn<run_on>(Tp,n_s);
                  }
                } else {
                  p_enthalpy_old = solids_parms.calc_h_s<run_on>(Tp);
                }
              }

              // Flag to stop computing particle's quantities if mass_new < 0,
              // i.e. the particle disappears because of chemical reactions
              int proceed = 1;
              Real coeff = 1.;

              //***************************************************************
              // First step: update particles' mass and density
              //***************************************************************
              if (update_mass) {

                part_mass_update(p, ptile_data, p_realarray, i, idx_X_sn,
                    nspecies_s, subdt, p_mass_old, p_mass_new, &p_oneOverI_old,
                    &p_oneOverI_new, proceed, coeff, idx_mass_txfr, 1);

              }

              if (proceed) {
                //***************************************************************
                // Second step: update particles' positions and velocities
                //***************************************************************
                if (update_momentum) {

                  part_momentum_update(p, ptile_data, p_realarray, p_intarray,
                      i, subdt, coeff, p_mass_new, p_oneOverI_new, ntot, fc_ptr,
                      tow_ptr, gravity, p_lo, p_hi, solve_reactions, idx_vel_txfr,
                      lo_hi_bc);

                }

                //***************************************************************
                // Third step: update particles' temperature
                //***************************************************************
                if (solve_enthalpy) {

                  part_enthalpy_update(ptile_data, p_realarray, i, idx_X_sn,
                      nspecies_s, solid_is_a_mixture, solids_parms, subdt, coeff,
                      p_enthalpy_old, p_mass_new, cond_ptr, enthalpy_source,
                      solve_reactions, idx_h_txfr, abstol, reltol, maxiter, 1);

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
