#include <mfix_des_K.H>

#include <mfix_solids_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_bc_parms.H>
#include <mfix_solvers.H>

using namespace amrex;
using namespace Solvers;


int  MFIXParticleContainer::domain_bc[6] {0};


MFIXParticleContainer::MFIXParticleContainer (AmrCore* amr_core,
                                              SolidsPhase& arg_solids,
                                              Reactions& arg_reactions)
    : NeighborParticleContainer<AoSrealData::count,AoSintData::count,
                                SoArealData::count,SoAintData::count>(amr_core->GetParGDB(), 1)
    , m_runtimeRealData(arg_solids.nspecies, arg_reactions.nreactions)
    , solids(arg_solids)
    , reactions(arg_reactions)
{
    ReadStaticParameters();

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
    setRealCommComp(20, false); // temperature
    setRealCommComp(21, false); // convection

#if defined(AMREX_DEBUG) || defined(AMREX_USE_ASSERTION)
    setIntCommComp(0, true); // id
    setIntCommComp(1, true); // cpu
#else
    setIntCommComp(0, false); // id
    setIntCommComp(1, false); // cpu
#endif
    setIntCommComp(2, true);  // phase
    setIntCommComp(3, false); // state

    nlev = amr_core->maxLevel() + 1;

    // Add solids.nspecies components
    for (int n(0); n < m_runtimeRealData.count; ++n) {
      AddRealComp(false); // Turn off ghost particle communication
    }
}

void MFIXParticleContainer::AllocData ()
{
    reserveData();
    resizeData();
}

void MFIXParticleContainer::PrintParticleCounts ()
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

void MFIXParticleContainer::printParticles ()
{
    const int lev = 0;
    auto& plevel = GetParticles(lev);

    for (auto& kv : plevel)
    {
       const auto& particles = kv.second.GetArrayOfStructs();
       auto& soa = kv.second.GetStructOfArrays();
       auto p_realarray = soa.realarray();
       auto p_intarray = soa.intarray();

       for (int i = 0; i < particles.numParticles(); ++i)
       {
          std::cout << "Particle ID  = " << i << " " << std::endl;
          std::cout << "X            = " << particles[i].pos(0) << " " << std::endl;
          std::cout << "Y            = " << particles[i].pos(1) << " " << std::endl;
          std::cout << "Z            = " << particles[i].pos(2) << " " << std::endl;
          std::cout << "state        = " << p_intarray[SoAintData::state][i] << " " << std::endl;
          std::cout << "phase        = " << p_intarray[SoAintData::phase][i] << " " << std::endl;
          std::cout << "Real properties = " << std::endl;

          for (int j = 0; j < SoArealData::count; j++)
            std::cout << "property " << j << "  = " << p_realarray[j][i] << " " << std::endl;

          std::cout << std::endl;
       }
    }
}

void MFIXParticleContainer::ReadStaticParameters ()
{
    static bool initialized = false;

    if (!initialized)
        initialized = true;
}

void MFIXParticleContainer::EvolveParticles (int lev,
                                             int nstep,
                                             Real dt,
                                             Real time,
                                             RealVect& gravity,
                                             EBFArrayBoxFactory* ebfactory,
                                             const MultiFab* ls_phi,
                                             const int ls_refinement,
                                             MultiFab* cost,
                                             std::string& knapsack_weight_type,
                                             int& nsubsteps,
                                             const int advect_enthalpy,
                                             const Real enthalpy_source,
                                             const int update_mass,
                                             const int update_momentum,
                                             const int update_enthalpy)
{
    BL_PROFILE_REGION_START("mfix_dem::EvolveParticles()");
    BL_PROFILE("mfix_dem::EvolveParticles()");

    Real eps = std::numeric_limits<Real>::epsilon();

    amrex::Print() << "Evolving particles on level: " << lev
                   << " ... with fluid dt " << dt << std::endl;

    /****************************************************************************
     * DEBUG flag toggles:                                                      *
     *   -> Print number of collisions                                          *
     *   -> Print max (over substeps) particle velocity at each time step       *
     *   -> Print max particle-wall and particle-particle forces                *
     ***************************************************************************/

    // Debug level controls the detail of debug output:
    //   -> debug_level = 0 : no debug output
    //   -> debug_level = 1 : debug output for every fluid step
    //   -> debug_level = 2 : debug output for every substep
    const int debug_level = 0;

    /****************************************************************************
     * Geometry                                                                 *
     ***************************************************************************/

    const Real* dx = Geom(lev).CellSize();

    /****************************************************************************
     * Init substeps                                                            *
     ***************************************************************************/

    Real subdt;
    // des_init_time_loop(&dt, &nsubsteps, &subdt);
    if ( dt >= DEM::dtsolid )
    {
       nsubsteps = static_cast<int>(amrex::Math::ceil(dt / static_cast<amrex::Real>(DEM::dtsolid)));
       subdt     = dt / nsubsteps;
    } else {
       nsubsteps = 1;
       subdt     = dt;
    }

    /****************************************************************************
     * Get particle EB geometric info
     ***************************************************************************/
    const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());

    /****************************************************************************
     * Init temporary storage:                                                  *
     *   -> particle-particle, and particle-wall forces                         *
     *   -> particle-particle, and particle-wall torques                        *
     ***************************************************************************/
    std::map<PairIndex, Gpu::DeviceVector<Real>> tow;
    std::map<PairIndex, Gpu::DeviceVector<Real>> fc, pfor, wfor;

    std::map<PairIndex, bool> tile_has_walls;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const Box& bx = pti.tilebox();
        PairIndex index(pti.index(), pti.LocalTileIndex());
        tow[index]  = Gpu::DeviceVector<Real>();
        fc[index]   = Gpu::DeviceVector<Real>();
        pfor[index] = Gpu::DeviceVector<Real>();
        wfor[index] = Gpu::DeviceVector<Real>();

        // Determine if this particle tile actually has any walls
        bool has_wall = false;

        if ((ebfactory != NULL)
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

            Real ls_min_over_box = ((*ls_phi)[pti]).min<RunOn::Gpu>(refined_box,0);

            if (ls_min_over_box < tol) has_wall = true;
        }

        tile_has_walls[index] = has_wall;
    }

    /****************************************************************************
     * Iterate over sub-steps                                                   *
     ***************************************************************************/

    int ncoll_total = 0;  // Counts total number of collisions
    loc_maxvel  = RealVect(0., 0., 0.);  // Tracks max (absolute) velocity
    loc_maxpfor = RealVect(0., 0., 0.);  // Tracks max particle-particle force
    loc_maxwfor = RealVect(0., 0., 0.);  // Tracks max particle-wall force
    int n = 0; // Counts sub-steps

    // sort particles by cell, this can significantly improve the locality
    if (sort_int > 0 && nstep % sort_int == 0) {
      SortParticlesByBin(m_sorting_bin);
      Print() << "   Sort particles at step " << nstep+1 << "\n";
    }

    while (n < nsubsteps)
    {
        int ncoll = 0;  // Counts number of collisions (over sub-steps)

        // Redistribute particles ever so often BUT always update the neighbour
        // list (Note that this fills the neighbour list after every
        // redistribute operation)
        if (update_momentum) {
          if (n % 25 == 0) {
              clearNeighbors();
              Redistribute(0, 0, 0, 1);
              fillNeighbors();
              // send in "false" for sort_neighbor_list option

              buildNeighborList(MFIXCheckPair(DEM::neighborhood), false);
          } else {
              updateNeighbors();
          }

          /********************************************************************
           * Compute number of Particle-Particle collisions
           *******************************************************************/
          if (debug_level > 0) 
          {
#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion())
            {
              // Reduce sum operation for ncoll
              ReduceOps<ReduceOpSum> reduce_op;
              ReduceData<int> reduce_data(reduce_op);
              using ReduceTuple = typename decltype(reduce_data)::Type;

              for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
              {
                PairIndex index(pti.index(), pti.LocalTileIndex());

                const int nrp = GetParticles(lev)[index].numRealParticles();

                auto& plev = GetParticles(lev);
                auto& ptile = plev[index];
                auto& aos   = ptile.GetArrayOfStructs();
                ParticleType* pstruct = aos().dataPtr();

                auto& soa = ptile.GetStructOfArrays();
                auto p_realarray = soa.realarray();

                auto nbor_data = m_neighbor_list[lev][index].data();

                constexpr Real small_number = 1.0e-15;

                reduce_op.eval(nrp, reduce_data, [nrp, pstruct,p_realarray,
                    nbor_data,small_number]
                  AMREX_GPU_DEVICE (int i) -> ReduceTuple
                {
                  int l_ncoll(0);

                  ParticleType p1 = pstruct[i];
                  const RealVect pos1 = p1.pos();
                  const Real radius1 = p_realarray[SoArealData::radius][i];

                  const auto neighbs = nbor_data.getNeighbors(i);
                  for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit)
                  {
                    const auto p2 = *mit;
                    const int j = mit.index();

                    const RealVect pos2 = p2.pos();
                    const Real radius2 = p_realarray[SoArealData::radius][j];

                    Real r2 = (pos1 - pos2).radSquared();
                    Real r_lm = radius1 + radius2;

                    if (r2 <= (r_lm-small_number)*(r_lm-small_number))
                      l_ncoll += (j < nrp ? 2 : 1);
                  }

                  return {l_ncoll};
                });
              }

              ReduceTuple host_tuple = reduce_data.value();
              ncoll += amrex::get<0>(host_tuple);
            }
            else
#endif
            {
#ifdef _OPENMP
#pragma omp parallel reduction(+:ncoll) if (Gpu::notInLaunchRegion())
#endif
              for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
              {
                PairIndex index(pti.index(), pti.LocalTileIndex());

                const int nrp = GetParticles(lev)[index].numRealParticles();

                auto& plev = GetParticles(lev);
                auto& ptile = plev[index];
                auto& aos   = ptile.GetArrayOfStructs();
                ParticleType* pstruct = aos().dataPtr();

                auto& soa = ptile.GetStructOfArrays();
                auto p_realarray = soa.realarray();

                auto nbor_data = m_neighbor_list[lev][index].data();

                constexpr Real small_number = 1.0e-15;

                for(int i(0); i < nrp; ++i)
                {
                  ParticleType p1 = pstruct[i];
                  const RealVect pos1 = p1.pos();
                  const Real radius1 = p_realarray[SoArealData::radius][i];

                  const auto neighbs = nbor_data.getNeighbors(i);
                  for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit)
                  {
                    const auto p2 = *mit;
                    const int j = mit.index();

                    const RealVect pos2 = p2.pos();
                    const Real radius2 = p_realarray[SoArealData::radius][j];

                    Real r2 = (pos1 - pos2).radSquared();

                    Real r_lm = radius1 + radius2;

                    if (r2 <= (r_lm-small_number)*(r_lm-small_number))
                    {
                      ncoll += 1;
                    }
                  }
                }
              }
            } 
          } // end if (debug_level > 0)
        }

        /********************************************************************
         * Particles routines                                               *
         *******************************************************************/
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        BL_PROFILE_VAR("particles_computation", particles_computation);
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            // Timer used for load-balancing
            Real wt = ParallelDescriptor::second();

            //const Box& bx = pti.tilebox(); // UNUSED_VARIABLE
            PairIndex index(pti.index(), pti.LocalTileIndex());

            auto& plev = GetParticles(lev);
            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            ParticleType* pstruct = aos().dataPtr();

            const int nrp = GetParticles(lev)[index].numRealParticles();

            auto& soa   = pti.GetStructOfArrays();
            auto p_realarray = soa.realarray();
            auto p_intarray = soa.intarray();

            // Number of particles including neighbor particles
            int ntot = nrp;

              // Particle-particle (and particle-wall) forces and torques. We need
              // these to be zero every time we start a new batch (i.e tile and
              // substep) of particles.
              tow[index].clear();
              fc[index].clear();
              tow[index].resize(3*ntot, 0.0);
              fc[index].resize(3*ntot, 0.0);

              Real* fc_ptr = fc[index].dataPtr();
              Real* tow_ptr = tow[index].dataPtr();

            if (update_momentum) {

              // For debugging: keep track of particle-particle (pfor) and
              // particle-wall (wfor) forces
              pfor[index].clear();
              wfor[index].clear();
              pfor[index].resize(3*ntot, 0.0);
              wfor[index].resize(3*ntot, 0.0);

              /********************************************************************
               * Particle-Wall collision forces (and torques)                     *
               *******************************************************************/

              if (tile_has_walls[index])
              {
                  // Calculate forces and torques from particle-wall collisions
                  BL_PROFILE_VAR("calc_wall_collisions()", calc_wall_collisions);

                  auto& geom = this->Geom(lev);
                  const auto dxi = geom.InvCellSizeArray();
                  const auto plo = geom.ProbLoArray();
                  const auto& phiarr = ls_phi->array(pti);

                  amrex::ParallelFor(nrp,
                    [pstruct,p_realarray,p_intarray,ls_refinement,phiarr,plo,dxi,subdt,ntot,fc_ptr,tow_ptr,
                     local_mew_w=DEM::mew_w,local_kn_w=DEM::kn_w,local_etan_w=DEM::etan_w]
                    AMREX_GPU_DEVICE (int i) noexcept
                    {
                      auto particle = pstruct[i];

                      Real rp = p_realarray[SoArealData::radius][i];

                      RealVect pos(particle.pos());

                      Real ls_value = interp_level_set(pos, ls_refinement, phiarr, plo, dxi);

                      Real overlap_n = rp - ls_value;

                      if (ls_value < rp)
                      {
                          RealVect normal(0.);
                          level_set_normal(pos, ls_refinement, normal, phiarr, plo, dxi);

                          normal[0] *= -1;
                          normal[1] *= -1;
                          normal[2] *= -1;

                          RealVect v_rot(0.);
                          v_rot[0] = ls_value * p_realarray[SoArealData::omegax][i];
                          v_rot[1] = ls_value * p_realarray[SoArealData::omegay][i];
                          v_rot[2] = ls_value * p_realarray[SoArealData::omegaz][i];

                          RealVect vreltrans(0.);
                          RealVect cprod(0.);

                          cross_product(v_rot, normal, cprod);
                          vreltrans[0] = p_realarray[SoArealData::velx][i] + cprod[0];
                          vreltrans[1] = p_realarray[SoArealData::vely][i] + cprod[1];
                          vreltrans[2] = p_realarray[SoArealData::velz][i] + cprod[2];

                          Real vreltrans_norm = dot_product(vreltrans, normal);

                          RealVect vrel_t(0.);
                          vrel_t[0] = vreltrans[0] - vreltrans_norm*normal[0];
                          vrel_t[1] = vreltrans[1] - vreltrans_norm*normal[1];
                          vrel_t[2] = vreltrans[2] - vreltrans_norm*normal[2];

                          int phase = p_intarray[SoAintData::phase][i];

                          Real kn_des_w   = local_kn_w;
                          Real etan_des_w = local_etan_w(phase-1);

                          // NOTE - we don't use the tangential components right now,
                          // but we might in the future
                          // Real kt_des_w = DEM::kt_w;
                          // Real etat_des_w = DEM::etat_w[phase-1];

                          RealVect fn(0.);
                          RealVect ft(0.);
                          RealVect overlap_t(0.);
                          Real mag_overlap_t(0.);

                          // calculate the normal contact force
                          fn[0] = -(kn_des_w*overlap_n*normal[0]
                                  + etan_des_w*vreltrans_norm*normal[0]);
                          fn[1] = -(kn_des_w*overlap_n*normal[1]
                                  + etan_des_w*vreltrans_norm*normal[1]);
                          fn[2] = -(kn_des_w*overlap_n*normal[2]
                                  + etan_des_w*vreltrans_norm*normal[2]);

                          // calculate the tangential displacement
                          overlap_t[0] = subdt*vrel_t[0];
                          overlap_t[1] = subdt*vrel_t[1];
                          overlap_t[2] = subdt*vrel_t[2];

                          mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));

                          if (mag_overlap_t > 0.0) {
                              Real fnmd = local_mew_w * sqrt(dot_product(fn, fn));
                              RealVect tangent(0.);
                              tangent[0] = overlap_t[0]/mag_overlap_t;
                              tangent[1] = overlap_t[1]/mag_overlap_t;
                              tangent[2] = overlap_t[2]/mag_overlap_t;
                              ft[0] = -fnmd * tangent[0];
                              ft[1] = -fnmd * tangent[1];
                              ft[2] = -fnmd * tangent[2];
                          } else {
                              ft[0] = 0.0;
                              ft[1] = 0.0;
                              ft[2] = 0.0;
                          }

#ifdef _OPENMP
#pragma omp critical
                          {
#endif
                            Gpu::Atomic::Add(&fc_ptr[i         ], fn[0] + ft[0]);
                            Gpu::Atomic::Add(&fc_ptr[i + ntot  ], fn[1] + ft[1]);
                            Gpu::Atomic::Add(&fc_ptr[i + 2*ntot], fn[2] + ft[2]);
#ifdef _OPENMP
                          }
#endif

                          RealVect tow_force(0.);

                          cross_product(normal, ft, tow_force);

#ifdef _OPENMP
#pragma omp critical
                          {
#endif
                            Gpu::Atomic::Add(&tow_ptr[i         ], ls_value*tow_force[0]);
                            Gpu::Atomic::Add(&tow_ptr[i + ntot  ], ls_value*tow_force[1]);
                            Gpu::Atomic::Add(&tow_ptr[i + 2*ntot], ls_value*tow_force[2]);
#ifdef _OPENMP
                          }
#endif
                      }
                  });

                  // Debugging: copy data from the fc (all forces) vector to
                  // the wfor (wall forces) vector.
                  if (debug_level > 0) {
                      Gpu::synchronize();
                      for (size_t i = 0; i < wfor[index].size(); i++ ) {
                          wfor[index][i] = fc[index][i];
                      }
                  }

                  BL_PROFILE_VAR_STOP(calc_wall_collisions);
              }

              /********************************************************************
               * Particle-Particle collision forces (and torques)                 *
               *******************************************************************/

              // BL_PROFILE_VAR("calc_particle_collisions()", calc_particle_collisions);

              auto nbor_data = m_neighbor_list[lev][index].data();

              constexpr Real small_number = 1.0e-15;

              // now we loop over the neighbor list and compute the forces
              amrex::ParallelFor(nrp,
                  [nrp,pstruct,p_realarray,p_intarray,fc_ptr,tow_ptr,nbor_data,
#if defined(AMREX_DEBUG) || defined(AMREX_USE_ASSERTION)
                   eps,
#endif
                   subdt,ntot,local_mew=DEM::mew,local_kn=DEM::kn,
                   local_etan=DEM::etan]
                AMREX_GPU_DEVICE (int i) noexcept
                {
                    auto particle = pstruct[i];

                    RealVect pos1(particle.pos());

                    const auto neighbs = nbor_data.getNeighbors(i);
                    for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit)
                    {
                        const auto p2 = *mit;
                        const int j = mit.index();

                        Real dist_x = p2.pos(0) - pos1[0];
                        Real dist_y = p2.pos(1) - pos1[1];
                        Real dist_z = p2.pos(2) - pos1[2];

                        Real r2 = dist_x*dist_x +
                                  dist_y*dist_y +
                                  dist_z*dist_z;

                        const Real p1radius = p_realarray[SoArealData::radius][i];
                        const Real p2radius = p_realarray[SoArealData::radius][j];

                        Real r_lm = p1radius + p2radius;

                        AMREX_ASSERT_WITH_MESSAGE(
                            !(particle.id() == p2.id() &&
                                 particle.cpu() == p2.cpu()),
                          "A particle should not be its own neighbor!");

                        if ( r2 <= (r_lm - small_number)*(r_lm - small_number) )
                        {
                            Real dist_mag = sqrt(r2);

                            AMREX_ASSERT(dist_mag >= eps);

                            Real dist_mag_inv = 1.e0/dist_mag;

                            RealVect normal(0.);
                            normal[0] = dist_x * dist_mag_inv;
                            normal[1] = dist_y * dist_mag_inv;
                            normal[2] = dist_z * dist_mag_inv;

                            Real overlap_n = r_lm - dist_mag;
                            Real vrel_trans_norm;
                            RealVect vrel_t(0.);

                            RealVect p1vel(p_realarray[SoArealData::velx][i],
                                           p_realarray[SoArealData::vely][i],
                                           p_realarray[SoArealData::velz][i]);

                            RealVect p2vel(p_realarray[SoArealData::velx][j],
                                           p_realarray[SoArealData::vely][j],
                                           p_realarray[SoArealData::velz][j]);

                            RealVect p1omega(p_realarray[SoArealData::omegax][i],
                                             p_realarray[SoArealData::omegay][i],
                                             p_realarray[SoArealData::omegaz][i]);

                            RealVect p2omega(p_realarray[SoArealData::omegax][j],
                                             p_realarray[SoArealData::omegay][j],
                                             p_realarray[SoArealData::omegaz][j]);

                            cfrelvel(p1vel, p2vel, p1radius, p2radius, p1omega,
                                p2omega, vrel_trans_norm, vrel_t, normal, dist_mag);

                            int phase1 = p_intarray[SoAintData::phase][i];
                            int phase2 = p_intarray[SoAintData::phase][j];

                            Real kn_des = local_kn;
                            Real etan_des = local_etan(phase1-1,phase2-1);

                            // NOTE - we don't use the tangential components right now,
                            // but we might in the future
                            // Real kt_des = DEM::kt;
                            // Real etat_des = DEM::etat[phase1-1][phase2-1];

                            RealVect fn(0.);
                            RealVect ft(0.);
                            RealVect overlap_t(0.);
                            Real mag_overlap_t(0.);

                            // calculate the normal contact force
                            fn[0] = -(kn_des*overlap_n*normal[0]
                                    + etan_des*vrel_trans_norm*normal[0]);
                            fn[1] = -(kn_des*overlap_n*normal[1]
                                    + etan_des*vrel_trans_norm*normal[1]);
                            fn[2] = -(kn_des*overlap_n*normal[2]
                                    + etan_des*vrel_trans_norm*normal[2]);

                            // calculate the tangential overlap
                            overlap_t[0] = subdt*vrel_t[0];
                            overlap_t[1] = subdt*vrel_t[1];
                            overlap_t[2] = subdt*vrel_t[2];
                            mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));

                            if (mag_overlap_t > 0.0) {
                                Real fnmd = local_mew * sqrt(dot_product(fn, fn));
                                RealVect tangent(0.);
                                tangent[0] = overlap_t[0]/mag_overlap_t;
                                tangent[1] = overlap_t[1]/mag_overlap_t;
                                tangent[2] = overlap_t[2]/mag_overlap_t;
                                ft[0] = -fnmd * tangent[0];
                                ft[1] = -fnmd * tangent[1];
                                ft[2] = -fnmd * tangent[2];
                            } else {
                                ft[0] = 0.0;
                                ft[1] = 0.0;
                                ft[2] = 0.0;
                            }


#ifdef _OPENMP
#pragma omp critical
                            {
#endif
                              Gpu::Atomic::Add(&fc_ptr[i         ], fn[0] + ft[0]);
                              Gpu::Atomic::Add(&fc_ptr[i + ntot  ], fn[1] + ft[1]);
                              Gpu::Atomic::Add(&fc_ptr[i + 2*ntot], fn[2] + ft[2]);

                              if (j < nrp)
                              {
                                Gpu::Atomic::Add(&fc_ptr[j         ], -(fn[0] + ft[0]));
                                Gpu::Atomic::Add(&fc_ptr[j + ntot  ], -(fn[1] + ft[1]));
                                Gpu::Atomic::Add(&fc_ptr[j + 2*ntot], -(fn[2] + ft[2]));
                              }
#ifdef _OPENMP
                            }
#endif

                            Real dist_cl1 = 0.5 * (dist_mag + (p1radius*p1radius - p2radius*p2radius) * dist_mag_inv);
                            dist_cl1 = dist_mag - dist_cl1;

                            Real dist_cl2 = 0.5 * (dist_mag + (p2radius*p2radius - p1radius*p1radius) * dist_mag_inv);
                            dist_cl2 = dist_mag - dist_cl2;

                            RealVect tow_force(0.);

                            cross_product(normal, ft, tow_force);

#ifdef _OPENMP
#pragma omp critical
                            {
#endif
                              Gpu::Atomic::Add(&tow_ptr[i         ], dist_cl1*tow_force[0]);
                              Gpu::Atomic::Add(&tow_ptr[i + ntot  ], dist_cl1*tow_force[1]);
                              Gpu::Atomic::Add(&tow_ptr[i + 2*ntot], dist_cl1*tow_force[2]);

                              if (j < nrp)
                              {
                                  Gpu::Atomic::Add(&tow_ptr[j         ], dist_cl2*tow_force[0]);
                                  Gpu::Atomic::Add(&tow_ptr[j + ntot  ], dist_cl2*tow_force[1]);
                                  Gpu::Atomic::Add(&tow_ptr[j + 2*ntot], dist_cl2*tow_force[2]);
                              }
#ifdef _OPENMP
                            }
#endif
                        }
                    }
              });

              Gpu::Device::synchronize();

              // Debugging: copy data from the fc (all forces) vector to the wfor
              // (wall forces) vector. Note that since fc already contains the
              // wall forces, these need to be subtracted here.
              if (debug_level > 0)
              {
                  for (size_t i = 0; i < pfor[index].size(); i++ ) {
                      pfor[index][i] = fc[index][i] - wfor[index][i];
                  }
              }

              // BL_PROFILE_VAR_STOP(calc_particle_collisions);

              // BL_PROFILE_VAR("des::update_particle_velocity_and_position()", des_time_march);

            }

            /********************************************************************
             * Move particles based on collision forces and torques             *
             *******************************************************************/

            const auto p_lo = Geom(lev).ProbLoArray();
            const auto p_hi = Geom(lev).ProbHiArray();

            int x_lo_bc = BC::domain_bc[0];
            int x_hi_bc = BC::domain_bc[1];
            int y_lo_bc = BC::domain_bc[2];
            int y_hi_bc = BC::domain_bc[3];
            int z_lo_bc = BC::domain_bc[4];
            int z_hi_bc = BC::domain_bc[5];

            //Access to added variables
            auto ptile_data = ptile.getParticleTileData();

            const int nspecies_s = solids.nspecies;

            const int idx_X_sn = m_runtimeRealData.X_sn;
            const int idx_mass_sn_txfr = m_runtimeRealData.mass_sn_txfr;
            const int idx_vel_s_txfr = m_runtimeRealData.vel_s_txfr;
            const int idx_h_s_txfr = m_runtimeRealData.h_s_txfr;

            const int local_update_mass = update_mass && solids.solve_species && reactions.solve;
            const int local_update_enthalpy = update_enthalpy && advect_enthalpy;
            const int local_solve_reactions = reactions.solve;

            const int solid_is_a_mixture = solids.is_a_mixture;

            auto& solids_parms = *solids.parameters;

            amrex::ParallelFor(nrp, [pstruct,p_realarray,p_intarray,subdt,
                ptile_data,nspecies_s,idx_X_sn,idx_mass_sn_txfr,idx_vel_s_txfr,
                idx_h_s_txfr,local_update_mass,fc_ptr,ntot,gravity,tow_ptr,eps,
                p_hi,p_lo,x_lo_bc,x_hi_bc,y_lo_bc,y_hi_bc,z_lo_bc,z_hi_bc,
                enthalpy_source,update_momentum,local_solve_reactions,time,
                solid_is_a_mixture,solids_parms,local_update_enthalpy]
              AMREX_GPU_DEVICE (int i) noexcept
            {
              ParticleType& p = pstruct[i];

              GpuArray<Real,SPECIES::NMAX> X_sn;

              // Get current particle's mass
              const Real p_mass_old = p_realarray[SoArealData::mass][i];
              Real p_mass_new(p_mass_old);

              // Get current particle's density
              const Real p_density_old = p_realarray[SoArealData::density][i];
              Real p_density_new(p_density_old);

              // Get current particle's oneOverI
              const Real p_oneOverI_old = p_realarray[SoArealData::oneOverI][i];
              Real p_oneOverI_new(p_oneOverI_old);

              // Get current particle's volume
              const Real p_vol = p_realarray[SoArealData::volume][i];

              // Flag to stop computing particle's quantities if mass_new < 0,
              // i.e. the particle disappears because of chemical reactions
              int proceed = 1;

              //***************************************************************
              // First step: update particles' mass and density
              //***************************************************************
              if(local_update_mass)
              {
                // Total particle density exchange rate
                Real total_mass_rate(0);

                for (int n_s(0); n_s < nspecies_s; ++n_s) {

                  // Current species mass fraction
                  X_sn[n_s] = ptile_data.m_runtime_rdata[idx_X_sn+n_s][i];

                  // Get the current reaction rate for species n_s
                  const Real mass_sn_rate = ptile_data.m_runtime_rdata[idx_mass_sn_txfr+n_s][i];

                  X_sn[n_s] = X_sn[n_s]*p_mass_old + subdt*mass_sn_rate;

                  // Update the total mass exchange rate
                  total_mass_rate += mass_sn_rate;
                }

                // Update the total mass of the particle
                p_mass_new = p_mass_old + subdt * total_mass_rate;

                if (p_mass_new > 0) {

                  Real total_X(0.);

                  // Normalize species mass fractions
                  for (int n_s(0); n_s < nspecies_s; n_s++) {
                    Real X_sn_new = X_sn[n_s] / p_mass_new;

                    if (X_sn_new < 0) X_sn_new = 0;
                    if (X_sn_new > 1) X_sn_new = 1;

                    total_X += X_sn_new;
                    X_sn[n_s] = X_sn_new;
                  }

                  for (int n_s(0); n_s < nspecies_s; n_s++) {
                    // Divide updated species mass fractions by total_X
                    X_sn[n_s] /= total_X;
                    ptile_data.m_runtime_rdata[idx_X_sn+n_s][i] = X_sn[n_s];
                  }

                  // Write out to global memory particle's mass and density
                  p_realarray[SoArealData::mass][i] = p_mass_new;
                  p_density_new = p_mass_new / p_vol;
                  p_realarray[SoArealData::density][i] = p_density_new;

                  // Write out to global memory particle's moment of inertia
                  p_oneOverI_new = (p_density_old/p_density_new)*p_oneOverI_old;
                  p_realarray[SoArealData::oneOverI][i] = p_oneOverI_new;

                } else {
                  p.id() = -1;
                  proceed = 0;
                }
              }

              if (proceed) {
                //***************************************************************
                // Second step: update particles' positions and velocities
                //***************************************************************
                if (update_momentum) {
                  const Real p_velx_old = p_realarray[SoArealData::velx][i];
                  const Real p_vely_old = p_realarray[SoArealData::vely][i];
                  const Real p_velz_old = p_realarray[SoArealData::velz][i];

                  const Real vel_coeff = local_update_mass ? p_mass_old/p_mass_new : 1.;

                  Real p_velx_new = vel_coeff*p_velx_old +
                    subdt*((p_realarray[SoArealData::dragx][i]+fc_ptr[i]) / p_mass_new + vel_coeff*gravity[0]);
                  Real p_vely_new = vel_coeff*p_vely_old +
                    subdt*((p_realarray[SoArealData::dragy][i]+fc_ptr[i+ntot]) / p_mass_new + vel_coeff*gravity[1]);
                  Real p_velz_new = vel_coeff*p_velz_old +
                    subdt*((p_realarray[SoArealData::dragz][i]+fc_ptr[i+2*ntot]) / p_mass_new + vel_coeff*gravity[2]);

                  if (local_solve_reactions) {
                    p_velx_new += subdt*(ptile_data.m_runtime_rdata[idx_vel_s_txfr+0][i] / p_mass_new);
                    p_vely_new += subdt*(ptile_data.m_runtime_rdata[idx_vel_s_txfr+1][i] / p_mass_new);
                    p_velz_new += subdt*(ptile_data.m_runtime_rdata[idx_vel_s_txfr+2][i] / p_mass_new);
                  }

                  const Real p_omegax_old = p_realarray[SoArealData::omegax][i];
                  const Real p_omegay_old = p_realarray[SoArealData::omegay][i];
                  const Real p_omegaz_old = p_realarray[SoArealData::omegaz][i];

                  const Real omega_coeff = local_update_mass ? p_oneOverI_new/p_oneOverI_old : 1.;

                  Real p_omegax_new = omega_coeff*p_omegax_old + subdt * p_oneOverI_new * tow_ptr[i];
                  Real p_omegay_new = omega_coeff*p_omegay_old + subdt * p_oneOverI_new * tow_ptr[i+ntot];
                  Real p_omegaz_new = omega_coeff*p_omegaz_old + subdt * p_oneOverI_new * tow_ptr[i+2*ntot];

                  const Real p_posx_old = p.pos(0);
                  const Real p_posy_old = p.pos(1);
                  const Real p_posz_old = p.pos(2);

                  Real p_posx_new = p_posx_old + subdt * p_velx_new;
                  Real p_posy_new = p_posy_old + subdt * p_vely_new;
                  Real p_posz_new = p_posz_old + subdt * p_velz_new;

                  if (x_lo_bc && p_posx_new < p_lo[0])
                  {
                      p_posx_new = p_lo[0] + eps;
                      p_velx_new = -p_velx_new;
                  }
                  if (x_hi_bc && p_posx_new > p_hi[0])
                  {
                      p_posx_new = p_hi[0] - eps;
                      p_velx_new = -p_velx_new;
                  }
                  if (y_lo_bc && p_posy_new < p_lo[1])
                  {
                      p_posy_new = p_lo[1] + eps;
                      p_vely_new = -p_vely_new;
                  }
                  if (y_hi_bc && p_posy_new > p_hi[1])
                  {
                      p_posy_new = p_hi[1] - eps;
                      p_vely_new = -p_vely_new;
                  }
                  if (z_lo_bc && p_posz_new < p_lo[2])
                  {
                      p_posz_new = p_lo[2] + eps;
                      p_velz_new = -p_velz_new;
                  }
                  if (z_hi_bc && p_posz_new > p_hi[2])
                  {
                      p_posz_new = p_hi[2] - eps;
                      p_velz_new = -p_velz_new;
                  }

                  // Update positions
                  p.pos(0) = p_posx_new;
                  p.pos(1) = p_posy_new;
                  p.pos(2) = p_posz_new;

                  // Update velocities
                  p_realarray[SoArealData::velx][i] = p_velx_new;
                  p_realarray[SoArealData::vely][i] = p_vely_new;
                  p_realarray[SoArealData::velz][i] = p_velz_new;

                  // Update angular velocities
                  p_realarray[SoArealData::omegax][i] = p_omegax_new;
                  p_realarray[SoArealData::omegay][i] = p_omegay_new;
                  p_realarray[SoArealData::omegaz][i] = p_omegaz_new;
                }

                //***************************************************************
                // Third step: update particles' temperature
                //***************************************************************
                if (local_update_enthalpy) {

                  const int phase = p_intarray[SoAintData::phase][i];

                  const Real Tp_old = p_realarray[SoArealData::temperature][i];

                  const Real coeff = local_update_mass ? (p_mass_old/p_mass_new) : 1.;

                  Real p_enthalpy_old(0);

                  if (solid_is_a_mixture) {
                    for (int n_s(0); n_s < nspecies_s; ++n_s) {
                      p_enthalpy_old += X_sn[n_s]*solids_parms.calc_h_sn<RunOn::Gpu>(Tp_old,n_s);
                    }
                  } else {
                    p_enthalpy_old = solids_parms.calc_h_s<RunOn::Gpu>(phase-1,Tp_old);
                  }

                  Real p_enthalpy_new = coeff*p_enthalpy_old +
                    subdt*((p_realarray[SoArealData::convection][i]+enthalpy_source) / p_mass_new);

                  if (local_solve_reactions) {
                    p_enthalpy_new -= subdt*(ptile_data.m_runtime_rdata[idx_h_s_txfr][i] / p_mass_new);
                  }

                  // ************************************************************
                  // Newton-Raphson solver for solving implicit equation for
                  // temperature
                  // ************************************************************
                  // Residual computation
                  auto R = [&] AMREX_GPU_DEVICE (Real Tp_arg)
                  {
                    Real hp_loc(0);

                    if (!solid_is_a_mixture) {

                      hp_loc = solids_parms.calc_h_s<RunOn::Gpu>(phase-1,Tp_arg);
                    } else {

                      for (int n_s(0); n_s < nspecies_s; ++n_s)
                        hp_loc += X_sn[n_s]*solids_parms.calc_h_sn<RunOn::Gpu>(Tp_arg,n_s);
                    }

                    return hp_loc - p_enthalpy_new;
                  };

                  // Partial derivative computation
                  auto partial_R = [&] AMREX_GPU_DEVICE (Real Tp_arg)
                  {
                    Real gradient(0);

                    if (!solid_is_a_mixture) {

                      gradient = solids_parms.calc_partial_h_s<RunOn::Gpu>(phase-1,Tp_arg);
                    } else {

                      for (int n_s(0); n_s < nspecies_s; ++n_s)
                        gradient += X_sn[n_s]*solids_parms.calc_partial_h_sn<RunOn::Gpu>(Tp_arg,n_s);
                    }

                    return gradient;
                  };

                  Real Tp_new(Tp_old);

                  const Real dumping_factor = 1.;

                  DampedNewton::solve(Tp_new, R, partial_R, dumping_factor, 1.e-6, 1.e-6);

                  p_realarray[SoArealData::temperature][i] = Tp_new;

                  // Update cp_s
                  Real cp_s_new(0);

                  if (solid_is_a_mixture) {
                    for (int n_s(0); n_s < nspecies_s; ++n_s)
                      cp_s_new += X_sn[n_s]*solids_parms.calc_cp_sn<RunOn::Gpu>(Tp_new,n_s);

                  } else {
                    cp_s_new = solids_parms.calc_cp_s<RunOn::Gpu>(phase-1,Tp_new);
                  }

                  AMREX_ASSERT(cp_s_new > 0.);
                  p_realarray[SoArealData::cp_s][i] = cp_s_new;
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

        /************************************************************************
         * DEBUG: output the number of collisions in current substep            *
         *        output the max velocity (and forces) in current substep       *
         *        update max velocities and forces                              *
         ***********************************************************************/

        if (debug_level > 0) ncoll_total += ncoll;

        if (debug_level > 1) {
            ParallelDescriptor::ReduceIntSum(ncoll, ParallelDescriptor::IOProcessorNumber());
            Print() << "Number of collisions: " << ncoll << " at step " << n << std::endl;
        }

        if (debug_level > 0){
            UpdateMaxVelocity();
            UpdateMaxForces(pfor, wfor);
        }

        if (debug_level > 1) {
            RealVect max_vel = GetMaxVelocity();
            Vector<RealVect> max_forces = GetMaxForces();

            const Real * dx_crse = Geom(0).CellSize();
            amrex::Print() << "Maximum distance traveled:"
                           << std::endl
                           <<  "x= " << max_vel[0] * dt
                           << " y= " << max_vel[1] * dt
                           << " z= " << max_vel[2] * dt
                           << " and note that "
                           << " dx= " << dx_crse[0] << std::endl;

            amrex::Print() << "Maximum particle-particle (pp) and particle-wall (pw) forces:"
                           << std::endl
                           <<  "ppx= " << max_forces[0][0]
                           << " ppy= " << max_forces[0][1]
                           << " ppz= " << max_forces[0][2] << std::endl
                           <<  "pwx= " << max_forces[1][0]
                           << " pwy= " << max_forces[1][1]
                           << " pwz= " << max_forces[1][2] << std::endl;

        }

    } // end of loop over substeps

    // Redistribute particles at the end of all substeps (note that the particle
    // neighbour list needs to be reset when redistributing).
    if (update_momentum) {
      clearNeighbors();
      Redistribute(0, 0, 0, 1);
    }

    /****************************************************************************
     * DEBUG: output the total number of collisions over all substeps           *
     *        output the maximum velocity and forces over all substeps          *
     ***************************************************************************/
    if (debug_level > 0) {
        ParallelDescriptor::ReduceIntSum(ncoll_total, ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "Number of collisions: " << ncoll_total << " in "
                       << nsubsteps << " substeps " << std::endl;
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

    if (debug_level > 0) {
        RealVect max_vel = GetMaxVelocity();
        Vector<RealVect> max_forces = GetMaxForces();

        const Real * dx_crse = Geom(0).CellSize();
        amrex::Print() << "Maximum possible distance traveled:" << std::endl
                       <<  "x= " << max_vel[0] * dt
                       << " y= " << max_vel[1] * dt
                       << " z= " << max_vel[2] * dt
                       << " and note that "
                       << " dx= " << dx_crse[0] << std::endl;

        amrex::Print() << "Maximum particle-particle (pp) and particle-wall (pw) forces:" << std::endl
                       <<  "ppx= " << max_forces[0][0]
                       << " ppy= " << max_forces[0][1]
                       << " ppz= " << max_forces[0][2] << std::endl
                       <<  "pwx= " << max_forces[1][0]
                       << " pwy= " << max_forces[1][1]
                       << " pwz= " << max_forces[1][2] << std::endl;
    }

    amrex::Print() << "done. \n";

    BL_PROFILE_REGION_STOP("mfix_dem::EvolveParticles()");
}

void MFIXParticleContainer::writeAllAtLevel (int lev)
{
    // Not threaded because its print to terminal
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        int np = pti.numParticles();
        Gpu::HostVector<ParticleType> host_particles(np);
        Gpu::copy(Gpu::deviceToHost, particles.begin(), particles.end(), host_particles.begin());

        for (const auto& p: host_particles)
        {
           const IntVect& iv = Index(p, lev);

           RealVect xyz(p.pos(0), p.pos(1), p.pos(2));
           std::cout << " id " << p.id()
                << " index " << iv
                << " position " << xyz << std::endl;
       }
    }
}

void
MFIXParticleContainer::WriteAsciiFileForInit (const std::string& filename)
{
    BL_ASSERT(!filename.empty());

    int lev = 0;
    long nparticles = NumberOfParticlesAtLevel(lev);

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Have I/O processor open file and write out particle metadata.
        //
        std::ofstream File;

        File.open(filename.c_str(), std::ios::out | std::ios::trunc);

        if (!File.good())
            amrex::FileOpenFailed(filename);

        File << nparticles  << '\n';

        File.flush();

        File.close();

        if (!File.good())
            amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::WriteAsciiFile(): problem writing file");
    }

    ParallelDescriptor::Barrier();

    const int MyProc = ParallelDescriptor::MyProc();

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (MyProc == i)
        {
            //
            // Each CPU opens the file for appending and adds its particles.
            //

            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            std::ofstream File;

            File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            File.open(filename.c_str(), std::ios::out|std::ios::app);

            File.precision(15);

            if (!File.good())
                amrex::FileOpenFailed(filename);

            for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

              auto& particles = pti.GetArrayOfStructs();
              int np = pti.numParticles();

              auto& soa = pti.GetStructOfArrays();
              auto p_realarray = soa.realarray();
              auto p_intarray = soa.intarray();

              std::array<Gpu::HostVector<Real>, SoArealData::count> host_realarrays;
              std::array<Gpu::HostVector<int>, SoAintData::count> host_intarrays;

              for (int comp(0); comp < SoArealData::count; ++comp)
                host_realarrays[comp].resize(np);

              for (int comp(0); comp < SoAintData::count; ++comp)
                host_intarrays[comp].resize(np);

              // Copy particles from device to host
              for (int comp(0); comp < SoArealData::count; ++comp) {
                Gpu::copyAsync(Gpu::deviceToHost, &(p_realarray[comp][0]),
                    &(p_realarray[comp][np]), host_realarrays[comp].begin());
              }

              // Copy particles from device to host
              for (int comp(0); comp < SoAintData::count; ++comp) {
                Gpu::copyAsync(Gpu::deviceToHost, &(p_intarray[comp][0]),
                    &(p_intarray[comp][np]), host_intarrays[comp].begin());
              }

              Gpu::HostVector<ParticleType> host_particles(np);
              Gpu::copy(Gpu::deviceToHost, particles.begin(), particles.end(), host_particles.begin());

              int index = 0;
              for (int ip(0); ip < np; ++ip)
              {
                  auto& p = host_particles[ip];

                  if (p.id() > 0) {
                      File << p_intarray[SoAintData::phase][ip] << ' ';
                      File << p.pos(0) << ' ';
                      File << p.pos(1) << ' ';
                      File << p.pos(2) << ' ';
                      File << p_realarray[SoArealData::radius][ip] << ' ';
                      File << p_realarray[SoArealData::density][ip] << ' ';
                      File << p_realarray[SoArealData::velx][ip] << ' ';
                      File << p_realarray[SoArealData::vely][ip] << ' ';
                      File << p_realarray[SoArealData::velz][ip] << ' ';

                      File << '\n';

                      index++;
                  }
              }
            }

            File.flush();

            File.close();

            if (!File.good())
                amrex::Abort("MFIXParticleContainer::WriteAsciiFileForInit(): problem writing file");

        }
        ParallelDescriptor::Barrier();
    }
}



void MFIXParticleContainer::UpdateMaxVelocity ()
{
  Real max_vel_x = loc_maxvel[0];
  Real max_vel_y = loc_maxvel[1];
  Real max_vel_z = loc_maxvel[2];

  for (int lev = 0; lev < nlev; lev++)
  {
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
      // Reduce max operation for velx, vely,velz
      ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpMax> reduce_op;
      ReduceData<Real, Real, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      for(MFIXParIter pti(* this, lev); pti.isValid(); ++ pti)
      {
        const int np = pti.numParticles();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();

        reduce_op.eval(np, reduce_data,
            [p_realarray] AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
        {
          Real l_vel_x = Math::abs(p_realarray[SoArealData::velx][p_id]);
          Real l_vel_y = Math::abs(p_realarray[SoArealData::vely][p_id]);
          Real l_vel_z = Math::abs(p_realarray[SoArealData::velz][p_id]);

          return {l_vel_x, l_vel_y, l_vel_z};
        });
      }

      ReduceTuple host_tuple = reduce_data.value();
      max_vel_x = amrex::max(max_vel_x, Math::abs(amrex::get<0>(host_tuple)));
      max_vel_y = amrex::max(max_vel_y, Math::abs(amrex::get<1>(host_tuple)));
      max_vel_z = amrex::max(max_vel_z, Math::abs(amrex::get<2>(host_tuple)));
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(max:max_vel_x,max_vel_y,max_vel_z) \
                              if (Gpu::notInLaunchRegion())
#endif
      for(MFIXParIter pti(* this, lev); pti.isValid(); ++ pti)
      {
        const int np = pti.numParticles();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();

        for(int p_id(0); p_id < np; ++p_id)
        {
          max_vel_x = amrex::max(Math::abs(p_realarray[SoArealData::velx][p_id]), max_vel_x);
          max_vel_y = amrex::max(Math::abs(p_realarray[SoArealData::vely][p_id]), max_vel_y);
          max_vel_z = amrex::max(Math::abs(p_realarray[SoArealData::velz][p_id]), max_vel_z);
        }
      }
    }
  }

  loc_maxvel = RealVect(max_vel_x, max_vel_y, max_vel_z);
}

void MFIXParticleContainer::UpdateMaxForces (std::map<PairIndex, Gpu::DeviceVector<Real>>& pfor,
                                             std::map<PairIndex, Gpu::DeviceVector<Real>>& wfor)
{
  Real max_pfor_x = loc_maxpfor[0];
  Real max_pfor_y = loc_maxpfor[1];
  Real max_pfor_z = loc_maxpfor[2];

  Real max_wfor_x = loc_maxwfor[0];
  Real max_wfor_y = loc_maxwfor[1];
  Real max_wfor_z = loc_maxwfor[2];

  for (int lev = 0; lev < nlev; lev++)
  {
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
      // Reduce max operation for pforx, pfory, pforz, wforx, wfory, wforz
      ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpMax,
                ReduceOpMax, ReduceOpMax, ReduceOpMax> reduce_op;
      ReduceData<Real, Real, Real, Real, Real, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      for(MFIXParIter pti(* this, lev); pti.isValid(); ++ pti)
      {
        PairIndex index(pti.index(), pti.LocalTileIndex());

        // Note the particle force data layout:
        //      p1_x, p2_x, ..., pn_x, p1_y, p2_y, ..., pn_y, p1_z, p2_z, ..., pn_z
        // Where n is the total number of particle and neighbor particles.
        const int nrp = GetParticles(lev)[index].numRealParticles();

        auto& plev = GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        int size_ng = aos.numNeighborParticles();

        // Number of particles including neighbor particles
        const int ntot = nrp + size_ng;

        Real* p_pfor = pfor[index].data();
        Real* p_wfor = wfor[index].data();

        // Find max (abs) of particle-particle forces:
        reduce_op.eval(ntot, reduce_data, [p_pfor,p_wfor,ntot]
          AMREX_GPU_DEVICE (int i) -> ReduceTuple
        {
          Real l_pfor_x = Math::abs(p_pfor[i]);
          Real l_pfor_y = Math::abs(p_pfor[i+ntot]);
          Real l_pfor_z = Math::abs(p_pfor[i+2*ntot]);

          Real l_wfor_x = Math::abs(p_wfor[i]);
          Real l_wfor_y = Math::abs(p_wfor[i+ntot]);
          Real l_wfor_z = Math::abs(p_wfor[i+2*ntot]);

          return {l_pfor_x, l_pfor_y, l_pfor_z, l_wfor_x, l_wfor_y, l_wfor_z};
        });
      }

      ReduceTuple host_tuple = reduce_data.value();
      
      max_pfor_x = amrex::max(max_pfor_x, Math::abs(amrex::get<0>(host_tuple)));
      max_pfor_y = amrex::max(max_pfor_y, Math::abs(amrex::get<1>(host_tuple)));
      max_pfor_z = amrex::max(max_pfor_z, Math::abs(amrex::get<2>(host_tuple)));
      
      max_pfor_x = amrex::max(max_pfor_x, Math::abs(amrex::get<3>(host_tuple)));
      max_pfor_y = amrex::max(max_pfor_y, Math::abs(amrex::get<4>(host_tuple)));
      max_pfor_z = amrex::max(max_pfor_z, Math::abs(amrex::get<5>(host_tuple)));
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(max:max_pfor_x,max_pfor_y,max_pfor_z, \
                                   max_wfor_x,max_wfor_y,max_wfor_z) \
                              if (Gpu::notInLaunchRegion())
#endif
      for(MFIXParIter pti(* this, lev); pti.isValid(); ++ pti)
      {
        PairIndex index(pti.index(), pti.LocalTileIndex());

        // Note the particle force data layout:
        //      p1_x, p2_x, ..., pn_x, p1_y, p2_y, ..., pn_y, p1_z, p2_z, ..., pn_z
        // Where n is the total number of particle and neighbor particles.
        const int nrp = GetParticles(lev)[index].numRealParticles();

        int size_ng = neighbors[lev][index].size();
        
        // Number of particles including neighbor particles
        const int ntot = nrp + size_ng;

        Real* p_pfor = pfor[index].data();
        Real* p_wfor = wfor[index].data();

        // Find max (abs) of particle-particle forces:
        for (int i(0); i < ntot; ++i)
        {
          max_pfor_x = amrex::max(Math::abs(p_pfor[i]), max_pfor_x);
          max_pfor_y = amrex::max(Math::abs(p_pfor[i+ntot]), max_pfor_y);
          max_pfor_z = amrex::max(Math::abs(p_pfor[i+2*ntot]), max_pfor_z);

          max_wfor_x = amrex::max(Math::abs(p_wfor[i]), max_wfor_x);
          max_wfor_y = amrex::max(Math::abs(p_wfor[i+ntot]), max_wfor_y);
          max_wfor_z = amrex::max(Math::abs(p_wfor[i+2*ntot]), max_wfor_z);
        }
      }
    }
  }

  loc_maxpfor = RealVect(max_pfor_x, max_pfor_y, max_pfor_z);
  loc_maxwfor = RealVect(max_wfor_x, max_wfor_y, max_wfor_z);
}

RealVect MFIXParticleContainer::GetMaxVelocity ()
{
    Real max_vel_x = loc_maxvel[0], max_vel_y = loc_maxvel[1], max_vel_z = loc_maxvel[2];

    ParallelDescriptor::ReduceRealMax({max_vel_x, max_vel_y, max_vel_z},
                                      ParallelDescriptor::IOProcessorNumber());

    RealVect max_vel(max_vel_x, max_vel_y, max_vel_z);

    return max_vel;
}

Vector<RealVect> MFIXParticleContainer::GetMaxForces ()
{
    Real max_pfor_x = loc_maxpfor[0], max_pfor_y = loc_maxpfor[1], max_pfor_z = loc_maxpfor[2];
    Real max_wfor_x = loc_maxwfor[0], max_wfor_y = loc_maxwfor[1], max_wfor_z = loc_maxwfor[2];


    ParallelDescriptor::ReduceRealMax({max_pfor_x, max_pfor_y, max_pfor_z,
                                       max_wfor_x, max_wfor_y, max_wfor_z},
                                      ParallelDescriptor::IOProcessorNumber());

    Vector<RealVect> max_forces(2);
    max_forces[0] = RealVect(max_pfor_x, max_pfor_y, max_pfor_z);
    max_forces[1] = RealVect(max_wfor_x, max_wfor_y, max_wfor_z);

    return max_forces;
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

    std::vector<long> region_np(nregions, 0);
    std::vector<Real> region_velx(nregions, 0.0);
    std::vector<Real> region_vely(nregions, 0.0);
    std::vector<Real> region_velz(nregions, 0.0);
    std::vector<Real> region_k_en(nregions, 0.0);

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

      long sum_np   = 0;    // Number of particle in avg region
      Real sum_velx = 0.;
      Real sum_vely = 0.;
      Real sum_velz = 0.;
      Real sum_k_en = 0.;

#ifdef AMREX_USE_GPU
      if (Gpu::inLaunchRegion())
      {
        // Reduce sum operation for np, velx, vely, velz, kinetic energy
        ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
        ReduceData<long, Real, Real, Real, Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

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
        sum_np   = amrex::get<0>(host_tuple);
        sum_velx = amrex::get<1>(host_tuple);
        sum_vely = amrex::get<2>(host_tuple);
        sum_velz = amrex::get<3>(host_tuple);
        sum_k_en = amrex::get<4>(host_tuple);
      }
      else
#endif
      {
#ifdef _OPENMP
#pragma omp parallel reduction(+:sum_np,sum_velx,sum_vely,sum_velz, \
                                 sum_k_en) if (Gpu::notInLaunchRegion())
#endif
        for (MFIXParIter pti(*this,lev); pti.isValid(); ++pti)
        {
          Box bx = pti.tilebox();
          RealBox tile_region(bx, Geom(lev).CellSize(), Geom(lev).ProbLo());

          if (tile_region.intersects(avg_region))
          {
            const int np         = NumberOfParticles(pti);
            const AoS &particles = pti.GetArrayOfStructs();
            const ParticleType* pstruct = particles().dataPtr();

            SoA& soa = pti.GetStructOfArrays();
            auto p_realarray = soa.realarray();

            for (int p_id(0); p_id < np; ++p_id)
            {
              const ParticleType p = pstruct[p_id];

              if (avg_region.contains(p.pos()))
              {
                const Real mass = p_realarray[SoArealData::mass][p_id];
                const Real velx = p_realarray[SoArealData::velx][p_id];
                const Real vely = p_realarray[SoArealData::vely][p_id];
                const Real velz = p_realarray[SoArealData::velz][p_id];
                const Real k_en = 0.5*mass*(velx*velx + vely*vely + velz*velz);

                sum_np   += static_cast<long>(1);
                sum_velx += velx;
                sum_vely += vely;
                sum_velz += velz;
                sum_k_en += k_en;
              }
            }
          }
        }
      }

      region_np[nr]   = sum_np;
      region_velx[nr] = sum_velx;
      region_vely[nr] = sum_vely;
      region_velz[nr] = sum_velz;
      region_k_en[nr] = sum_k_en;
    }

    // Compute parallel reductions
    ParallelDescriptor::ReduceLongSum(region_np.data(),   nregions);
    ParallelDescriptor::ReduceRealSum(region_velx.data(), nregions);
    ParallelDescriptor::ReduceRealSum(region_vely.data(), nregions);
    ParallelDescriptor::ReduceRealSum(region_velz.data(), nregions);
    ParallelDescriptor::ReduceRealSum(region_k_en.data(), nregions);

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

    std::vector<long> region_np(nregions, 0);
    std::vector<Real> region_Tp(nregions, 0.0);

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

      long sum_np = 0;    // Number of particle in avg region
      Real sum_Tp = 0.;

#ifdef AMREX_USE_GPU
      if (Gpu::inLaunchRegion())
      {
        // Reduce sum operation for np, Tp
        ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
        ReduceData<long, Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

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
        sum_np = amrex::get<0>(host_tuple);
        sum_Tp = amrex::get<1>(host_tuple);
      }
      else
#endif
      {
#ifdef _OPENMP
#pragma omp parallel reduction(+:sum_np,sum_Tp) if (Gpu::notInLaunchRegion())
#endif
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++ pti)
        {
          Box bx = pti.tilebox();
          RealBox tile_region(bx, Geom(lev).CellSize(), Geom(lev).ProbLo());

          if (tile_region.intersects(avg_region))
          {
            const int np          = NumberOfParticles(pti);
            const AoS& particles  = pti.GetArrayOfStructs();
            const ParticleType* pstruct = particles().dataPtr();

            auto& soa = pti.GetStructOfArrays();
            auto p_realarray = soa.realarray();

            for (int p_id(0); p_id < np; ++p_id)
            {
              const ParticleType& p = pstruct[p_id];

              if (avg_region.contains(p.pos()))
              {
                sum_np += static_cast<long>(1);
                sum_Tp += p_realarray[SoArealData::temperature][p_id];
              }
            }
          }
        }
      }

      region_np[nr] = sum_np;
      region_Tp[nr] = sum_Tp;
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

void MFIXParticleContainer::set_particle_properties (int /*pstate*/,
                                                     Real pradius,
                                                     Real pdensity,
                                                     Real& pvol,
                                                     Real& pmass,
                                                     Real& omoi,
                                                     Real& omega)
{
    pvol  = (4.0/3.0)*M_PI*(pradius*pradius*pradius);
    pmass = pvol * pdensity;
    omoi  = 2.5/(pmass * (pradius*pradius));
    omega = 0.0;
}


void MFIXParticleContainer::checkParticleBoxSize(int      lev, 
                                                 IntVect& loc_max_grid_size, 
                                                 Real     frac_particle_bin)
{
  // count total # particles
  if (m_total_numparticle <= 0) {
    m_total_numparticle = 0;
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      m_total_numparticle += pti.numRealParticles();
    ParallelDescriptor::ReduceIntSum(m_total_numparticle);
  }

  // average # particles per process
  Real avg_np = static_cast<Real>(m_total_numparticle)
              / ParallelDescriptor::NProcs()
              * frac_particle_bin;
  // Print() << "total # particles " << m_total_numparticle 
          // << " average # particles " << avg_np << std::endl;

  // 
  IntVect ncut(AMREX_D_DECL(4, 4, 4));
  IntVect bin_size(AMREX_D_DECL(
            amrex::max(loc_max_grid_size[0]/ncut[0], 1), 
            amrex::max(loc_max_grid_size[1]/ncut[1], 1),
            amrex::max(loc_max_grid_size[2]/ncut[2], 1)));

  // dictoinary for # particles per bin
  std::map<PairIndex, Gpu::DeviceVector<int> > np_bin;
  std::map<PairIndex, IntVect> tile_max_bin_size;

  // particle tiles and geometry of this level
  const auto& geom = Geom(lev);
  const auto  dx_inv = geom.InvCellSizeArray();
  const auto  prob_lo = geom.ProbLoArray();

  auto& plev = GetParticles(lev);

  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {
    PairIndex   index(pti.index(), pti.LocalTileIndex());

    // particle box, range
    const auto& pbox  = pti.tilebox();
    const auto  lo    = lbound(pbox);
    Real nx  = static_cast<Real>(pbox.length(0));
    Real ny  = static_cast<Real>(pbox.length(1));
    Real nz  = static_cast<Real>(pbox.length(2));
    IntVect nbin(AMREX_D_DECL(
              static_cast<int>(Math::ceil(nx/bin_size[0])),
              static_cast<int>(Math::ceil(ny/bin_size[1])),
              static_cast<int>(Math::ceil(nz/bin_size[2]))));

    const auto& ptile = plev[index];
    int np = ptile.numRealParticles();
    auto& aos = ptile.GetArrayOfStructs();
    const auto* pstruct = aos().dataPtr();

    np_bin[index].clear();
    np_bin[index].resize(nbin[0]*nbin[1]*nbin[2], 0);
    auto np_bin_ptr = np_bin[index].dataPtr();

    ParallelFor(np, [pstruct, prob_lo, dx_inv, lo, np_bin_ptr, bin_size, nbin]
      AMREX_GPU_DEVICE (int i) noexcept
      {
        const auto p = pstruct[i];
        int bin_i = ( static_cast<int>(Math::floor((p.pos(0) - prob_lo[0]) * dx_inv[0])) 
                    - lo.x) / bin_size[0];
        int bin_j = ( static_cast<int>(Math::floor((p.pos(1) - prob_lo[1]) * dx_inv[1])) 
                    - lo.y) / bin_size[1];
        int bin_k = ( static_cast<int>(Math::floor((p.pos(2) - prob_lo[2]) * dx_inv[2])) 
                    - lo.z) / bin_size[2];
        Gpu::Atomic::Add(np_bin_ptr + bin_i + bin_j*nbin[0] + bin_k*nbin[0]*nbin[1], 1);
      });// end parallel for

    Gpu::streamSynchronize();

    IntVect incr(IntVect::TheUnitVector());
    int nincr = 0;            // # increments done
    int nfail = 0;            // # fails to increment bin size
    while (incr[0] < nbin[0] || incr[1] < nbin[1] || incr[2] < nbin[2]) {
      // max # particles among merged bins
      int max_np_bin_merge = 0; 
      for (int bin_i=0; bin_i<nbin[0]; bin_i+=incr[0]) {
        for (int bin_j=0; bin_j<nbin[1]; bin_j+=incr[1]) {
          for (int bin_k=0; bin_k<nbin[2]; bin_k+=incr[2]) {
            // # particles in merged bins
            int np_bin_merge = 0;
            for (int ii=bin_i; ii<min(bin_i+incr[0], nbin[0]); ++ii)
              for (int jj=bin_j; jj<min(bin_j+incr[1], nbin[1]); ++jj)
                for (int kk=bin_k; kk<min(bin_k+incr[2], nbin[2]); ++kk)
                  np_bin_merge += np_bin_ptr[ii+jj*nbin[0]+kk*nbin[0]*nbin[1]];
            max_np_bin_merge = max(max_np_bin_merge, np_bin_merge);
          }
        }
      } // end for bin_i

      // if exceeds the allowed average, then undo the increment 
      // and try it in another direction
      if (max_np_bin_merge > avg_np) {
        nfail += 1;
        incr[2-(nincr-1)%3] /= 2;
      }

      if (nfail == 3)  break;

      // increment by 2x, exclude incr[idim] == nbin[idim]
      int idim = 2 - nincr % 3;
      if (incr[idim] < nbin[idim])  incr[idim] *= 2;
      ++nincr;
    } // end while

    for (int i=0; i<AMREX_SPACEDIM; ++i)
      tile_max_bin_size[index][i] = bin_size[i] * min(incr[i], nbin[i]);

  }// end for

  IntVect max_bin_size(AMREX_D_DECL(1024, 1024, 1024));
  for (auto key_val: tile_max_bin_size)
    max_bin_size.min(key_val.second);

  ParallelDescriptor::ReduceIntMin(max_bin_size.begin(), AMREX_SPACEDIM);

  //for (auto key_val: np_bin) {
    //Print() << "index pair " << key_val.first << " # particles per bin ";
    //for (auto n: key_val.second)  Print() << n << " ";
    //Print() << "\n";
  //}
  // Print() << "new max grid size " << max_bin_size << std::endl;
}


void MFIXParticleContainer::downsizeParticleBoxes(int lev, 
                                                  const IntVect& new_grid_size)
{
  const auto         ba   = ParticleBoxArray(lev);
  const auto         dmap = ParticleDistributionMap(lev);
  const Vector<int>& pmap = dmap.ProcessorMap();

  // create boxarray with input sizes
  const auto& geom   = Geom(lev);
  const auto& domain = geom.Domain();
  IntVect domainlen(domain.size());
  BoxArray new_ba(domain);
  new_ba.maxSize(new_grid_size);

  // number of new boxes in each dimension
  IntVect new_nbox(1);
  for (int i=0; i<AMREX_SPACEDIM; ++i)
    new_nbox[i] = domainlen[i] % new_grid_size[i] == 0 ? \
                  domainlen[i] / new_grid_size[i] : 
                  domainlen[i] / new_grid_size[i] + 1;

  // decomposition of current grids
  IntVect nbox(1);
  for (int i=0; i<AMREX_SPACEDIM; ++i)
    nbox[i] = domainlen[i] % max_grid_size[i] == 0 ? \
              domainlen[i] / max_grid_size[i] : 
              domainlen[i] / max_grid_size[i] + 1;

  // property of first box and first new box
  int new_box_idx = 0;
  Box     new_box    = new_ba[new_box_idx];
  Dim3    new_boxhi  = ubound(new_box);
  IntVect new_boxlen = new_box.size();
  int box_idx = 0, box_k = 0, box_j = 0;
  Dim3    boxhi      = ubound(ba[box_idx]);

  // set the mapping of new boxes 
  // note that the boxes are stored following i->j->k order
  Vector<int> new_pmap;
  for (int k=0; k<new_nbox[2]; ++k) {
    if (new_boxhi.z >= boxhi.z + new_boxlen[2]/2)  
      ++box_k;
    for (int j=0; j<new_nbox[1]; ++j) {
      if (new_boxhi.y >= boxhi.y + new_boxlen[1]/2)  
        ++box_j;
      box_idx = box_j*nbox[0] + box_k*nbox[0]*nbox[1]; 
      boxhi = ubound(ba[box_idx]);
      for (int i=0; i<new_nbox[0]; ++i) {
        if (new_boxhi.x >= boxhi.x + new_boxlen[0]/2) {
          ++box_idx;
          boxhi = ubound(ba[box_idx]);
        }
        new_pmap.push_back(pmap[box_idx]);
        ++new_box_idx;
        new_box    = new_ba[new_box_idx];
        new_boxhi  = ubound(new_box);
        new_boxlen = new_box.size();
      }// end new box i
    }// end new box j
  }// end new box k

  // new distribution mapping
  DistributionMapping new_dmap(new_pmap);

  // ba and dmap to particle container
  SetParticleBoxArray(lev, new_ba);
  SetParticleDistributionMap(lev, new_dmap);

  max_grid_size = new_grid_size;
}


void MFIXParticleContainer::resetCostByCount(int lev, 
                                             Vector<MultiFab*>& cost)
{
  // properties of current particle boxes
  // domain length
  const auto&  geom   = Geom(lev); 
  const auto&  domain = geom.Domain();
  IntVect domainlen(domain.size());
  // # boxes in each dimension
  IntVect nbox(1);
  for (int i=0; i<AMREX_SPACEDIM; ++i)
    nbox[i] = domainlen[i] % max_grid_size[i] == 0 ? \
              domainlen[i] / max_grid_size[i] : 
              domainlen[i] / max_grid_size[i] + 1;
  // box sizes in each dimension
  // Each sizes has two values, nlarge large values followed by 
  // small values
  const auto&  ba = ParticleBoxArray(lev);
  IntVect  nlarge(1);
  IntVect  boxlenlarge = ba[0].size();
  IntVect  boxlensmall(boxlenlarge);
  IntVect  stride(AMREX_D_DECL(1, nbox[0], nbox[0]*nbox[1]));
  // find nlarge in each dimension
  for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
    for (int i=1; i<nbox[idim]; ++i) {
      IntVect next_boxlen = ba[i*stride[idim]].size();
      if (next_boxlen[idim] < boxlenlarge[idim]) {
        boxlensmall[idim] = next_boxlen[idim];
        break;
      }
      boxlenlarge = next_boxlen;
      ++nlarge[idim];
    }
  }
  // offset for small box, i.e., total length of large boxes
  // IntVect offset_boxsmall(nlarge * boxlenlarge);
  
  // vector to hold # particles per box
  Gpu::DeviceVector<int> np_box(ba.size(), 0);
  auto np_box_ptr = np_box.dataPtr();

  // find # particles for each box
  const auto  dx_inv  = geom.InvCellSizeArray();
  const auto  prob_lo = geom.ProbLoArray();
  for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {
    const auto& aos     = pti.GetArrayOfStructs();
    const auto* pstruct = aos().dataPtr();
    const int   np      = pti.numRealParticles();
    ParallelFor(np, [pstruct, prob_lo, dx_inv, nbox, boxlenlarge, boxlensmall, 
                     nlarge, np_box_ptr]
      AMREX_GPU_DEVICE (int i) noexcept
      {
        const auto p = pstruct[i];
        IntVect cell_ijk(0), box_ijk(0);
        for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
          int offset_boxsmall = nlarge[idim] * boxlenlarge[idim];
          cell_ijk[idim] = static_cast<int>(Math::floor((p.pos(idim) - prob_lo[idim]) * dx_inv[idim]));
          box_ijk[idim]  = cell_ijk[idim] > offset_boxsmall ?
                           nlarge[idim] + (cell_ijk[idim] - offset_boxsmall) / boxlensmall[idim] :
                           cell_ijk[idim] / boxlenlarge[idim];
        }
        int box_idx = box_ijk[0] + box_ijk[1]*nbox[0] + box_ijk[2]*nbox[0]*nbox[1];
        Gpu::Atomic::Add(np_box_ptr + box_idx, 1);
      });// end parallel for
  }// end pariter for

  // collect # particle per box from all processes
  ParallelDescriptor::ReduceIntSum(np_box_ptr, np_box.size());

  // free current cost and construst new
  if (cost[lev] != nullptr)  
    delete cost[lev];
  cost[lev] = new MultiFab(ParticleBoxArray(lev), 
                           ParticleDistributionMap(lev), 1, 0);
  
  // set value by particle counts
  for (MFIter mfi(*cost[lev], false); mfi.isValid(); ++mfi) {
    Real wt = np_box_ptr[mfi.index()] / ba[mfi.index()].d_numPts();
    (*cost[lev])[mfi].setVal<RunOn::Device>(wt);
  }

}


Real MFIXParticleContainer::particleImbalance()
{
  // # paricles on this process
  int local_count = 0;
  for (MFIXParIter pti(*this, 0); pti.isValid(); ++pti)
      local_count += pti.numParticles();
  // count total # particles if not counted
  if (m_total_numparticle <= 0) {
    m_total_numparticle = local_count;
    ParallelDescriptor::ReduceIntSum(m_total_numparticle);
  }
  // max # particles per process
  ParallelDescriptor::ReduceIntMax(local_count,
                                   ParallelDescriptor::IOProcessorNumber());

  return ( static_cast<Real>(m_total_numparticle)
         / ParallelDescriptor::NProcs()
         / (static_cast<Real>(local_count) + 1e-10));
}
