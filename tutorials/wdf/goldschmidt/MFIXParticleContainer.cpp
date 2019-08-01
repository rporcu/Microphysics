#include <AMReX.H>
#include "AMReX_Particles.H"
#include "AMReX_RealVect.H"
#include <iostream>
#include <MFIXParticleContainer.H>
#include <AMReX_LoadBalanceKD.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB_F.H>

#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBSupport.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

#include <math.h>

#include "mfix_F.H"
#include "mfix_des_F.H"
#include "mfix_eb_F.H"
#include "mfix_util_F.H"
#include "mfix_des_K.H"
#include "bc_mod_F.H"
#include "constant_mod_F.H"
#include "MFIX_DEM_Parms.H"

using namespace amrex;
using namespace std;

Real MFIXParticleContainer::gravity[3] {0.0};

int  MFIXParticleContainer::domain_bc[6] {0};

MFIXParticleContainer::MFIXParticleContainer (AmrCore* amr_core)
    : NeighborParticleContainer<realData::count,intData::count>
      (amr_core->GetParGDB(), 1)
{
    ReadStaticParameters();

    this->SetVerbose(0);

    // turn off certain components for ghost particle communication
    setRealCommComp(4, false);
    setRealCommComp(5, false);
    setRealCommComp(6, false);
    setRealCommComp(7, false);
    setRealCommComp(14, false);
    setRealCommComp(15, false);
    setRealCommComp(16, false);

    setIntCommComp(0, false);
    setIntCommComp(1, false);
    setIntCommComp(3, false);

    nlev = amr_core->maxLevel() + 1;

    get_domain_bc(domain_bc);
}

void MFIXParticleContainer::AllocData ()
{
    reserveData();
    resizeData();

}

void MFIXParticleContainer::PrintParticleCounts() {

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

void MFIXParticleContainer:: printParticles()
{
    const int lev = 0;
    const auto& plevel = GetParticles(lev);

    for (const auto& kv : plevel)
    {
       const auto& particles = kv.second.GetArrayOfStructs();

       for (unsigned i = 0; i < particles.numParticles(); ++i)
       {
          std::cout << "Particle ID  = " << i << " " << std::endl;
          std::cout << "X            = " << particles[i].pos(0) << " " << std::endl;
          std::cout << "Y            = " << particles[i].pos(1) << " " << std::endl;
          std::cout << "Z            = " << particles[i].pos(2) << " " << std::endl;
          std::cout << "state        = " << particles[i].idata(intData::state) << " " << std::endl;
          std::cout << "phase        = " << particles[i].idata(intData::phase) << " " << std::endl;
          std::cout << "Real properties = " << std::endl;

          for (int j = 0; j < realData::count; j++)
            std::cout << "property " << j << "  = " << particles[i].rdata(j) << " " << std::endl;

          std::cout << std::endl;
       }
    }
}

void MFIXParticleContainer::ReadStaticParameters ()
{
    static bool initialized = false;

    get_gravity(gravity);

    if (!initialized)
        initialized = true;
}

void MFIXParticleContainer::EvolveParticles(int lev, int nstep, Real dt, Real time,
                                            EBFArrayBoxFactory * ebfactory,
                                            const MultiFab * ls_phi, const iMultiFab * ls_valid,
                                            const int ls_refinement,
                                            MultiFab * cost, std::string & knapsack_weight_type)
{
    BL_PROFILE_REGION_START("mfix_dem::EvolveParticles()");
    BL_PROFILE("mfix_dem::EvolveParticles()");

    Real eps = std::numeric_limits<Real>::epsilon();

    amrex::Print() << "Evolving particles on level: " << lev << " ... with fluid dt " << dt << std::endl;

    /****************************************************************************
     * DEBUG flag toggles:                                                      *
     *   -> Print number of collisions                                          *
     *   -> Print max (over substeps) particle velocity at each time step       *
     *   -> Print max particle-wall and particle-particle forces                *
     ***************************************************************************/

    // Debug level controls the detail of debug outut:
    //   -> debug_level = 0 : no debug output
    //   -> debug_level = 1 : debug output for every fluid step
    //   -> debug_level = 2 : debug output for every substep
    const int debug_level = 0;

    /****************************************************************************
     * Geometry                                                                 *
     ***************************************************************************/

    Box domain(Geom(lev).Domain());

    const Real* dx = Geom(lev).CellSize();

    /****************************************************************************
     * Init substeps                                                            *
     ***************************************************************************/

    int   nsubsteps;
    Real  subdt;
    des_init_time_loop( &time, &dt, &nsubsteps, &subdt );

    /****************************************************************************
     * Get particle EB geometric info
     ***************************************************************************/
    const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());

    /****************************************************************************
     * Init temporary storage:                                                  *
     *   -> particle-particle, and particle-wall forces                         *
     *   -> particle-particle, and particle-wall torques                        *
     ***************************************************************************/
    std::map<PairIndex, Gpu::ManagedDeviceVector<Real>> tow;
    std::map<PairIndex, Gpu::ManagedDeviceVector<Real>> fc, pfor, wfor;

    std::map<PairIndex, bool> tile_has_walls;
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const Box& bx = pti.tilebox();
        PairIndex index(pti.index(), pti.LocalTileIndex());
        tow[index]  = Gpu::ManagedDeviceVector<Real>();
        fc[index]   = Gpu::ManagedDeviceVector<Real>();
        pfor[index] = Gpu::ManagedDeviceVector<Real>();
        wfor[index] = Gpu::ManagedDeviceVector<Real>();

        // Only call the routine for wall collisions if we actually have walls
        BL_PROFILE_VAR("ls_has_walls", has_wall);
        bool has_wall = false;
        if ((ebfactory != NULL)
            && ((*flags)[pti].getType(amrex::grow(bx,1)) == FabType::singlevalued))
        {
            has_wall = true;
        }
        else
        {
            int int_has_wall = 0;
            Real tol = std::min(dx[0], std::min(dx[1], dx[2])) / 2;
            ls_has_walls(& int_has_wall, BL_TO_FORTRAN_3D((* ls_phi)[pti]), & tol);
            has_wall = (int_has_wall > 0);
        }

        tile_has_walls[index] = has_wall;

        BL_PROFILE_VAR_STOP(has_wall);
    }

    /****************************************************************************
     * Iterate over sub-steps                                                   *
     ***************************************************************************/

    int ncoll_total = 0;  // Counts total number of collisions
    loc_maxvel  = RealVect(0., 0., 0.);  // Tracks max (absolute) velocity
    loc_maxpfor = RealVect(0., 0., 0.);  // Tracks max particle-particle force
    loc_maxwfor = RealVect(0., 0., 0.);  // Tracks max particle-wall force
    int n = 0; // Counts sub-steps

    while (n < nsubsteps)
    {
        int ncoll = 0;  // Counts number of collisions (over sub-steps)

        // Redistribute particles ever so often BUT always update the neighbour
        // list (Note that this fills the neighbour list after every
        // redistribute operation)
        if (n % 25 == 0) {
            clearNeighbors();
            Redistribute(0, 0, 0, 1);
            fillNeighbors();
            // send in "false" for sort_neighbor_list option

            buildNeighborList(MFIXCheckPair(), false);
        } else {
            updateNeighbors();
        }

#ifdef _OPENMP
#pragma omp parallel reduction(+:ncoll) if (Gpu::notInLaunchRegion())
#endif
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            // Timer used for load-balancing
            Real wt = ParallelDescriptor::second();

            const Box& bx = pti.tilebox();
            PairIndex index(pti.index(), pti.LocalTileIndex());

            const int nrp = GetParticles(lev)[index].numRealParticles();
            RealType* particles  = pti.GetArrayOfStructs().data();

            auto& plev = GetParticles(lev);
            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            ParticleType* pstruct = aos().dataPtr();

            // Neighbor particles
#ifdef AMREX_USE_CUDA
            int size_ng = aos.numNeighborParticles();
#else
            int size_ng = neighbors[lev][index].size();
            int size_nl = neighbor_list[lev][index].size();
#endif

            // Number of particles including neighbor particles
            int ntot = nrp + size_ng;

            // Particle-particle (and particle-wall) forces and torques. We need
            // these to be zero every time we start a new batch (i.e tile and
            // substep) of particles.
            tow[index].clear();
            fc[index].clear();
            tow[index].resize(ntot*3,0.0);
            fc[index].resize(ntot*3,0.0);

            Real* fc_ptr = fc[index].dataPtr();
            Real* tow_ptr = tow[index].dataPtr();

            // For debugging: keep track of particle-particle (pfor) and
            // particle-wall (wfor) forces
            pfor[index].clear();
            wfor[index].clear();
            pfor[index].resize(3 * ntot, 0.0);
            wfor[index].resize(3 * ntot, 0.0);

            /********************************************************************
             * Particle-Wall collision forces (and torques)                     *
             *******************************************************************/

            if (tile_has_walls[index])
            {
                // Calculate forces and torques from particle-wall collisions
#ifndef AMREX_USE_CUDA
                BL_PROFILE_VAR("calc_wall_collisions()", calc_wall_collisions);
#endif
                auto& geom = this->Geom(lev);
                const auto dxi = geom.InvCellSizeArray();
                const auto plo = geom.ProbLoArray();
                const auto phiarr = ls_phi->array(pti);

                AMREX_FOR_1D ( nrp, i,
                {
                    ParticleType& p = pstruct[i];
                    Real rp = p.rdata(realData::radius);

                    Real ls_value = interp_level_set(p, ls_refinement, phiarr, plo, dxi);

                    Real overlap_n = rp - ls_value;

                    if (ls_value < rp)
                    {
                        Real normal[3];
                        level_set_normal(p, ls_refinement, &normal[0], phiarr, plo, dxi);

                        normal[0] *= -1;
                        normal[1] *= -1;
                        normal[2] *= -1;

                        Real v_rot[3];
                        v_rot[0] = ls_value * p.rdata(realData::omegax);
                        v_rot[1] = ls_value * p.rdata(realData::omegay);
                        v_rot[2] = ls_value * p.rdata(realData::omegaz);

                        Real vreltrans[3];
                        Real cprod[3];

                        cross_product(v_rot, normal, cprod);
                        vreltrans[0] = p.rdata(realData::velx) + cprod[0];
                        vreltrans[1] = p.rdata(realData::vely) + cprod[1];
                        vreltrans[2] = p.rdata(realData::velz) + cprod[2];

                        Real vreltrans_norm = dot_product(vreltrans, normal);

                        Real vrel_t[3];
                        vrel_t[0] = vreltrans[0] - vreltrans_norm*normal[0];
                        vrel_t[1] = vreltrans[1] - vreltrans_norm*normal[1];
                        vrel_t[2] = vreltrans[2] - vreltrans_norm*normal[2];

                        int phase = p.idata(intData::phase);

                        Real kn_des_w   = DEMParams::kn_w;
                        Real etan_des_w = DEMParams::etan_w[phase-1];

                        // NOTE - we don't use the tangential components right now,
                        // but we might in the future
                        // Real kt_des_w = DEMParams::kt_w;
                        // Real etat_des_w = DEMParams::etat_w[phase-1];

                        Real fn[3];
                        Real ft[3];
                        Real overlap_t[3];
                        Real mag_overlap_t;

                        // calculate the normal contact force
                        fn[0] = -(kn_des_w*overlap_n*normal[0] + etan_des_w*vreltrans_norm*normal[0]);
                        fn[1] = -(kn_des_w*overlap_n*normal[1] + etan_des_w*vreltrans_norm*normal[1]);
                        fn[2] = -(kn_des_w*overlap_n*normal[2] + etan_des_w*vreltrans_norm*normal[2]);

                        // calculate the tangential displacement
                        overlap_t[0] = subdt*vrel_t[0];
                        overlap_t[1] = subdt*vrel_t[1];
                        overlap_t[2] = subdt*vrel_t[2];

                        mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));

                        if (mag_overlap_t > 0.0) {
                            Real fnmd = DEMParams::mew * sqrt(dot_product(fn, fn));
                            Real tangent[3];
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

                        // each particle updates its force (no need for atomics)
                        fc_ptr[i         ] += fn[0] + ft[0];
                        fc_ptr[i + ntot  ] += fn[1] + ft[1];
                        fc_ptr[i + 2*ntot] += fn[2] + ft[2];

                        Real tow_force[3];

                        cross_product(normal, ft, tow_force);

                        tow_ptr[i         ] += ls_value*tow_force[0];
                        tow_ptr[i + ntot  ] += ls_value*tow_force[1];
                        tow_ptr[i + 2*ntot] += ls_value*tow_force[2];
                    }
                });

                // Debugging: copy data from the fc (all forces) vector to
                // the wfor (wall forces) vector.
                if (debug_level > 0) {
                    Gpu::Device::streamSynchronize();
                    for (int i = 0; i < wfor[index].size(); i++ ) {
                        wfor[index][i] = fc[index][i];
                    }
                }
#ifndef AMREX_USE_CUDA
                BL_PROFILE_VAR_STOP(calc_wall_collisions);
#endif
            }

            /********************************************************************
             * Particle-Particle collision forces (and torques)                 *
             *******************************************************************/

#ifndef AMREX_USE_CUDA
            BL_PROFILE_VAR("calc_particle_collisions()", calc_particle_collisions);
#endif

#ifdef AMREX_USE_CUDA
            auto nbor_data = m_neighbor_list[index].data();

            constexpr Real small_number = 1.0e-15;
            Gpu::DeviceScalar<int> ncoll_gpu(ncoll);
            int* pncoll = ncoll_gpu.dataPtr();

#if defined(AMREX_DEBUG) || defined(AMREX_USE_ASSERTION)
            Real eps = std::numeric_limits<Real>::epsilon();
#endif
            // now we loop over the neighbor list and compute the forces
            AMREX_FOR_1D ( nrp, i,
            {
                ParticleType& p1 = pstruct[i];
                
                for (const auto& p2 : nbor_data.getNeighbors(i))
                {
                    Real dx = p2.pos(0) - p1.pos(0);
                    Real dy = p2.pos(1) - p1.pos(1);
                    Real dz = p2.pos(2) - p1.pos(2);

                    Real r2 = dx*dx + dy*dy + dz*dz;
                    Real r_lm = p1.rdata(realData::radius) + p2.rdata(realData::radius);

                    if ( r2 <= (r_lm - small_number)*(r_lm - small_number) )
                    {
                        Gpu::Atomic::Add(pncoll, 1);
                        Real dist_mag = sqrt(r2);
                        AMREX_ASSERT(dist_mag >= eps);

                        Real normal[3];
                        normal[0] = dx / dist_mag;
                        normal[1] = dy / dist_mag;
                        normal[2] = dz / dist_mag;

                        Real overlap_n = r_lm - dist_mag;
                        Real vrel_trans_norm;
                        Real vrel_t[3];

                        cfrelvel(p1, p2, vrel_trans_norm, vrel_t, normal, dist_mag);

                        int phase1 = p1.idata(intData::phase);
                        int phase2 = p2.idata(intData::phase);

                        Real kn_des = DEMParams::kn;
                        Real etan_des = DEMParams::etan[phase1-1][phase2-1];

                        // NOTE - we don't use the tangential components right now,
                        // but we might in the future
                        // Real kt_des = DEMParams::kt;
                        // Real etat_des = DEMParams::etat[phase1-1][phase2-1];

                        Real fn[3];
                        Real ft[3];
                        Real overlap_t[3];
                        Real mag_overlap_t;

                        // calculate the normal contact force
                        fn[0] = -(kn_des*overlap_n*normal[0] + etan_des*vrel_trans_norm*normal[0]);
                        fn[1] = -(kn_des*overlap_n*normal[1] + etan_des*vrel_trans_norm*normal[1]);
                        fn[2] = -(kn_des*overlap_n*normal[2] + etan_des*vrel_trans_norm*normal[2]);

                        // calculate the tangential overlap
                        overlap_t[0] = subdt*vrel_t[0];
                        overlap_t[1] = subdt*vrel_t[1];
                        overlap_t[2] = subdt*vrel_t[2];
                        mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));

                        if (mag_overlap_t > 0.0) {
                            Real fnmd = DEMParams::mew * sqrt(dot_product(fn, fn));
                            Real tangent[3];
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

                        // each particle updates its force (no need for atomics)
                        fc_ptr[i         ] += fn[0] + ft[0];
                        fc_ptr[i + ntot  ] += fn[1] + ft[1];
                        fc_ptr[i + 2*ntot] += fn[2] + ft[2];

                        Real r1 = p1.rdata(realData::radius);
                        Real r2 = p2.rdata(realData::radius);

                        Real dist_cl = 0.5 * (dist_mag + (r1*r1 - r2*r2) / dist_mag);
                        dist_cl = dist_mag - dist_cl;

                        Real tow_force[3];

                        cross_product(normal, ft, tow_force);

                        tow_ptr[i         ] += dist_cl*tow_force[0];
                        tow_ptr[i + ntot  ] += dist_cl*tow_force[1];
                        tow_ptr[i + 2*ntot] += dist_cl*tow_force[2];
                        }
                    }
            });
#else
            calc_particle_collisions ( particles                          , &nrp,
                                       neighbors[lev][index].dataPtr()    , &size_ng,
                                       neighbor_list[lev][index].dataPtr(), &size_nl,
                                       tow[index].dataPtr(), fc[index].dataPtr(),
                                       &subdt, &ncoll);
#endif

            // Debugging: copy data from the fc (all forces) vector to the wfor
            // (wall forces) vector. Note that since fc already contains the
            // wall forces, these need to be subtracted here.
            if (debug_level > 0)
            {
                for (int i = 0; i < pfor[index].size(); i++ ) {
                    pfor[index][i] = fc[index][i] - wfor[index][i];
                }
            }

#ifndef AMREX_USE_CUDA
            BL_PROFILE_VAR_STOP(calc_particle_collisions);
#endif

#ifndef AMREX_USE_CUDA
            BL_PROFILE_VAR("des_time_march()", des_time_march);
#endif
            /********************************************************************
             * Move particles based on collision forces and torques             *
             *******************************************************************/

            GpuArray<Real, 3> grav;
            grav[0] = gravity[0];
            grav[1] = gravity[1];
            grav[2] = gravity[2];

            const auto p_lo = Geom(lev).ProbLoArray();
            const auto p_hi = Geom(lev).ProbHiArray();

            int x_lo_bc = domain_bc[0];
            int x_hi_bc = domain_bc[1];
            int y_lo_bc = domain_bc[2];
            int y_hi_bc = domain_bc[3];
            int z_lo_bc = domain_bc[4];
            int z_hi_bc = domain_bc[5];

            AMREX_FOR_1D ( nrp, i,
            {
                ParticleType& p = pstruct[i];

                p.rdata(realData::velx) += subdt * (
                    (p.rdata(realData::dragx) + fc_ptr[i       ]) /  p.rdata(realData::mass) + grav[0]);
                p.rdata(realData::vely) += subdt * (
                    (p.rdata(realData::dragy) + fc_ptr[i + ntot]) /  p.rdata(realData::mass) + grav[1]);
                p.rdata(realData::velz) += subdt * (
                    (p.rdata(realData::dragz) + fc_ptr[i+2*ntot]) /  p.rdata(realData::mass) + grav[2]);

                p.rdata(realData::omegax) += subdt * p.rdata(realData::oneOverI) * tow_ptr[i       ];
                p.rdata(realData::omegay) += subdt * p.rdata(realData::oneOverI) * tow_ptr[i+  ntot];
                p.rdata(realData::omegaz) += subdt * p.rdata(realData::oneOverI) * tow_ptr[i+2*ntot];

                p.pos(0) += subdt * p.rdata(realData::velx);
                p.pos(1) += subdt * p.rdata(realData::vely);
                p.pos(2) += subdt * p.rdata(realData::velz);

                if (x_lo_bc && p.pos(0) < p_lo[0])
                {
                    p.pos(0) = p_lo[0] + eps;
                    p.rdata(realData::velx) = -p.rdata(realData::velx);
                }
                if (x_hi_bc && p.pos(0) > p_hi[0])
                {
                   p.pos(0) = p_hi[0] - eps;
                   p.rdata(realData::velx) = -p.rdata(realData::velx);
                }
                if (y_lo_bc && p.pos(1) < p_lo[1])
                {
                    p.pos(1) = p_lo[1] + eps;
                    p.rdata(realData::vely) = -p.rdata(realData::vely);
                }
                if (y_hi_bc && p.pos(1) > p_hi[1])
                {
                   p.pos(1) = p_hi[1] - eps;
                   p.rdata(realData::vely) = -p.rdata(realData::vely);
                }
                if (z_lo_bc && p.pos(2) < p_lo[2])
                {
                   p.pos(2) = p_lo[2] + eps;
                   p.rdata(realData::velz) = -p.rdata(realData::velz);
                }
                if (z_hi_bc && p.pos(2) > p_hi[2])
                {
                   p.pos(2) = p_hi[2] - eps;
                   p.rdata(realData::velz) = -p.rdata(realData::velz);
                }
            });

            Gpu::Device::streamSynchronize();

#ifndef AMREX_USE_CUDA
            BL_PROFILE_VAR_STOP(des_time_march);
#endif

#ifdef AMREX_USE_CUDA
            ncoll = ncoll_gpu.dataValue();
#endif
            call_usr2_des(&nrp, pstruct);

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
                (*cost)[pti].plus(wt, tbx);
            }
        }

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
    clearNeighbors();
    Redistribute(0, 0, 0, 1);

    /****************************************************************************
     * DEBUG: output the total number of collisions over all substeps           *
     *        output the maximum velocity and forces over all substeps          *
     ***************************************************************************/
    if (debug_level > 0) {
        ParallelDescriptor::ReduceIntSum(ncoll_total, ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "Number of collisions: " << ncoll_total << " in " << nsubsteps << " substeps " << std::endl;
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int nrp   = NumberOfParticles(pti);
        void* particles = pti.GetArrayOfStructs().data();

        call_usr3_des( &nrp, particles );
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

void MFIXParticleContainer::
CalcVolumeFraction(const Vector<std::unique_ptr<MultiFab>> & mf_to_be_filled,
                   const Vector<std::unique_ptr<EBFArrayBoxFactory>> & ebfactory,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_ilo,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_ihi,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_jlo,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_jhi,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_klo,
                   const Vector<std::unique_ptr<IArrayBox>> & bc_khi,
                   const BcList & bc_list,
                   int nghost )
{
    // NOTE: ebfactory HAS to be the PARTICLE EB factory!
    int fortran_volume_comp = 5;
    PICDeposition(mf_to_be_filled, ebfactory,
                  bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo,bc_khi,
                  bc_list, fortran_volume_comp,nghost);

    for (int lev = 0; lev < nlev; lev++)
    {
        // Now define this mf = (1 - particle_vol)
        mf_to_be_filled[lev]->mult(-1.0,mf_to_be_filled[lev]->nGrow());
        mf_to_be_filled[lev]->plus( 1.0,mf_to_be_filled[lev]->nGrow());

        // We set ep_g to 1 rather than 0 in covered cells so that when we divide by ep_g
        //    following the projection we don't have to protect against divide by 0.
        EB_set_covered(*mf_to_be_filled[lev],1.0);

        // Impose a lower bound on volume fraction
        CapSolidsVolFrac(*mf_to_be_filled[lev]);
    }

    // HACK -- we really should average down (ep_g * volfrac) not ep_g.
    for (int lev = nlev - 1; lev > 0; lev --)
    {
        amrex::EB_average_down(* mf_to_be_filled[lev], * mf_to_be_filled[lev - 1],
                               0, 1, m_gdb->refRatio(lev - 1));
    }
}

void MFIXParticleContainer::
CalcDragOnFluid(const Vector<std::unique_ptr<MultiFab>> & drag_mf,
                const Vector<std::unique_ptr<EBFArrayBoxFactory>> & ebfactory,
                const Vector<std::unique_ptr<IArrayBox>> & bc_ilo,
                const Vector<std::unique_ptr<IArrayBox>> & bc_ihi,
                const Vector<std::unique_ptr<IArrayBox>> & bc_jlo,
                const Vector<std::unique_ptr<IArrayBox>> & bc_jhi,
                const Vector<std::unique_ptr<IArrayBox>> & bc_klo,
                const Vector<std::unique_ptr<IArrayBox>> & bc_khi,
                int nghost )
{
    int fortran_beta_comp = 15;
    int fortran_vel_comp  =  9;
    PICMultiDeposition(drag_mf, ebfactory,
                       bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo, bc_khi,
                       fortran_beta_comp, fortran_vel_comp, nghost);
}

void MFIXParticleContainer::
PICDeposition(const amrex::Vector< std::unique_ptr<MultiFab> >& mf_to_be_filled,
              const amrex::Vector< std::unique_ptr<EBFArrayBoxFactory>  >& ebfactory,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
              const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
              const BcList & bc_list,
              int fortran_particle_comp, int nghost )
{
    BL_PROFILE("MFIXParticleContainer::PICDeposition()");

    int ncomp = 1;

    MultiFab* mf_pointer[nlev];

    // Start the timers ...
    const Real strttime = ParallelDescriptor::second();

    if (nlev > 2)
      amrex::Abort("For right now MFIXParticleContainer::PICDeposition can only handle up to 2 levels");

    for (int lev = 0; lev < nlev; lev++)
    {
       if (lev == 0 && OnSameGrids(lev, *mf_to_be_filled[lev])) {
          // If we are already working with the internal mf defined on the
          // particle_box_array, then we just work with this.
          mf_pointer[lev] = mf_to_be_filled[lev].get();

       } else if (lev == 0 && !OnSameGrids(lev, *mf_to_be_filled[lev]))  {
          // If mf_to_be_filled is not defined on the particle_box_array, then we need
          // to make a temporary here and copy into mf_to_be_filled at the end.
          mf_pointer[lev] = new MultiFab(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                                         ncomp, mf_to_be_filled[lev]->nGrow());

       } else {
          // If lev > 0 we make a temporary at the coarse resolution
          BoxArray ba_crse(amrex::coarsen(ParticleBoxArray(lev),this->m_gdb->refRatio(0)));
          mf_pointer[lev] = new MultiFab(ba_crse, ParticleDistributionMap(lev),
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

    // We always use the coarse dx
    const Geometry& gm          = Geom(0);
    const auto      plo         = gm.ProbLoArray();
    const auto      dx          = gm.CellSizeArray();
    const auto      dxi         = gm.InvCellSizeArray();

    using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;

    const FabArray<EBCellFlagFab>* flags;

    for (int lev = 0; lev < nlev; lev++)
    {
        Vector<int> ngrow = {1,1,1};
        std::unique_ptr<EBFArrayBoxFactory> crse_factory;

        if (lev == 0) {
            flags   = &(ebfactory[lev]->getMultiEBCellFlagFab());
        } else {
            // std::unique_ptr<EBFArrayBoxFactory>
            crse_factory = makeEBFabFactory(
                  gm, mf_pointer[lev]->boxArray(), mf_pointer[lev]->DistributionMap(),
                  ngrow, EBSupport::volume);
            flags   = &(crse_factory->getMultiEBCellFlagFab());
        }

        const MultiFab* volfrac = (lev == 0) ? &(ebfactory[lev]->getVolFrac()) : &(crse_factory->getVolFrac());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       {
        FArrayBox local_vol;
        for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {
            const auto& particles = pti.GetArrayOfStructs();
            const ParticleType* pstruct = particles().dataPtr();

            const long nrp = pti.numParticles();
            FArrayBox& fab = (*mf_pointer[lev])[pti];
#ifdef _OPENMP
            Box tile_box = pti.tilebox();
            tile_box.grow(1);
            local_vol.resize(tile_box,ncomp);
            local_vol = 0.0;
#else
            const Box& box = fab.box();
#endif

            const Box& bx  = pti.tilebox(); // I need a box without ghosts

            if ((*flags)[pti].getType(bx) != FabType::covered )
            {
                auto volarr = fab.array();
                auto flagsarr = (*flags)[pti].array();
                auto vratioarr = (*volfrac)[pti].array();

                AMREX_FOR_1D ( nrp, ip,
                {
                    const ParticleType& p = pstruct[ip];

                    amrex::Real reg_cell_vol = dx[0]*dx[1]*dx[2];

                    amrex::Real x = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
                    amrex::Real y = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
                    amrex::Real z = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

                    int i = std::floor(x);
                    int j = std::floor(y);
                    int k = std::floor(z);

                    amrex::Real wx_hi = x - i;
                    amrex::Real wy_hi = y - j;
                    amrex::Real wz_hi = z - k;

                    amrex::Real wx_lo = 1.0 - wx_hi;
                    amrex::Real wy_lo = 1.0 - wy_hi;
                    amrex::Real wz_lo = 1.0 - wz_hi;

                    amrex::Real weights[2][2][2];

                    weights[0][0][0] = vratioarr(i-1, j-1, k-1) * wx_lo * wy_lo * wz_lo;
                    weights[0][0][1] = vratioarr(i-1, j-1, k  ) * wx_lo * wy_lo * wz_hi;
                    weights[0][1][0] = vratioarr(i-1, j  , k-1) * wx_lo * wy_hi * wz_lo;
                    weights[0][1][1] = vratioarr(i-1, j  , k  ) * wx_lo * wy_hi * wz_hi;
                    weights[1][0][0] = vratioarr(i  , j-1, k-1) * wx_hi * wy_lo * wz_lo;
                    weights[1][0][1] = vratioarr(i  , j-1, k  ) * wx_hi * wy_lo * wz_hi;
                    weights[1][1][0] = vratioarr(i  , j  , k-1) * wx_hi * wy_hi * wz_lo;
                    weights[1][1][1] = vratioarr(i  , j  , k  ) * wx_hi * wy_hi * wz_hi;

                    amrex::Real total_weight = 0.0;
                    for (int ii = 0; ii <= 1; ++ii)
                        for (int jj = 0; jj <= 1; ++jj)
                            for (int kk = 0; kk <= 1; ++kk)
                                total_weight += weights[ii][jj][kk];

                    for (int ii = 0; ii <= 1; ++ii)
                        for (int jj = 0; jj <= 1; ++jj)
                            for (int kk = 0; kk <= 1; ++kk)
                                weights[ii][jj][kk] /= total_weight;

                    amrex::Real pvol = p.rdata(realData::volume);

                    for (int ii = -1; ii <= 0; ++ii) {
                        for (int jj = -1; jj <= 0; ++jj) {
                            for (int kk = -1; kk <= 0; ++kk) {
                                if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                                    continue;
                                amrex::Real this_cell_vol = vratioarr(i+ii,j+jj,k+kk) * reg_cell_vol;
                                amrex::Gpu::Atomic::Add(&volarr(i+ii,j+jj,k+kk),
                                                        weights[ii+1][jj+1][kk+1]*pvol/this_cell_vol);
                            }
                        }
                    }
                });
            }

#ifdef _OPENMP
            fab.atomicAdd(local_vol, tile_box, tile_box, 0, 0, ncomp);
#endif
        }
       }

       // Move any field deposited outside the domain back into the domain
       // when BC is pressure inlet and mass inflow.
       Box domain(Geom(lev).Domain());

       const int minf = bc_list.get_minf();
       const int pinf = bc_list.get_pinf();
       const int pout = bc_list.get_pout();

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
           AMREX_HOST_DEVICE_FOR_3D(sbx_yz, i, j, k,
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
           AMREX_HOST_DEVICE_FOR_3D(sbx_yz, i, j, k,
           {
             if(bc_ihi_type(dom_hi.x+1,j,k,0) == pinf or
                bc_ihi_type(dom_hi.x+1,j,k,0) == minf)
             {
               vol(ihi,j,k) += vol(ihi+1,j,k);
               vol(ihi+1,j,k) = 0;
             }
           });
         }

#ifdef AMREX_USE_CUDA
         Gpu::Device::synchronize();
#endif

         if(sbx_lo[1] < dom_lo.y)
         {
           const int jlo = dom_lo.y;
           AMREX_HOST_DEVICE_FOR_3D(sbx_xz, i, j, k,
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
           AMREX_HOST_DEVICE_FOR_3D(sbx_xz, i, j, k,
           {
             if(bc_jhi_type(i,dom_hi.y+1,k,0) == pinf or
                bc_jhi_type(i,dom_hi.y+1,k,0) == minf)
             {
               vol(i,jhi,k) += vol(i,jhi+1,k);
               vol(i,jhi+1,k) = 0;
             }
           });
         }

#ifdef AMREX_USE_CUDA
         Gpu::Device::synchronize();
#endif

         if(sbx_lo[2] < dom_lo.z)
         {
           const int klo = dom_lo.z;
           AMREX_HOST_DEVICE_FOR_3D(sbx_xy, i, j, k,
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
           AMREX_HOST_DEVICE_FOR_3D(sbx_xy, i, j, k,
           {
             if(bc_khi_type(i,j,dom_hi.z+1,0) == pinf or
                bc_khi_type(i,j,dom_hi.z+1,0) == minf)
             {
               vol(i,j,khi) += vol(i,j,khi+1);
               vol(i,j,khi+1) = 0;
             }
           });
         }
       }
    }

#ifdef AMREX_USE_CUDA
    Gpu::Device::synchronize();
#endif

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
            mf_to_be_filled[lev]->setVal(0.0);
            amrex::InterpFromCoarseLevel(*mf_to_be_filled[lev], time, *mf_pointer[lev-1],
                                         0, 0, 1, Geom(lev-1), Geom(lev),
                                         cphysbc, 0, fphysbc, 0,
                                         ref_ratio, mapper,
                                         bcs, 0);
        }
    }

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.

    if (mf_pointer[0] != mf_to_be_filled[0].get())
       mf_to_be_filled[0]->copy(*mf_pointer[0],0,0,ncomp);

    for (int lev = 0; lev < nlev; lev++)
       if (mf_pointer[lev] != mf_to_be_filled[lev].get())
          delete mf_pointer[lev];

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICDeposition time: " << stoptime << '\n';
    }
}

void MFIXParticleContainer::
PICMultiDeposition(const amrex::Vector< std::unique_ptr<MultiFab> >& drag_mf,
                   const amrex::Vector< std::unique_ptr<EBFArrayBoxFactory>  >& ebfactory,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ilo,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_ihi,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jlo,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_jhi,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_klo,
                   const amrex::Vector< std::unique_ptr<amrex::IArrayBox> >& bc_khi,
                   int fortran_beta_comp, int fortran_vel_comp, int nghost )
{
    BL_PROFILE("MFIXParticleContainer::PICMultiDeposition()");

    const Real strttime = ParallelDescriptor::second();

    if (nlev > 2)
        amrex::Abort("For right now MFIXParticleContainer::PICMultiDeposition can only handle up to 2 levels");

    MultiFab*  drag_ptr[nlev];

    for (int lev = 0; lev < nlev; lev++)
    {
       if (lev == 0 && OnSameGrids(lev, *drag_mf[lev])) {
          // If we are already working with the internal mf defined on the
          // particle_box_array, then we just work with this.
          drag_ptr[lev] = drag_mf[lev].get();

       } else if (lev == 0 && !OnSameGrids(lev, *drag_mf[lev]))  {
          // If beta_mf is not defined on the particle_box_array, then we need
          // to make a temporary here and copy into beta_mf at the end.
          drag_ptr[lev] = new MultiFab(ParticleBoxArray(lev), ParticleDistributionMap(lev),
                                       drag_mf[lev]->nComp(), drag_mf[lev]->nGrow());

       } else {
          // If lev > 0 we make a temporary at the coarse resolution
          BoxArray ba_crse(amrex::coarsen(ParticleBoxArray(lev),this->m_gdb->refRatio(0)));
          drag_ptr[lev] = new MultiFab(ba_crse, ParticleDistributionMap(lev),drag_mf[lev]->nComp(),1);
       }

       // We must have ghost cells for each FAB so that a particle in one grid can spread
       // its effect to an adjacent grid by first putting the value into ghost cells of its
       // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
       // to another grid's valid region.
       if (drag_ptr[lev]->nGrow() < 1)
          amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

       drag_ptr[lev]->setVal(0.0,0,4,drag_ptr[lev]->nGrow());
    }

    // We always use the coarse dx
    const Geometry& gm          = Geom(0);
    const auto      plo         = gm.ProbLoArray();
    const auto      dx          = gm.CellSizeArray();
    const auto      dxi         = gm.InvCellSizeArray();

    using ParConstIter = ParConstIter<realData::count,intData::count,0,0>;
    const FabArray<EBCellFlagFab>* flags;

    for (int lev = 0; lev < nlev; lev++)
    {
        Vector<int> ngrow = {1,1,1};
        std::unique_ptr<EBFArrayBoxFactory> crse_factory;

        if (lev == 0) {
           flags   = &(ebfactory[lev]->getMultiEBCellFlagFab());
        } else {
           // std::unique_ptr<EBFArrayBoxFactory>
           crse_factory = makeEBFabFactory(
                  gm, drag_ptr[lev]->boxArray(), drag_ptr[lev]->DistributionMap(),
                  ngrow, EBSupport::volume);
           flags   = &(crse_factory->getMultiEBCellFlagFab());
        }

        const MultiFab* volfrac = (lev == 0) ? &(ebfactory[lev]->getVolFrac()) : &(crse_factory->getVolFrac());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {

        const int* lo;
        const int* hi;

        FArrayBox local_vol;
         for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {

            const auto& particles = pti.GetArrayOfStructs();
            const ParticleType* pstruct = particles().dataPtr();
            const long nrp = pti.numParticles();

            FArrayBox& drag_fab = (*drag_ptr[lev])[pti];

#ifdef _OPENMP
            // Note that we actually grow the tilebox rather than calling growntilebox
            //     because we need the overlap even in the interior.
            Box grown_tilebox = pti.tilebox();
            grown_tilebox.grow(1);

            int ncomp = 4;
            local_vol.resize(grown_tilebox,ncomp);

            local_vol = 0.0;

            lo = grown_tilebox.loVect();
            hi = grown_tilebox.hiVect();
#else

            const Box& bx  = drag_fab.box();

            lo = bx.loVect();
            hi = bx.hiVect();
#endif
            const Box& box = pti.tilebox(); // I need a box without ghosts

            if ((*flags)[pti].getType(box) != FabType::covered )
            {
                auto drag_arr = drag_fab.array();
                auto flagsarr = (*flags)[pti].array();
                auto vratioarr = (*volfrac)[pti].array();

                AMREX_FOR_1D ( nrp, ip,
                {
                    const ParticleType& p = pstruct[ip];

                    amrex::Real reg_cell_vol = dx[0]*dx[1]*dx[2];

                    amrex::Real x = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
                    amrex::Real y = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
                    amrex::Real z = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

                    int i = std::floor(x);
                    int j = std::floor(y);
                    int k = std::floor(z);

                    amrex::Real wx_hi = x - i;
                    amrex::Real wy_hi = y - j;
                    amrex::Real wz_hi = z - k;

                    amrex::Real wx_lo = 1.0 - wx_hi;
                    amrex::Real wy_lo = 1.0 - wy_hi;
                    amrex::Real wz_lo = 1.0 - wz_hi;

                    amrex::Real weights[2][2][2];

                    weights[0][0][0] = vratioarr(i-1, j-1, k-1) * wx_lo * wy_lo * wz_lo;
                    weights[0][0][1] = vratioarr(i-1, j-1, k  ) * wx_lo * wy_lo * wz_hi;
                    weights[0][1][0] = vratioarr(i-1, j  , k-1) * wx_lo * wy_hi * wz_lo;
                    weights[0][1][1] = vratioarr(i-1, j  , k  ) * wx_lo * wy_hi * wz_hi;
                    weights[1][0][0] = vratioarr(i  , j-1, k-1) * wx_hi * wy_lo * wz_lo;
                    weights[1][0][1] = vratioarr(i  , j-1, k  ) * wx_hi * wy_lo * wz_hi;
                    weights[1][1][0] = vratioarr(i  , j  , k-1) * wx_hi * wy_hi * wz_lo;
                    weights[1][1][1] = vratioarr(i  , j  , k  ) * wx_hi * wy_hi * wz_hi;

                    amrex::Real total_weight = 0.0;
                    for (int ii = 0; ii <= 1; ++ii)
                        for (int jj = 0; jj <= 1; ++jj)
                            for (int kk = 0; kk <= 1; ++kk)
                                total_weight += weights[ii][jj][kk];

                    for (int ii = 0; ii <= 1; ++ii)
                        for (int jj = 0; jj <= 1; ++jj)
                            for (int kk = 0; kk <= 1; ++kk)
                                weights[ii][jj][kk] /= total_weight;

                    amrex::Real pbeta = p.rdata(realData::dragx);
                    amrex::Real pvx   = p.rdata(realData::velx) * pbeta;
                    amrex::Real pvy   = p.rdata(realData::vely) * pbeta;
                    amrex::Real pvz   = p.rdata(realData::velz) * pbeta;

                    for (int ii = -1; ii <= 0; ++ii) {
                        for (int jj = -1; jj <= 0; ++jj) {
                            for (int kk = -1; kk <= 0; ++kk) {
                                if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                                    continue;

                                amrex::Real this_cell_vol = vratioarr(i+ii,j+jj,k+kk) * reg_cell_vol;
                                amrex::Real weight_vol = weights[ii+1][jj+1][kk+1]/this_cell_vol;

                                amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,0),weight_vol*pvx);
                                amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,1),weight_vol*pvy);
                                amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,2),weight_vol*pvz);
                                amrex::Gpu::Atomic::Add(&drag_arr(i+ii,j+jj,k+kk,3),weight_vol*pbeta);
                            }
                        }
                    }
                });
            }

#ifdef _OPENMP
            drag_fab.atomicAdd(local_vol, grown_tilebox, grown_tilebox, 0, 0, drag_fab.nComp());
#endif

         }
        }
    }

    int  src_nghost = 1;
    int dest_nghost = 0;
    for (int lev = 1; lev < nlev; lev++)
    {
        drag_ptr[0]->copy(*drag_ptr[lev],0,0,drag_ptr[0]->nComp(),src_nghost,dest_nghost,gm.periodicity(),FabArrayBase::ADD);
    }

    drag_ptr[0]->SumBoundary(gm.periodicity());

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
            drag_mf[lev]->setVal(0.0);
            amrex::InterpFromCoarseLevel(*drag_mf[lev], time, *drag_ptr[lev-1],
                                         0, 0, 1, Geom(lev-1), Geom(lev),
                                         cphysbc, 0, fphysbc, 0,
                                         ref_ratio, mapper,
                                         bcs, 0);
        }
    }

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.

    if (drag_ptr[0] != drag_mf[0].get())
    {
       drag_mf[0]->copy(*drag_ptr[0],0,0,drag_mf[0]->nComp());
    }

    for (int lev = 0; lev < nlev; lev++)
    {
       if (drag_ptr[lev] != drag_mf[lev].get())
          delete drag_ptr[lev];
    }

    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      amrex::Print() << "MFIXParticleContainer::PICMultiDeposition time: " << stoptime << '\n';
    }
}


void MFIXParticleContainer::writeAllAtLevel(int lev)
{
    // Not threaded because its print to terminal
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        int np = pti.numParticles();
        Gpu::HostVector<ParticleType> host_particles(np);
        Cuda::thrust_copy(particles.begin(), particles.end(), host_particles.begin());

        for (const auto& p: host_particles)
        {
           const IntVect& iv = Index(p, lev);

           RealVect xyz(p.pos(0), p.pos(1), p.pos(2));
           cout << " id " << p.id()
                << " index " << iv
                << " position " << xyz << endl;
       }
    }
}

void
MFIXParticleContainer::WriteAsciiFileForInit (const std::string& filename)
{
    BL_ASSERT(!filename.empty());

    int lev = 0;
    long nparticles = NumberOfParticlesAtLevel (lev);

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Have I/O processor open file and write out particle metadata.
        //
        std::ofstream File;

        File.open(filename.c_str(), std::ios::out|std::ios::trunc);

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

              Gpu::HostVector<ParticleType> host_particles(np);
              Cuda::thrust_copy(particles.begin(), particles.end(), host_particles.begin());

              int index = 0;
              for (const auto& p: host_particles)
              {
                  if (p.id() > 0) {
                      File << p.idata(intData::phase) << ' ';
                      File << p.pos(0) << ' ';
                      File << p.pos(1) << ' ';
                      File << p.pos(2) << ' ';
                      File << p.rdata(realData::radius) << ' ';
                      File << p.rdata(realData::density) << ' ';
                      File << p.rdata(realData::velx) << ' ';
                      File << p.rdata(realData::vely) << ' ';
                      File << p.rdata(realData::velz) << ' ';

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

void MFIXParticleContainer::GetParticleAvgProp(Real (&min_dp)[10], Real (&min_ro)[10],
                                               Real (&max_dp)[10], Real (&max_ro)[10],
                                               Real (&avg_dp)[10], Real (&avg_ro)[10])
{
   // The number of phases was previously hard set at 10, however lowering
   //  this number would make this code faster.
   int num_of_phases_in_use = 10; //Number of different phases being simulated

   // Cycle through the different phases, starting from 1
   for (int phse = 1; phse <= num_of_phases_in_use; ++phse)
   {
     Real p_num  = 0.0; //number of particle
     Real p_diam = 0.0; //particle diameters
     Real p_dens = 0.0; //particle density

     Real min_diam =  1.0e32;
     Real min_den  =  1.0e32;

     Real max_diam = -1.0e32;
     Real max_den  = -1.0e32;

     for (int lev = 0; lev < nlev; lev++)
     {
#ifdef _OPENMP
#pragma omp parallel reduction(+:p_num, p_diam, p_dens) if (Gpu::notInLaunchRegion())
#endif
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& particles = pti.GetArrayOfStructs();

            Gpu::HostVector<ParticleType> host_particles(pti.numParticles());
            Cuda::thrust_copy(particles.begin(), particles.end(), host_particles.begin());

            for (const auto& p: host_particles){
                if ( phse==p.idata(intData::phase) )
                {
                    p_num  += 1.0;
                    p_diam += p.rdata(realData::radius) * 2.0;
                    p_dens += p.rdata(realData::density);

                    min_diam = amrex::min(min_diam, p.rdata(realData::radius) * 2.0 );
                    min_den  = amrex::min(min_den,  p.rdata(realData::density) );

                    max_diam = amrex::max(max_diam, p.rdata(realData::radius) * 2.0 );
                    max_den  = amrex::max(max_den,  p.rdata(realData::density) );

                }
            }
        }
     }

     // A single MPI call passes all three variables
     ParallelDescriptor::ReduceRealSum({p_num,p_diam,p_dens});
     ParallelDescriptor::ReduceRealMin({min_diam, min_den});
     ParallelDescriptor::ReduceRealMax({max_diam, max_den});

     //calculate averages or set = zero if no particles of that phase
     if (p_num==0){
       avg_dp[phse-1] = 0.0;
       avg_ro[phse-1] = 0.0;

       min_dp[phse-1] = 0.0;
       min_ro[phse-1] = 0.0;

       max_dp[phse-1] = 0.0;
       max_ro[phse-1] = 0.0;

     } else {
       avg_dp[phse-1] = p_diam/p_num;
       avg_ro[phse-1] = p_dens/p_num;

       min_dp[phse-1] = min_diam;
       min_ro[phse-1] = min_den;

       max_dp[phse-1] = max_diam;
       max_ro[phse-1] = max_den;
     }
   }
}

void MFIXParticleContainer::UpdateMaxVelocity ()
{
    Real max_vel_x = loc_maxvel[0], max_vel_y = loc_maxvel[1], max_vel_z = loc_maxvel[2];

#ifdef _OPENMP
#pragma omp parallel reduction(max:max_vel_x,max_vel_y,max_vel_z) if (Gpu::notInLaunchRegion())
#endif
    for (int lev = 0; lev < nlev; lev++)
    {
       for(MFIXParIter pti(* this, lev); pti.isValid(); ++ pti)
       {
           auto & particles = pti.GetArrayOfStructs();
           int np = pti.numParticles();
           Gpu::HostVector<ParticleType> host_particles(np);
           Cuda::thrust_copy(particles.begin(), particles.end(), host_particles.begin());

           for(const auto & particle : host_particles)
           {
              max_vel_x = std::max(Real(std::fabs(particle.rdata(realData::velx))), max_vel_x);
              max_vel_y = std::max(Real(std::fabs(particle.rdata(realData::vely))), max_vel_y);
              max_vel_z = std::max(Real(std::fabs(particle.rdata(realData::velz))), max_vel_z);
           }
       }
    }
    loc_maxvel = RealVect(max_vel_x, max_vel_y, max_vel_z);
}

void MFIXParticleContainer::UpdateMaxForces( std::map<PairIndex, Gpu::ManagedDeviceVector<Real>> pfor,
                                             std::map<PairIndex, Gpu::ManagedDeviceVector<Real>> wfor)
{
    Real max_pfor_x = loc_maxpfor[0], max_pfor_y = loc_maxpfor[1], max_pfor_z = loc_maxpfor[2];
    Real max_wfor_x = loc_maxwfor[0], max_wfor_y = loc_maxwfor[1], max_wfor_z = loc_maxwfor[2];

    for (int lev = 0; lev < nlev; lev++)
    {
#ifdef _OPENMP
#pragma omp parallel reduction(max:max_pfor_x,max_pfor_y,max_pfor_z,max_wfor_x,max_wfor_y,max_wfor_z) if (Gpu::notInLaunchRegion())
#endif
        for(MFIXParIter pti(* this, lev); pti.isValid(); ++ pti)
        {
            PairIndex index(pti.index(), pti.LocalTileIndex());

            // Note the particle force data layout:
            //      p1_x, p2_x, ..., pn_x, p1_y, p2_y, ..., pn_y, p1_z, p2_z, ..., pn_z
            // Where n is the total number of particle and neighbor particles.
            const int nrp     = NumberOfParticles(pti);
#ifdef AMREX_USE_CUDA
            auto& plev = GetParticles(lev);
            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            int size_ng = aos.numNeighborParticles();
#else
            int size_ng = neighbors[lev][index].size();
#endif
            // Number of particles including neighbor particles
            const int ntot = nrp + size_ng;

            // Find max (abs) of particle-particle forces:
            for(int i = 0; i < ntot; i++ )
                max_pfor_x = std::max(Real(std::fabs(pfor[index][i])), max_pfor_x);
            for(int i = ntot; i < 2 * ntot; i++ )
                max_pfor_y = std::max(Real(std::fabs(pfor[index][i])), max_pfor_y);
            for(int i = 2 * ntot; i < 3 * ntot; i++ )
                max_pfor_z = std::max(Real(std::fabs(pfor[index][i])), max_pfor_z);

            // Find max (abs) of particle-wall forces:
            for(int i = 0; i < ntot; i++ )
                max_wfor_x = std::max(Real(std::fabs(wfor[index][i])), max_wfor_x);
            for(int i = ntot; i < 2 * ntot; i++ )
                max_wfor_y = std::max(Real(std::fabs(wfor[index][i])), max_wfor_y);
            for(int i = 2 * ntot; i < 3 * ntot; i++ )
                max_wfor_z = std::max(Real(std::fabs(wfor[index][i])), max_wfor_z);
        }
    }

    loc_maxpfor = RealVect(max_pfor_x, max_pfor_y, max_pfor_z);
    loc_maxwfor = RealVect(max_wfor_x, max_wfor_y, max_wfor_z);
}

RealVect MFIXParticleContainer::GetMaxVelocity()
{
    Real max_vel_x = loc_maxvel[0], max_vel_y = loc_maxvel[1], max_vel_z = loc_maxvel[2];

    ParallelDescriptor::ReduceRealMax({max_vel_x, max_vel_y, max_vel_z},
                                      ParallelDescriptor::IOProcessorNumber());

    RealVect max_vel(max_vel_x, max_vel_y, max_vel_z);

    return max_vel;
};

Vector<RealVect> MFIXParticleContainer::GetMaxForces()
{
    Real max_pfor_x = loc_maxpfor[0], max_pfor_y = loc_maxpfor[1], max_pfor_z = loc_maxpfor[2];
    Real max_wfor_x = loc_maxwfor[0], max_wfor_y = loc_maxwfor[1], max_wfor_z = loc_maxwfor[2];


    ParallelDescriptor::ReduceRealMax({ max_pfor_x, max_pfor_y, max_pfor_z,
                                        max_wfor_x, max_wfor_y, max_wfor_z      },
                                      ParallelDescriptor::IOProcessorNumber());

    Vector<RealVect> max_forces(2);
    max_forces[0] = RealVect(max_pfor_x, max_pfor_y, max_pfor_z);
    max_forces[1] = RealVect(max_wfor_x, max_wfor_y, max_wfor_z);

    return max_forces;
}

void
MFIXParticleContainer::BalanceParticleLoad_KDTree()
{
  int lev = 0;
  bool verbose = true;
  BoxArray old_ba = ParticleBoxArray(lev);

  if (NumberOfParticlesAtLevel(lev) == 0)
  {
     amrex::Print() << "No particles so can't use KDTree approach " << std::endl;
     return;
  }

  if (verbose)
  {
     Vector<long> num_part;
     num_part = NumberOfParticlesInGrid(0);
     long min_number = num_part[0];
     long max_number = num_part[0];
     for (int i = 0; i < old_ba.size(); i++)
     {
        max_number = std::max(max_number, num_part[i]);
        min_number = std::min(min_number, num_part[i]);
     }
     amrex::Print() << "Before KDTree: BA had " << old_ba.size() << " GRIDS " << std::endl;
     amrex::Print() << "Before KDTree: MIN/MAX NUMBER OF PARTICLES PER GRID  " <<
                        min_number << " " << max_number << std::endl;
  }

  Vector<Real> box_costs;

  BoxArray new_ba;
  Real cell_weight = 0.;
  loadBalanceKD::balance<MFIXParticleContainer>(*this, new_ba, ParallelDescriptor::NProcs(), cell_weight, box_costs);

  // Create a new DM to go with the new BA
  DistributionMapping new_dm = DistributionMapping::makeKnapSack(box_costs);

  Regrid(new_dm, new_ba);

  if (verbose)
  {
     Vector<long> num_part;
     num_part = NumberOfParticlesInGrid(0);
     long min_number = num_part[0];
     long max_number = num_part[0];
     for (int i = 0; i < new_ba.size(); i++)
     {
        max_number = std::max(max_number, num_part[i]);
        min_number = std::min(min_number, num_part[i]);
     }
     amrex::Print() << "After  KDTree: BA had " << new_ba.size() << " GRIDS " << std::endl;
     amrex::Print() << "After  KDTree: MIN/MAX NUMBER OF PARTICLES PER GRID  " <<
                        min_number << " " << max_number << std::endl;
  }
}

void MFIXParticleContainer::
ComputeAverageVelocities ( const int lev,
                           const amrex::Real time,
                           const string&  basename,
                           const Vector<int>& avg_vel_p,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_x_w,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_x_e,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_y_s,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_y_n,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_z_b,
                           const Gpu::ManagedDeviceVector<Real>& avg_region_z_t )
{

  // Count number of calls -- Used to determin when to create file from scratch
  static int ncalls = 0;
  ++ncalls;

  int  nregions = avg_region_x_w.size();

  if(avg_vel_p.size() > 0)
    {

      //
      // Check the regions are defined correctly
      //
      if (  ( avg_region_x_e.size() != nregions ) ||
            ( avg_region_y_s.size() != nregions ) ||
            ( avg_region_y_n.size() != nregions ) ||
            ( avg_region_z_b.size() != nregions ) ||
            ( avg_region_z_t.size() != nregions )  )
        {
          amrex::Print () << "ComputeAverageVelocities: some regions are not properly defined: skipping.";
          return;
        }

      vector<long> region_np (nregions, 0);
      vector<Real> region_velx (nregions, 0.0);
      vector<Real> region_vely (nregions, 0.0);
      vector<Real> region_velz (nregions, 0.0);

      for ( int nr = 0; nr < nregions; ++nr )
        {

          amrex::Print() << "size of avg_vel_p " << avg_vel_p[nr] << "\n";

          // This region isn't needed for particle data.
          if( avg_vel_p[nr] == 0) continue;

          // Create Real box for this region
          RealBox avg_region ( {AMREX_D_DECL(avg_region_x_w[nr],avg_region_y_s[nr],avg_region_z_b[nr])},
                               {AMREX_D_DECL(avg_region_x_e[nr],avg_region_y_n[nr],avg_region_z_t[nr])} );

          // Jump to next iteration if this averaging region is not valid
          if ( !avg_region.ok () )
            {
              amrex::Print() << "ComputeAverageVelocities: region "<< nr <<" is invalid: skipping\n";
              continue;
            }

          long sum_np     = 0;    // Number of particle in avg region
          Real sum_velx   = 0.;
          Real sum_vely   = 0.;
          Real sum_velz   = 0.;

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum_np,sum_velx,sum_vely,sum_velz) if (Gpu::notInLaunchRegion())
#endif
          for ( MFIXParIter pti(*this, lev); pti.isValid(); ++ pti)
            {
              Box bx       = pti.tilebox ();
              RealBox tile_region ( bx, Geom(lev).CellSize (), Geom(lev).ProbLo() );

              if ( tile_region.intersects ( avg_region ) )
                {
                  const int np         = NumberOfParticles(pti);
                  const AoS &particles = pti.GetArrayOfStructs();

                  for (int p = 0; p < np; ++p )
                    {
                      if ( avg_region.contains ( &(particles[p].m_rdata.pos[0]) ) )
                        {
                          sum_np++;
//wdf switch U/V/W to X/Y/Z
//                        sum_velx += particles[p].rdata(realData::velx);
//                        sum_vely += particles[p].rdata(realData::vely);
//                        sum_velz += particles[p].rdata(realData::velz);
                          sum_velx += particles[p].pos(0);
                          sum_vely += particles[p].pos(1);
                          sum_velz += particles[p].pos(2);
//wdf
                        }

                    }

                }

            }

          region_np[nr]    = sum_np;
          region_velx[nr]  = sum_velx;
          region_vely[nr]  = sum_vely;
          region_velz[nr]  = sum_velz;
        }

      // Compute parallel reductions
      ParallelDescriptor::ReduceLongSum ( region_np.data(),   nregions );
      ParallelDescriptor::ReduceRealSum ( region_velx.data(), nregions );
      ParallelDescriptor::ReduceRealSum ( region_vely.data(), nregions );
      ParallelDescriptor::ReduceRealSum ( region_velz.data(), nregions );

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
              region_velx[nr] /= region_np[nr];
              region_vely[nr] /= region_np[nr];
              region_velz[nr] /= region_np[nr];

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
                  // Use ios:trunc to delete previous contect
                  ofs.open ( fname.c_str(), ios::out | ios::trunc );
                }
              else
                {
                  // If this is not the first time we write to this file
                  // we append to it
                  ofs.open ( fname.c_str(), ios::out | ios::app );
                }

              // Check if file is good
              if ( !ofs.good() )
                amrex::FileOpenFailed ( fname );

              // Print header if first access
              if ( ncalls == 1 )
                ofs << "#  Time   NP  U  V  W" << std::endl;

              ofs << time << " "
                  << region_np[nr] << " "
                  << region_velx[nr] << " "
                  << region_vely[nr] << " "
                  << region_velz[nr] << std::endl;

              ofs.close();
            }
        }

    }

}

void MFIXParticleContainer::CapSolidsVolFrac(amrex::MultiFab& mf_to_be_filled)
{
    for (MFIter mfi(mf_to_be_filled); mfi.isValid(); ++mfi) {
       const Box& sbx = mf_to_be_filled[mfi].box();

//       const Real max_pack = 0.42;
       const Real max_pack = 0.21;

       Array4<Real> const& ep_g = mf_to_be_filled.array(mfi);

// These lines are commented because this code represents a functionality which
// may be added in future developments
//       AMREX_HOST_DEVICE_FOR_3D(sbx, i, j, k,
//       {
//         ep_g(i,j,k) = std::max(max_pack, ep_g(i,j,k));
//       });
    }
}

void MFIXParticleContainer::set_particle_properties(int pstate, Real pradius, Real pdensity,
                                                    Real& pvol, Real& pmass, Real& omoi, Real& omega)
{
    pvol  = (4.0/3.0)*M_PI*(pradius*pradius*pradius);
    pmass = pvol * pdensity;
    omoi  = 2.5/(pmass * (pradius*pradius));
    omega = 0.0;
}
