#include <mfix.H>
#include <mfix_bc.H>
#include <mfix_algorithm.H>
#include <mfix_solvers.H>

using namespace amrex;

void MFIXParticleContainer::MFIX_PC_AdvanceParcels (Real dt,
                                                    Vector< MultiFab* >& cost,
                                                    std::string& knapsack_weight_type)
{

  BL_PROFILE("MFIXParticleContainer::MFIX_PC_AdvanceParcels()");

  const int is_IOProc = int(ParallelDescriptor::IOProcessor());

  const Real abstol = newton_abstol;
  const Real reltol = newton_reltol;
  const int maxiter = newton_maxiter;

  for (int lev = 0; lev < nlev; lev ++ )
  {

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {

      // Timer used for load-balancing
      Real wt = ParallelDescriptor::second();

      PairIndex index(pti.index(), pti.LocalTileIndex());

      const int nrp = GetParticles(lev)[index].numRealParticles();

      auto& plev = GetParticles(lev);
      auto& ptile = plev[index];
      auto& aos   = ptile.GetArrayOfStructs();
      ParticleType* pstruct = aos().dataPtr();

      auto ptile_data = ptile.getParticleTileData();

      auto& soa = ptile.GetStructOfArrays();
      auto p_realarray = soa.realarray();
      auto p_intarray  = soa.intarray();

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR("pic_time_march()", pic_time_march);
#endif
      /********************************************************************
       * Move particles based on collision forces and torques             *
       *******************************************************************/

      const int nspecies_s = solids.nspecies();
      const int nreactions = reactions.nreactions();

      // Particles SoA starting indexes for mass fractions and rate of
      // formations
      const SoAruntimerealData runtimedata_idxs = m_runtimedata_idxs;

      const int update_mass     = solids.update_mass() && solids.solve_species() && reactions.solve();
      const int solve_enthalpy  = solids.solve_enthalpy();
      const int solve_reactions = reactions.solve();

      const Real enthalpy_source = solids.enthalpy_source();

      const int solid_is_a_mixture = solids.isMixture();

      const auto& bc_parms = m_boundary_conditions.parameters();
      const auto& solids_parms = solids.parameters();

      MFIXSolvers::Newton nonlinear_solver(abstol, reltol, maxiter, is_IOProc);

      ParticleUpdates particle_updates(pstruct, p_realarray, p_intarray,
          ptile_data, bc_parms, solids_parms);

      amrex::ParallelFor(nrp,
          [nrp,p_realarray,dt,particle_updates,ptile_data,nspecies_s,nreactions,
           runtimedata_idxs,update_mass,solve_reactions,solid_is_a_mixture,
           solve_enthalpy,enthalpy_source,solids_parms,nonlinear_solver]
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

        // Flag to stop computing particle's quantities if mass_new < 0, i.e.
        // the particle disappears because of chemical reactions
        int proceed = 1;

        //*********************************************************************
        // First step: update parcels' mass and density
        //*********************************************************************
        if (update_mass) {

          const int solve_dem = 0;

          particle_updates.update_mass(i, nspecies_s, X_sn, dt,
              runtimedata_idxs, proceed, solve_dem);
        }

        if (proceed) {

          const Real p_mass_new = p_realarray[SoArealData::mass][i];
          const Real mass_coeff = update_mass ? p_mass_old/p_mass_new : 1.;

          //*********************************************************************
          // Second step: update parcels' temperature
          //*********************************************************************
          if (solve_enthalpy) {

            particle_updates.update_enthalpy(i, nspecies_s, X_sn, solve_reactions,
                solid_is_a_mixture, mass_coeff, nullptr, dt, runtimedata_idxs,
                enthalpy_source, nonlinear_solver);

          }
        }
      });

      /********************************************************************
       * Update runtime cost (used in load-balancing)                     *
       *******************************************************************/
      if (cost[lev])
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
        (*cost[lev])[pti].plus<RunOn::Device>(wt, tbx);
      }

      Gpu::synchronize();

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR_STOP(pic_time_march);
#endif

    } // particle-tile iterator

  } // loop over levels

}
