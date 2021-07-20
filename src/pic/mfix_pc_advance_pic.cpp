#include <mfix.H>
#include <mfix_bc_parms.H>
#include <mfix_algorithm.H>
#include <mfix_solvers.H>

using namespace amrex;
using namespace Solvers;

void MFIXParticleContainer::MFIX_PC_AdvanceParcels (Real dt,
                                                    const int advect_enthalpy,
                                                    const Real enthalpy_source,
                                                    Vector< MultiFab* >& cost,
                                                    std::string& knapsack_weight_type)
{

  BL_PROFILE("MFIXParticleContainer::MFIX_PC_AdvanceParcels()");

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

      const int nspecies_s = solids.nspecies;
      const int nreactions = REACTIONS::nreactions;

      // Particles SoA starting indexes for mass fractions and rate of
      // formations
      const int idx_X_sn        = m_runtimeRealData.X_sn;
      const int idx_ro_sn_txfr  = m_runtimeRealData.ro_sn_txfr;
      const int idx_vel_s_txfr  = m_runtimeRealData.vel_s_txfr;
      const int idx_h_s_txfr    = m_runtimeRealData.h_s_txfr;

      const int update_mass           = solids.solve_species && REACTIONS::solve;
      const int update_temperature    = advect_enthalpy;
      const int local_advect_enthalpy = advect_enthalpy;
      const int solve_reactions = REACTIONS::solve;

      const int solid_is_a_mixture = solids.is_a_mixture;

      auto& solids_parms = *solids.parameters;

      amrex::ParallelFor(nrp,
        [pstruct,p_realarray,p_intarray,ptile_data,dt,
         nspecies_s,nreactions,idx_X_sn,
         idx_ro_sn_txfr,idx_vel_s_txfr,update_mass,update_temperature,
         solve_reactions,idx_h_s_txfr,solid_is_a_mixture,
         local_advect_enthalpy,enthalpy_source,solids_parms]
        AMREX_GPU_DEVICE (int lp) noexcept
      {
        auto& p = pstruct[lp];

        GpuArray<Real,SPECIES::NMAX> X_sn;

        // Get current particle's mass
        Real p_mass_old = p_realarray[SoArealData::mass][lp];
        Real p_mass_new(p_mass_old);

        // Get current particle's density
        const Real p_density_old = p_realarray[SoArealData::density][lp];
        Real p_density_new(p_density_old);

        // Get current particle's volume
        Real p_vol = p_realarray[SoArealData::volume][lp];

        // Flag to stop computing particle's quantities if mass_new < 0, i.e.
        // the particle disappears because of chemical reactions
        int proceed = 1;

        //*********************************************************************
        // First step: update parcels' mass and density
        //*********************************************************************
        if (update_mass) {
          // Total particle density exchange rate
          Real total_ro_rate(0);

          // Loop over species
          for (int n_s(0); n_s < nspecies_s; n_s++)
          {
            // Current species mass fraction
            X_sn[n_s] = ptile_data.m_runtime_rdata[idx_X_sn+n_s][lp];

            // Get the current reaction rate for species n_s
            const Real ro_sn_rate = ptile_data.m_runtime_rdata[idx_ro_sn_txfr+n_s][lp];

            X_sn[n_s] = X_sn[n_s]*p_density_old + dt*ro_sn_rate;

            // Update the total mass exchange rate
            total_ro_rate += ro_sn_rate;
          }

          // Update the total mass of the particle
          p_density_new = p_density_old + dt * total_ro_rate;

          if (p_density_new > 0) {

            Real total_X(0.);

            // Normalize species mass fractions
            for (int n_s(0); n_s < nspecies_s; n_s++) {
              Real X_sn_new = X_sn[n_s] / p_density_new;

              if (X_sn_new < 0) X_sn_new = 0;
              if (X_sn_new > 1) X_sn_new = 1;

              total_X += X_sn_new;
              X_sn[n_s] = X_sn_new;
            }

            for (int n_s(0); n_s < nspecies_s; n_s++) {
              // Divide updated species mass fractions by total_X
              ptile_data.m_runtime_rdata[idx_X_sn+n_s][lp] = X_sn[n_s] / total_X;
            }

            // Write out to global memory particle's mass and density
            p_realarray[SoArealData::density][lp] = p_density_new;
            p_mass_new = p_density_new * p_vol;
            p_realarray[SoArealData::mass][lp] = p_mass_new;
          } else {
            p.id() = -1;
            proceed = 0;
          }
        }

        if (proceed) {
          //*********************************************************************
          // Second step: update parcels' positions and velocities
          //*********************************************************************

          //*********************************************************************
          // Third step: update parcels' temperature
          //*********************************************************************
          if(local_advect_enthalpy) {
            const int phase = p_intarray[SoAintData::phase][lp];

            const Real Tp_old = p_realarray[SoArealData::temperature][lp];

            const Real cp_s_old = solids_parms.calc_cp_s<RunOn::Gpu>(phase-1,Tp_old);
            Real cp_s_new(0);

            if (solid_is_a_mixture) {
              for (int n_s(0); n_s < nspecies_s; ++n_s)
                cp_s_new += solids_parms.calc_cp_sn<RunOn::Gpu>(Tp_old,n_s) *
                            ptile_data.m_runtime_rdata[idx_X_sn+n_s][lp];

              p_realarray[SoArealData::cp_s][lp] = cp_s_new;
            } else {
              cp_s_new = cp_s_old;
            }

            AMREX_ASSERT(cp_s_new > 0.);

            const Real coeff = update_mass ? (p_mass_old/p_mass_new) : 1.;

            Real p_enthalpy_new =
              coeff*solids_parms.calc_h_s<RunOn::Gpu>(phase-1,Tp_old) +
              dt*((p_realarray[SoArealData::convection][lp]+enthalpy_source)/p_mass_new);

            if (solve_reactions)
              p_enthalpy_new -= dt*p_density_new*ptile_data.m_runtime_rdata[idx_h_s_txfr][lp];

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

                for (int n(0); n < nspecies_s; ++n)
                  hp_loc += X_sn[n]*solids_parms.calc_h_sn<RunOn::Gpu>(Tp_arg,n);
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

                for (int n(0); n < nspecies_s; ++n)
                  gradient += X_sn[n]*solids_parms.calc_partial_h_sn<RunOn::Gpu>(Tp_arg,n);
              }

              return gradient;
            };

            Real Tp_new(Tp_old);

            const Real dumping_factor = 1.;

            DumpedNewton::solve(Tp_new, R, partial_R, dumping_factor, 1.e-5, 1.e-5);

            p_realarray[SoArealData::temperature][lp] = Tp_new;

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
