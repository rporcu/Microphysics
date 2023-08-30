#include <mfix.H>
#include <mfix_pc_updates_K.H>
#include <mfix_bc.H>
#include <mfix_algorithm.H>
#include <mfix_solvers.H>

using namespace amrex;
using namespace Solvers;

void MFIXParticleContainer::MFIX_PC_AdvanceParcels (Real dt,
                                                    Vector< MultiFab* >& cost,
                                                    std::string& knapsack_weight_type)
{

  BL_PROFILE("MFIXParticleContainer::MFIX_PC_AdvanceParcels()");

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
      const int idx_X_sn = m_runtimeRealData.X_sn;
      const int idx_mass_txfr = m_runtimeRealData.mass_txfr;
      const int idx_h_txfr = m_runtimeRealData.h_txfr;

      const int update_mass = solids.update_mass() && solids.solve_species() && reactions.solve();
      const int solve_enthalpy = solids.solve_enthalpy();
      const int solve_reactions = reactions.solve();

      const Real enthalpy_source = solids.enthalpy_source();

      const int solid_is_a_mixture = solids.isMixture();

      const auto& solids_parms = solids.parameters();

      amrex::ParallelFor(nrp,
          [pstruct,p_realarray,p_intarray,ptile_data,dt,nspecies_s,nreactions,
           idx_X_sn,idx_mass_txfr,update_mass,solve_reactions,idx_h_txfr,
           solid_is_a_mixture,solve_enthalpy,enthalpy_source,solids_parms,
           abstol,reltol,maxiter]
        AMREX_GPU_DEVICE (int i) noexcept
      {
        auto& p = pstruct[i];

        // Get current particle's mass
        const Real p_mass_old = p_realarray[SoArealData::mass][i];
        Real p_mass_new(p_mass_old);

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

        // Flag to stop computing particle's quantities if mass_new < 0, i.e.
        // the particle disappears because of chemical reactions
        int proceed = 1;
        Real coeff = 1.;

        //*********************************************************************
        // First step: update parcels' mass and density
        //*********************************************************************
        if (update_mass) {

          part_mass_update(p, ptile_data, p_realarray, i, idx_X_sn, nspecies_s,
              dt, p_mass_old, p_mass_new, nullptr, nullptr, proceed, coeff,
              idx_mass_txfr, 0);
        }

        if (proceed) {

          //*********************************************************************
          // Second step: update parcels' temperature
          //*********************************************************************
          if (solve_enthalpy) {

            part_enthalpy_update(ptile_data, p_realarray, i, idx_X_sn, nspecies_s,
                solid_is_a_mixture, solids_parms, dt, coeff, p_enthalpy_old,
                p_mass_new, nullptr, enthalpy_source, solve_reactions, idx_h_txfr,
                abstol, reltol, maxiter, 0);
          }
        }
      });

      /********************************************************************
       * Update runtime cost (used in load-balancing)                     *
       *******************************************************************/
      if (cost[lev]) {
        UpdateCost(cost[lev], pti, knapsack_weight_type, wt);
      }

      Gpu::synchronize();

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR_STOP(pic_time_march);
#endif

    } // particle-tile iterator

  } // loop over levels

}
