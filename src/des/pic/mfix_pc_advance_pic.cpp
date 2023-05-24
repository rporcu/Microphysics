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
      const int idx_X_sn = m_runtimeRealData.X_sn;
      const int idx_cp_s = m_runtimeRealData.cp_s;
      const int idx_temperature = m_runtimeRealData.temperature;
      const int idx_convection = m_runtimeRealData.convection;
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
           is_IOProc,abstol,reltol,maxiter,idx_cp_s,idx_temperature,idx_convection]
        AMREX_GPU_DEVICE (int i) noexcept
      {
        auto& p = pstruct[i];

        // Get current particle's mass
        const Real p_mass_old = p_realarray[SoArealData::mass][i];
        Real p_mass_new(p_mass_old);

        // Flag to stop computing particle's quantities if mass_new < 0, i.e.
        // the particle disappears because of chemical reactions
        int proceed = 1;
        Real coeff = 1.;

        //*********************************************************************
        // First step: update parcels' mass and density
        //*********************************************************************
        if (update_mass) {

          part_mass_update(p, ptile_data, p_realarray, i, idx_X_sn, nspecies_s,
              dt, p_mass_old, p_mass_new, proceed, coeff, idx_mass_txfr);
        }

        if (proceed) {

          //*********************************************************************
          // Second step: update parcels' temperature
          //*********************************************************************
          if (solve_enthalpy) {

            part_enthalpy_update(ptile_data, p_realarray, i, idx_X_sn, nspecies_s,
                solid_is_a_mixture, solids_parms, dt, coeff, p_mass_new, nullptr,
                enthalpy_source, solve_reactions, idx_cp_s, idx_temperature,
                idx_convection, idx_h_txfr, abstol, reltol, maxiter, is_IOProc);
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
