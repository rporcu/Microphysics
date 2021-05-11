#include <AMReX.H>
#include <mfix.H>

#include <mfix_pic_parms.H>

using namespace amrex;

void mfix::EvolveParcels (Real dt,
                          Real /*time*/,
                          RealVect& gravity_in,
                          const int ls_refinement_in,
                          Vector< MultiFab* >& cost,
                          std::string& knapsack_weight_type_in,
                          const int advect_enthalpy_in,
                          const Real enthalpy_source)

{
  BL_PROFILE_REGION_START("MFIX_PIC::EvolveParcels()");
  BL_PROFILE("mfix::EvolveParcels()");

  amrex::Print() << "Evolving PIC parcels\n" << std::endl;

  // MultiFabs that hold PIC specific field data all on particle grids.

  const int pic_nghost = 2;
  const int pic_ncomp  = 2;  // (vel_s + mass)

  // NOTE:Ideally we want all/most of these in the same MultiFab to
  // minimize calls to FillBoundary. However, we are limited given
  // the structure of the mfix_redistribute_deposition.
  Vector< MultiFab* > ep_s(nlev);
  Vector< Array<MultiFab*,3> > vel_s(nlev);

  for (int lev(0); lev < nlev; lev ++ ){

    const BoxArray&            pba = pc->ParticleBoxArray(lev);
    const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

    // Averaged particle velocity
    vel_s[lev][0] = new MultiFab(BoxArray(pba).surroundingNodes(0), pdm, pic_ncomp,
                          pic_nghost, MFInfo(), *particle_ebfactory[lev]);
    vel_s[lev][1] = new MultiFab(BoxArray(pba).surroundingNodes(1), pdm, pic_ncomp,
                          pic_nghost, MFInfo(), *particle_ebfactory[lev]);
    vel_s[lev][2] = new MultiFab(BoxArray(pba).surroundingNodes(2), pdm, pic_ncomp,
                          pic_nghost, MFInfo(), *particle_ebfactory[lev]);

    vel_s[lev][0]->setVal(0.0, 0, pic_ncomp, pic_nghost);
    vel_s[lev][1]->setVal(0.0, 0, pic_ncomp, pic_nghost);
    vel_s[lev][2]->setVal(0.0, 0, pic_ncomp, pic_nghost);

    // Solids volume fraction.
    ep_s[lev] = new MultiFab(pba, pdm, 1, pic_nghost, MFInfo(), *particle_ebfactory[lev]);

  }

  // Calculate the averaged solids velocity. This is used to assess how
  // a parcel is moving relative to the bulk solids motion.  A consequence
  // of the deposition is that we get ep_s on the correct boxes.
  MFIX_CalcAvgSolidsVel(vel_s);

  {
    // This might not be correct for more than one level, so let's
    // make sure that there is only one for now. Another note, if
    // the step the advance_eps_step is zero, we could use the
    // existing volume fraction field as eps = 1-ep_g.
    AMREX_ALWAYS_ASSERT(nlev == 1);

    const Real advance_eps_step = PIC::advance_eps_step;

    const int ls_refinement = 1;
    const MultiFab* ls_data = level_sets[0].get();

    amrex::Print() << "Predicting new volume frac at time level n+" << advance_eps_step << "\n";

    Real loc_dt = advance_eps_step*dt;
    mfix_predict_volume_fraction(loc_dt, gravity_in,
                                 particle_ebfactory[0].get(),
                                 ls_refinement, ls_data,
                                 ep_s);
  }

  // Calculate the solids stress gradient using the local solids
  // volume fraction. The solids stress is a field variable that
  // we will later interpolate to the parcel's position.
  MFIX_CalcSolidsStress(ep_s, cost, knapsack_weight_type_in);




  // Move the parcels.
  pc->MFIX_PC_AdvanceParcels(dt, gravity_in, vel_s,
                             cost, knapsack_weight_type_in,
                             advect_enthalpy_in,
                             enthalpy_source);

  // Account for cross process movements.
  pc->Redistribute(0, 0, 0, 0);

  // Now account for solid walls.
  for (int lev = 0; lev < nlev; lev ++ ) {

    if (nlev == 1){
      const int ls_refinement = 1;
      const MultiFab* ls_data = level_sets[lev].get();

      pc->MFIX_PC_ImposeWalls(lev, particle_ebfactory[lev].get(),
                              ls_refinement, ls_data,
                              cost[lev], knapsack_weight_type_in);

    } else {

      const int ls_refinement = ls_refinement_in;
      const MultiFab* ls_data = level_sets[1].get();

      pc->MFIX_PC_ImposeWalls(nlev, particle_ebfactory[lev].get(),
                              ls_refinement, ls_data,
                              cost[lev], knapsack_weight_type_in);
    }

  }

  // This redistribute might not be needed, but it's done "just in case"
  // that accounting for a boundary collision doesn't position a parcel
  // on another process.
  pc->Redistribute(0, 0, 0, 1);

  // Clean up the temporary fluid volume fraction
  for (int lev = 0; lev < nlev; lev++){
    delete ep_s[lev];
    delete vel_s[lev][0];
    delete vel_s[lev][1];
    delete vel_s[lev][2];
  }

    BL_PROFILE_REGION_STOP("MFIX_PIC::EvolveParcels()");
}
