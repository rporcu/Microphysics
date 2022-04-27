#include <AMReX.H>
#include <mfix.H>

#include <mfix_pic_parms.H>

using namespace amrex;

void mfix::EvolveParcels (Real dt,
                          Real /*time*/,
                          RealVect& gravity_in,
                          const int /*ls_refinement_in*/,
                          Vector< MultiFab* >& cost,
                          std::string& knapsack_weight_type_in)

{
  BL_PROFILE_REGION_START("MFIX_PIC::EvolveParcels()");
  BL_PROFILE("mfix::EvolveParcels()");

  amrex::Print() << "\nEvolving PIC parcels.";

  if ( PIC::verbose > 0 ) {
#if MFP_DISABLED
    amrex::Print() << "\nMean free path limiting is disabled.\n";
#elif SCALE_MFP_DISABLED
    amrex::Print() << "\nMean free path limiting is enabled without scaling.\n";
#elif SCALE_MFP_CHAPMAN
    amrex::Print() << "\nMean free path limiting is enabled with Chapman scaling.\n";
#elif SCALE_MFP_CUTCHIS
    amrex::Print() << "\nMean free path limiting is enabled with Cutchis scaling.\n";
#elif SCALE_MFP_MUSSER
    amrex::Print() << "\nMean free path limiting is enabled with Musser scaling.\n";
#else
    amrex::Abort("Invalid MFP -- check setting in mfix_pic_K.H");
#endif
    amrex::Print() << "Bulk velocity frame of reference: " << PIC::vel_ref_frame << "\n";
    amrex::Print() << "Predicting average mixture velocity at time level n+" << PIC::advance_vel_p << "\n";
  }
  // MultiFabs that hold PIC specific field data all on particle grids.

  const int pic_nghost = 2;
  const int pic_ncomp  = 2;  // (vel_s + mass)

  // NOTE:Ideally we want all/most of these in the same MultiFab to
  // minimize calls to FillBoundary. However, we are limited given
  // the structure of the mfix_redistribute_deposition.
  Vector< MultiFab* >  ep_s(nlev, nullptr);
  Vector< MultiFab* > resid(nlev, nullptr);
  Vector< Array<MultiFab*,3> > vel_s(nlev, {nullptr, nullptr, nullptr});

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
    ep_s[lev]->setVal(0.0, 0, 1, pic_nghost);

    if ( PIC::verbose > 0 ) {
      resid[lev] = new MultiFab(pba, pdm, 1, pic_nghost, MFInfo(), *particle_ebfactory[lev]);
    }
  }

  // This might not be correct for more than one level, so let's
  // make sure that there is only one for now.
  AMREX_ALWAYS_ASSERT(nlev == 1);

  // Parcels should be on level 0 and we don't want
  // them seeing a finer reconstruction of the boundary
  // because it may lead to overpacking.
  const int ls_refinement = 1;
  const MultiFab* ls_data = level_sets[0].get();


  if ( PIC::verbose > 0 ) {
    if ( PIC::initial_step == PIC::InitialStepType::zero_eps ) {
      amrex::Print() << "Initial PIC step using zero ep_s (no stress).\n";
    } else if ( PIC::initial_step == PIC::InitialStepType::nth_eps ) {
      amrex::Print() << "Initial PIC step using n^th step ep_s.\n";
    } else if ( PIC::initial_step == PIC::InitialStepType::taylor_approx ) {
      amrex::Print() << "Initial PIC step using Taylor approximate for ep_s.\n";
    }
  }

  { // Compute the initial guess for the solids volume fraction and
    // deposit the parcel velocity to the grid.
    const bool update_parcels(false);
    const bool apply_forces(false);
    const bool use_taylor_approx( PIC::initial_step == PIC::InitialStepType::taylor_approx);
    const Real advance_vel_p(0.0);

    mfix::pic_iteration(apply_forces, update_parcels, use_taylor_approx,
             advance_vel_p, dt, gravity_in, vel_s, ep_s,
             particle_ebfactory[0].get(), ls_refinement, ls_data);

    // Clear out the deposited volume fraction if initial setp is zero_eps.
    if ( PIC::initial_step == PIC::InitialStepType::zero_eps ) {
      for (int lev = 0; lev < nlev; ++lev) {
        ep_s[lev]->setVal(0.0, 0, 1, pic_nghost);
      }
    }
  } // End initial step block


  // A little bit of information about max solids volume fraction.
  if ( PIC::verbose > 0 ) {
    for (int lev = 0; lev < nlev; ++lev) {
      const Real max_eps = ep_s[lev]->norm0(0,0,false,true);
      amrex::Print().SetPrecision(8) << "\nmax(eps(0)) = " << max_eps;
      if ( PIC::verbose > 1 ) {
        const IntVect max_eps_cell = mfix_locate_max_eps(ep_s, max_eps);
        amrex::Print() << "  at cell  " << max_eps_cell;
      }
      amrex::Print() << "\n";
      MultiFab::Copy(*resid[lev], *ep_s[lev], 0, 0, 1, pic_nghost);
    }
  }

  const bool apply_forces(true);
  const bool use_taylor_approx(false);
  const Real advance_vel_p(PIC::advance_vel_p);

  for(int iter(1); iter <= PIC::max_iter; ++iter) {

    // Calculate the particle normal stress using the local solids
    // volume fraction. This stress is a field variable whose gradient
    // is interpolated to a parcel's position.
    MFIX_CalcSolidsStress(ep_s, gravity_in, cost, knapsack_weight_type_in);

    const bool update_parcels = (iter == PIC::max_iter) ? true : false;

    mfix::pic_iteration(apply_forces, update_parcels, use_taylor_approx,
             advance_vel_p, dt, gravity_in, vel_s, ep_s,
             particle_ebfactory[0].get(), ls_refinement, ls_data);

    if ( PIC::verbose > 0 ) {
      for (int lev = 0; lev < nlev; ++lev) {
        const Real max_eps = ep_s[lev]->norm0(0,0,false,true);
        amrex::Print().SetPrecision(8) << "max(eps(" << iter << ")) = " << max_eps;
        if ( PIC::verbose > 1 ) {
          const IntVect max_eps_cell = mfix_locate_max_eps(ep_s, max_eps);
          amrex::Print() << "  at cell  " << max_eps_cell;
        }
        MultiFab::Xpay(*resid[lev], amrex::Real(-1.0), *ep_s[lev], 0, 0, 1, 0);
        amrex::Print().SetPrecision(4)
          << "    max(eps(" << iter << ") - eps(" << iter-1 << ")) = "
          << resid[lev]->norm0(0,0,false,true) << "\n";
        MultiFab::Copy(*resid[lev], *ep_s[lev], 0, 0, 1, pic_nghost);
      }
    } // PIC verbose
  }


  if ((solids.solve_species && reactions.solve) || fluid.solve_enthalpy) {

    pc->MFIX_PC_AdvanceParcels(dt, cost, knapsack_weight_type_in);
  }

  // exit(0);

  // Account for cross process movements.
  pc->Redistribute(0, 0, 0, 1);

  // Clean up the temporary fluid volume fraction
  for (int lev = 0; lev < nlev; lev++){
    if( resid[lev]    != nullptr ) delete resid[lev];
    if( ep_s[lev]     != nullptr ) delete ep_s[lev];
    if( vel_s[lev][0] != nullptr ) delete vel_s[lev][0];
    if( vel_s[lev][1] != nullptr ) delete vel_s[lev][1];
    if( vel_s[lev][2] != nullptr ) delete vel_s[lev][2];
  }
  amrex::Print() << std::endl;
  BL_PROFILE_REGION_STOP("MFIX_PIC::EvolveParcels()");
}
