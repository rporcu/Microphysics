#include <AMReX.H>
#include <mfix.H>

using namespace amrex;
using namespace std;

void mfix::EvolveParcels (int nstep,
                          Real dt,
                          Real time,
                          RealVect& gravity,
                          const int ls_refinement_in,
                          Vector< MultiFab* >& cost,
                          std::string& knapsack_weight_type)
{
  BL_PROFILE_REGION_START("MFIX_PIC::EvolveParcels()");
  BL_PROFILE("mfix::EvolveParcels()");

  amrex::Print() << "Evolving PIC parcels" << std::endl;

  // MultiFabs that hold PIC specific field data all on particle grids.

  const int pic_nghost = 1;
  const int pic_ncomp  = 4;  // Four components (3 vel_s + 1 mass)


  // NOTE:Ideally we want all/most of these in the same MultiFab to
  // minimize calls to FillBoundary. However, we are limited given
  // the structure of the mfix_redistribute_deposition.
  Vector< MultiFab* > ep_s(nlev);
  Vector< MultiFab* > avg_prop(nlev);

  for (int lev(0); lev < nlev; lev ++ ){

    const amrex::BoxArray&            pba = pc->ParticleBoxArray(lev);
    const amrex::DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

    // Averaged particle properties
    // 0-2: Mass-averaged solids velocity
    //   3: Volume-averaged solids mass (macroscopic density)
    avg_prop[lev] = new amrex::MultiFab(pba, pdm, pic_ncomp, pic_nghost, MFInfo(), *particle_ebfactory[lev]);
    avg_prop[lev]->setVal(0.0, 0, pic_ncomp, avg_prop[lev]->nGrow());

    // Solids volume fraction.
    ep_s[lev] = new amrex::MultiFab(pba, pdm, 1, pic_nghost, MFInfo(), *particle_ebfactory[lev]);

    // Copy gas phase volume fraction
    ep_s[lev]->copy(*m_leveldata[lev]->ep_g,  0, 0,
                      m_leveldata[lev]->ep_g->nComp(),
                      pic_nghost, pic_nghost);

    // Now define this ep_s = (1 - ep_g)
    ep_s[lev]->mult(-1.0, ep_s[lev]->nGrow());
    ep_s[lev]->plus( 1.0, ep_s[lev]->nGrow());

    ep_s[lev]->FillBoundary(geom[lev].periodicity());
  }

  // Calculate the averaged solids velocity. This is used to assess how
  // a parcel is moving relative to the bulk solids motion.  A consequence
  // of the deposition is that we get ep_s on the correct boxes.
  MFIX_CalcAvgSolidsProp(avg_prop);

  // Calculate the solids stress gradient using the local solids
  // volume fraction. The solids stress is a field variable that
  // we will later interpolate to the parcel's position.
  MFIX_CalcSolidsStress(ep_s, avg_prop);


  // Move the parcels.
  pc->MFIX_PC_AdvanceParcels(dt, gravity, particle_ebfactory, ep_s, avg_prop,
                             cost, knapsack_weight_type);


  // Now account for solid walls.
  for (int lev = 0; lev < nlev; lev ++ ) {

    if (nlev == 1){
      const int ls_refinement = 1;
      const MultiFab* ls_data = level_sets[lev];

      pc->MFIX_PC_ImposeWalls(lev, particle_ebfactory[lev],
                              ls_refinement, ls_data);

    } else {

      const int ls_refinement = ls_refinement_in;
      const MultiFab* ls_data = level_sets[1];

      pc->MFIX_PC_ImposeWalls(nlev, particle_ebfactory[lev],
                              ls_refinement, ls_data);
    }

  }

  pc->Redistribute(0, 0, 0, 0);

  // Clean up the temporary fluid volume fraction
  for (int lev = 0; lev < nlev; lev++){
    delete ep_s[lev];
    delete avg_prop[lev];
  }


    amrex::Print() << "done. \n";

    BL_PROFILE_REGION_STOP("MFIX_PIC::EvolveParcels()");
}

