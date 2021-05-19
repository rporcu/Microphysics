#include <mfix.H>
#include <mfix_filcc.H>

#include <AMReX_FillPatchUtil.H>
#include <mfix_pc.H>
#include <mfix_mf_helpers.H>

#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>

void mfix::mfix_predict_volume_fraction (Real dt_in,
                                         RealVect& gravity_in,
                                         EBFArrayBoxFactory* ebfactory_in,
                                         const int ls_refinement,
                                         const MultiFab* ls_phi,
                                         Vector< MultiFab* >& ep_s_out)
{
  BL_PROFILE("mfix::mfix_predict_volume_fraction()");

  if (nlev > 2)
    amrex::Abort("For right now mfix::mfix_calc_volume_fraction can only handle up to 2 levels");

  for (int lev = 0; lev < nlev; lev++) {

    // We must have ghost cells for each FAB so that a particle in one grid can spread
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (ep_s_out[lev]->nGrow() < 1)
      amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

    ep_s_out[lev]->setVal(0.0,0,1,ep_s_out[lev]->nGrow());
  }

  const Geometry& gm  = Geom(0);
  const FabArray<EBCellFlagFab>* flags = nullptr;
  const MultiFab* volfrac = nullptr;

  for (int lev = 0; lev < nlev; lev++) {

    // Use level 0 to define the EB factory
    if (lev == 0) {

      flags   = &(particle_ebfactory[lev]->getMultiEBCellFlagFab());
      volfrac = &(particle_ebfactory[lev]->getVolFrac());

    } else {

      Vector<int> ngrow = {1,1,1};
      EBFArrayBoxFactory* crse_factory;

      crse_factory = (makeEBFabFactory(gm, ep_s_out[lev]->boxArray(),
                                      ep_s_out[lev]->DistributionMap(),
                                      ngrow, EBSupport::volume)).release();

      flags   = &(crse_factory->getMultiEBCellFlagFab());
      volfrac = &(crse_factory->getVolFrac());

      delete crse_factory;
    }

    // Deposit particle volume to the grid
    pc->PredictPICVolumeDeposition(lev, dt_in, gravity_in,
            ebfactory_in, ls_refinement, ls_phi,
            *ep_s_out[lev], volfrac, flags);
  }

  {
    // The deposition occurred on level 0, thus the next few operations
    // only need to be carried out on level 0.
    int lev(0);

    // Move any volume deposited outside the domain back into the domain
    // when BC is either a pressure inlet or mass inflow.
    mfix_deposition_bcs(lev, *ep_s_out[lev]);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    ep_s_out[lev]->SumBoundary(gm.periodicity());

    // Fill the boundaries so we calculate the correct average
    // solids volume fraction for periodic boundaries.
    ep_s_out[lev]->FillBoundary(gm.periodicity());

    // Create a copy of the solids volume fraction to use
    // in the redistribution.
    MultiFab* eps_tmp;
    eps_tmp = (MFHelpers::createFrom(*ep_s_out[lev])).release();

    // Clear the grid boundaries to prepare for redistribution which
    // could put material in your ghost cells. Do this AFTER
    // we created the copy so the copy has the correct boundary info.
    ep_s_out[lev]->setBndry(0.0);


    // Move excessive solids volume from small cells to neighboring cells. A copy
    // of the deposition field is made so that when an average is calc
    mfix_redistribute_deposition(lev, *eps_tmp, *ep_s_out[lev], volfrac, flags,
                                  mfix::m_max_solids_volume_fraction);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    ep_s_out[lev]->SumBoundary(gm.periodicity());
    ep_s_out[lev]->FillBoundary(gm.periodicity());

    // we no longer need the copy.
    delete eps_tmp;

  }


}
