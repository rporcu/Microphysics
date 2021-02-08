#include <mfix_leveldata.H>
#include <mfix_eb_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_reactions_parms.H>

using namespace amrex;

LevelData::LevelData (BoxArray const& ba,
                      DistributionMapping const& dmap,
                      const int nghost,
                      FabFactory<FArrayBox> const& factory)
  : ep_g(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , ep_go(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , p_g(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , p_go(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , ro_g(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , ro_go(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , MW_g(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , trac(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , trac_o(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , vel_g(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , vel_go(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , p0_g(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , gp(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , mu_g(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , T_g(nullptr)
  , T_go(nullptr)
  , cp_g(nullptr)
  , k_g(nullptr)
  , h_g(nullptr)
  , h_go(nullptr)
  , X_gk(nullptr)
  , X_gko(nullptr)
  , D_gk(nullptr)
  , cp_gk(nullptr)
  , h_gk(nullptr)
  , vort(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , txfr(new MultiFab(ba, dmap, Transfer::count, nghost, MFInfo(), factory))
  , ro_gk_txfr(nullptr)
  , diveu(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , mac_phi(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , divtau_o(new MultiFab(ba, dmap, 3, 0, MFInfo(), factory))
{

  if (FLUID::solve_enthalpy) {
    T_g  = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    T_go = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    cp_g = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    k_g  = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    h_g  = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    h_go = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);

    if (EB::fix_temperature) {
      T_g_on_eb = new MultiFab(ba, dmap, 1, 1, MFInfo(), factory);
      k_g_on_eb = new MultiFab(ba, dmap, 1, 1, MFInfo(), factory);
    }
  }

  if (FLUID::solve_species) {
    X_gk  = new MultiFab(ba, dmap, FLUID::nspecies, nghost, MFInfo(), factory);
    X_gko = new MultiFab(ba, dmap, FLUID::nspecies, nghost, MFInfo(), factory);
    D_gk  = new MultiFab(ba, dmap, FLUID::nspecies, nghost, MFInfo(), factory);
  }

  if (FLUID::solve_enthalpy and FLUID::solve_species) {
    cp_gk  = new MultiFab(ba, dmap, FLUID::nspecies, nghost, MFInfo(), factory);
    h_gk  = new MultiFab(ba, dmap, FLUID::nspecies, nghost, MFInfo(), factory);
  }

  if (REACTIONS::solve and FLUID::solve_species) {
    ro_gk_txfr = new MultiFab(ba, dmap, FLUID::nspecies, nghost, MFInfo(), factory);
  }
}

void LevelData::resetValues (const amrex::Real covered_val)
{
  ep_g->setVal(1);
  ep_go->setVal(1);
  p_g->setVal(0);
  p_go->setVal(0);
  ro_g->setVal(0);
  ro_go->setVal(0);
  MW_g->setVal(0);
  trac->setVal(0);
  trac_o->setVal(0);
  vel_g->setVal(0);
  vel_go->setVal(0);
  p0_g->setVal(0);
  gp->setVal(0);
  mu_g->setVal(0);
  vort->setVal(0);
  txfr->setVal(0);
  diveu->setVal(0);
  mac_phi->setVal(0);
  divtau_o->setVal(0);

  if (FLUID::solve_enthalpy) {
    T_g->setVal(0);
    T_go->setVal(0);
    cp_g->setVal(0);
    k_g->setVal(0);
    h_g->setVal(0);
    h_go->setVal(0);

    if (EB::fix_temperature) {
      T_g_on_eb->setVal(0);
      k_g_on_eb->setVal(0);
    }
  }

  if (FLUID::solve_species) {
    X_gk->setVal(0);
    X_gko->setVal(0);
    D_gk->setVal(0);
  }

  if (FLUID::solve_enthalpy and FLUID::solve_species) {
    cp_gk->setVal(0);
    h_gk->setVal(0);
  }

  if (REACTIONS::solve and FLUID::solve_species) {
    ro_gk_txfr->setVal(0);
  }
}

LevelData::~LevelData ()
{
  delete ep_g;
  delete ep_go;
  delete p_g;
  delete p_go;
  delete ro_g;
  delete ro_go;
  delete MW_g;
  delete trac;
  delete trac_o;
  delete vel_g;
  delete vel_go;
  delete p0_g;
  delete gp;
  delete mu_g;
  delete vort;
  delete txfr;
  delete diveu;
  delete mac_phi;
  delete divtau_o;

  if (FLUID::solve_enthalpy) {
    delete T_g;
    delete T_go;
    delete cp_g;
    delete k_g;
    delete h_g;
    delete h_go;

    if (EB::fix_temperature) {
      delete T_g_on_eb;
      delete k_g_on_eb;
    }
  }

  if (FLUID::solve_species) {
    delete X_gk;
    delete X_gko;
    delete D_gk;
  }

  if (FLUID::solve_enthalpy and FLUID::solve_species) {
    delete cp_gk;
    delete h_gk;
  }

  if (REACTIONS::solve and FLUID::solve_species) {
    delete ro_gk_txfr;
  }
}
