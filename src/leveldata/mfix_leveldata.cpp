#include <mfix_leveldata.H>
#include <mfix_eb_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_reactions_parms.H>

using namespace amrex;

LevelData::LevelData (BoxArray const& ba,
                      DistributionMapping const& dmap,
                      const int nghost,
                      FabFactory<FArrayBox> const& factory,
                      const int _solve_enthalpy,
                      const int _solve_species,
                      const int _nspecies_g)
  : ep_g(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , p_g(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , p_go(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , ro_g(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , ro_go(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , trac(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , trac_o(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , vel_g(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , vel_go(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , p0_g(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , gp(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , pressure_g(nullptr)
  , pressure_go(nullptr)
  , T_g(nullptr)
  , T_go(nullptr)
  , h_g(nullptr)
  , h_go(nullptr)
  , X_gk(nullptr)
  , X_gko(nullptr)
  , vort(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , txfr(new MultiFab(ba, dmap, Transfer::count, nghost, MFInfo(), factory))
  , chem_txfr(nullptr)
  , diveu(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , mac_phi(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , divtau_o(new MultiFab(ba, dmap, 3, 0, MFInfo(), factory))
  , level_allocated(1)
  , solve_enthalpy(_solve_enthalpy)
  , solve_species(_solve_species)
  , nspecies_g(_nspecies_g)
{
  if (solve_enthalpy) {
    pressure_g  = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    pressure_go = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    T_g  = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    T_go = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    h_g  = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    h_go = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);

    if (EB::fix_temperature) {
      T_g_on_eb = new MultiFab(ba, dmap, 1, 1, MFInfo(), factory);
    }
  }

  if (solve_species) {
    X_gk  = new MultiFab(ba, dmap, nspecies_g, nghost, MFInfo(), factory);
    X_gko = new MultiFab(ba, dmap, nspecies_g, nghost, MFInfo(), factory);
  }

  if (REACTIONS::solve && solve_species) {
    ChemTransfer chem_txfr_idxs(nspecies_g, REACTIONS::nreactions);
    chem_txfr = new MultiFab(ba, dmap, chem_txfr_idxs.count, nghost, MFInfo(), factory);
  }
}

void LevelData::resetValues (const amrex::Real init_value)
{
  ep_g->setVal(1);
  p_g->setVal(init_value);
  p_go->setVal(init_value);
  ro_g->setVal(init_value);
  ro_go->setVal(init_value);
  trac->setVal(0);
  trac_o->setVal(init_value);
  vel_g->setVal(init_value);
  vel_go->setVal(init_value);
  p0_g->setVal(init_value);
  gp->setVal(0);
  vort->setVal(init_value);
  txfr->setVal(0);
  diveu->setVal(init_value);
  mac_phi->setVal(init_value);
  divtau_o->setVal(init_value);

  if (solve_enthalpy) {
    pressure_g->setVal(init_value);
    pressure_go->setVal(init_value);
    T_g->setVal(init_value);
    T_go->setVal(init_value);
    h_g->setVal(init_value);
    h_go->setVal(init_value);

    if (EB::fix_temperature) {
      T_g_on_eb->setVal(init_value);
    }
  }

  if (solve_species) {
    X_gk->setVal(init_value);
    X_gko->setVal(init_value);
  }

  if (REACTIONS::solve && solve_species) {
    chem_txfr->setVal(init_value);
  }
}

LevelData::~LevelData ()
{
  if (level_allocated) {
    delete ep_g;
    delete p_g;
    delete p_go;
    delete ro_g;
    delete ro_go;
    delete trac;
    delete trac_o;
    delete vel_g;
    delete vel_go;
    delete p0_g;
    delete gp;
    delete vort;
    delete txfr;
    delete diveu;
    delete mac_phi;
    delete divtau_o;

    if (solve_enthalpy) {
      delete pressure_g;
      delete pressure_go;
      delete T_g;
      delete T_go;
      delete h_g;
      delete h_go;

      if (EB::fix_temperature) {
        delete T_g_on_eb;
      }
    }

    if (solve_species) {
      delete X_gk;
      delete X_gko;
    }

    if (REACTIONS::solve && solve_species) {
      delete chem_txfr;
    }
  }
}
