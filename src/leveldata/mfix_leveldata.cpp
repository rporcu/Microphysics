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
  , ba_proc(new MultiFab(ba, dmap, 1, 0, MFInfo(), factory))
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

void LevelData::resetValues (const amrex::Real covered_val)
{
  ep_g->setVal(1);
  p_g->setVal(0);
  p_go->setVal(0);
  ro_g->setVal(0);
  ro_go->setVal(0);
  trac->setVal(0);
  trac_o->setVal(0);
  vel_g->setVal(0);
  vel_go->setVal(0);
  p0_g->setVal(0);
  gp->setVal(0);
  vort->setVal(0);
  txfr->setVal(0);
  diveu->setVal(0);
  mac_phi->setVal(0);
  divtau_o->setVal(0);
  ba_proc->setVal(0);

  if (solve_enthalpy) {
    pressure_g->setVal(0);
    pressure_go->setVal(0);
    T_g->setVal(0);
    T_go->setVal(0);
    h_g->setVal(0);
    h_go->setVal(0);

    if (EB::fix_temperature) {
      T_g_on_eb->setVal(0);
    }
  }

  if (solve_species) {
    X_gk->setVal(0);
    X_gko->setVal(0);
  }

  if (REACTIONS::solve && solve_species) {
    chem_txfr->setVal(0);
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
    delete ba_proc;

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
