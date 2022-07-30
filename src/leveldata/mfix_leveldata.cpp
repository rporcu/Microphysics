#include <mfix_leveldata.H>
#include <mfix_eb.H>
#include <mfix_fluid.H>
#include <mfix_species.H>
#include <mfix_reactions.H>

using namespace amrex;

LevelData::LevelData (BoxArray const& ba,
                      DistributionMapping const& dmap,
                      const int nghost,
                      FabFactory<FArrayBox> const& factory,
                      MFIXEmbeddedBoundaries* embedded_boundaries,
                      MFIXFluidPhase* fluid_in,
                      MFIXReactions* reactions_in,
                      bool only_ep_g)
  : level_allocated(1)
  , m_embedded_boundaries(embedded_boundaries)
  , fluid(fluid_in)
  , reactions(reactions_in)
  , ep_g_only(only_ep_g)
{
  ep_g = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);

  if (!only_ep_g) {

    p_g = new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory);
    p_go = new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory);

    ro_g = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    ro_go = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);

    trac = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    trac_o = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);

    vel_g = new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory);
    vel_go = new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory);

    p0_g = new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory);
    gp = new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory);

    vort = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);

    Transfer txfr_idxs(fluid->nspecies(), reactions->nreactions());
    txfr = new MultiFab(ba, dmap, txfr_idxs.count, nghost, MFInfo(), factory);

    diveu = new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory);

    mac_phi = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);

    divtau_o = new MultiFab(ba, dmap, 3, 0, MFInfo(), factory);

    if (reactions->solve() ||
        (fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem ||
         fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem)) {
      thermodynamic_p_g  = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
      thermodynamic_p_go = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    }

    if (fluid->solve_enthalpy() ||
        (fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem ||
         fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem)) {
      T_g  = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
      T_go = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    }

    if (fluid->solve_enthalpy()) {
      h_g  = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
      h_go = new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory);
    }

    if (fluid->solve_enthalpy() && m_embedded_boundaries->fix_temperature()) {
      T_g_on_eb = new MultiFab(ba, dmap, 1, 1, MFInfo(), factory);
    }

    if (fluid->solve_species()) {
      X_gk  = new MultiFab(ba, dmap, fluid->nspecies(), nghost, MFInfo(), factory);
      X_gko = new MultiFab(ba, dmap, fluid->nspecies(), nghost, MFInfo(), factory);
    }

  }
}


void LevelData::resetValues (const amrex::Real init_value)
{
  ep_g->setVal(1);

  if (!ep_g_only) {
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

     if (reactions->solve() ||
         (fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem ||
          fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem)) {
       thermodynamic_p_g->setVal(init_value);
       thermodynamic_p_go->setVal(init_value);
     }

     if (fluid->solve_enthalpy() ||
         (fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem ||
          fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem)) {
       T_g->setVal(init_value);
       T_go->setVal(init_value);
     }

     if (fluid->solve_enthalpy()) {
       h_g->setVal(init_value);
       h_go->setVal(init_value);

       if (m_embedded_boundaries->fix_temperature()) {
         T_g_on_eb->setVal(init_value);
       }
     }

     if (fluid->solve_species()) {
       X_gk->setVal(init_value);
       X_gko->setVal(init_value);
     }
  }

}

LevelData::~LevelData ()
{
  if (level_allocated) {
    delete ep_g;

    if (!ep_g_only) {
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

       if (reactions->solve() ||
           (fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem ||
            fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem)) {
         delete thermodynamic_p_g;
         delete thermodynamic_p_go;
       }

       if (fluid->solve_enthalpy() ||
           (fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasOpenSystem ||
            fluid->constraint_type() == MFIXFluidPhase::ConstraintType::IdealGasClosedSystem)) {
         delete T_g;
         delete T_go;
       }

       if (fluid->solve_enthalpy()) {
         delete h_g;
         delete h_go;

         if (m_embedded_boundaries->fix_temperature()) {
           delete T_g_on_eb;
         }
       }

       if (fluid->solve_species()) {
         delete X_gk;
         delete X_gko;
       }
    }
  }
}
