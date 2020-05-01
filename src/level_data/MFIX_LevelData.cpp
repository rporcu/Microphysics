#include <MFIX_LevelData.H>

using namespace amrex;

LevelData::LevelData (BoxArray const& ba,
                      DistributionMapping const& dmap,
                      const int nghost,
                      FabFactory<FArrayBox> const& factory)
  : ep_g(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , ep_go(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , p_g(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , p_go(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , T_g(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , T_go(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , ro_g(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , ro_go(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , trac(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , trac_o(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , vel_g(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , vel_go(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , p0_g(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , gp(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , mu_g(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , vort(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , drag(new MultiFab(ba, dmap, 4, nghost, MFInfo(), factory))
  , xslopes_u(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , yslopes_u(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , zslopes_u(new MultiFab(ba, dmap, 3, nghost, MFInfo(), factory))
  , xslopes_s(new MultiFab(ba, dmap, 2, nghost, MFInfo(), factory))
  , yslopes_s(new MultiFab(ba, dmap, 2, nghost, MFInfo(), factory))
  , zslopes_s(new MultiFab(ba, dmap, 2, nghost, MFInfo(), factory))
  , diveu(new MultiFab(amrex::convert(ba, IntVect{1,1,1}), dmap, 1, nghost, MFInfo(), factory))
  , mac_phi(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , mac_rhs(new MultiFab(ba, dmap, 1, nghost, MFInfo(), factory))
  , u_mac(new MultiFab(BoxArray(ba).surroundingNodes(0), dmap, 1, 2, MFInfo(), factory))
  , v_mac(new MultiFab(BoxArray(ba).surroundingNodes(1), dmap, 1, 2, MFInfo(), factory))
  , w_mac(new MultiFab(BoxArray(ba).surroundingNodes(2), dmap, 1, 2, MFInfo(), factory))
{}

void LevelData::resetValues (const amrex::Real covered_val)
{
  ep_g->setVal(1);
  ep_go->setVal(1);
  T_g->setVal(0);
  T_go->setVal(0);
  ro_g->setVal(0);
  ro_go->setVal(0);
  trac->setVal(0);
  trac_o->setVal(0);
  p0_g->setVal(0);
  p_g->setVal(0);
  p_go->setVal(0);
  diveu->setVal(0);
  gp->setVal(0);
  mu_g->setVal(0);
  vel_g->setVal(0);
  vel_go->setVal(0);
  vort->setVal(0);
  drag->setVal(0);
  mac_rhs->setVal(0);
  mac_phi->setVal(0);
  xslopes_u->setVal(0);
  xslopes_s->setVal(0);
  yslopes_u->setVal(0);
  yslopes_s->setVal(0);
  zslopes_u->setVal(0);
  zslopes_s->setVal(0);
  u_mac->setVal(covered_val);
  v_mac->setVal(covered_val);
  w_mac->setVal(covered_val);
}

LevelData::~LevelData ()
{
  delete ep_g;
  delete ep_go;
  delete p_g;
  delete p_go;
  delete T_g;
  delete T_go;
  delete ro_g;
  delete ro_go;
  delete trac;
  delete trac_o;
  delete vel_g;
  delete vel_go;
  delete p0_g;
  delete gp;
  delete mu_g;
  delete vort;
  delete drag;
  delete xslopes_u;
  delete yslopes_u;
  delete zslopes_u;
  delete xslopes_s;
  delete yslopes_s;
  delete zslopes_s;
  delete diveu;
  delete mac_phi;
  delete mac_rhs;
  delete u_mac;
  delete v_mac;
  delete w_mac;
}
