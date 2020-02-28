#include <MFIX_LevelData.H>

using namespace amrex;

LevelData::LevelData ()
  : ep_g(nullptr)
  , ep_go(nullptr)
  , p_g(nullptr)
  , p_go(nullptr)
  , ro_g(nullptr)
  , ro_go(nullptr)
  , trac(nullptr)
  , trac_o(nullptr)
  , vel_g(nullptr)
  , vel_go(nullptr)
  , p0_g(nullptr)
  , gp(nullptr)
  , mu_g(nullptr)
  , vort(nullptr)
  , drag(nullptr)
  , xslopes_u(nullptr)
  , yslopes_u(nullptr)
  , zslopes_u(nullptr)
  , xslopes_s(nullptr)
  , yslopes_s(nullptr)
  , zslopes_s(nullptr)
  , diveu(nullptr)
  , mac_phi(nullptr)
  , u_mac(nullptr)
  , v_mac(nullptr)
  , w_mac(nullptr)
  , particle_cost(nullptr)
  , fluid_cost(nullptr)
{}

LevelData::~LevelData ()
{
  delete ep_g;
  delete ep_go;
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
  delete u_mac;
  delete v_mac;
  delete w_mac;
  delete particle_cost;
  delete fluid_cost;
}
