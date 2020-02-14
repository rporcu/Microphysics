#include <MFIX_LevelData.H>

using namespace amrex;

LevelData::LevelData ()
  : ep_g(new MultiFab())
  , ep_go(new MultiFab())
{}
