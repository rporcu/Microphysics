#include <MFIX_LevelData.H>

using namespace amrex;

LevelData::LevelData ()
  : ep_g(nullptr)
  , ep_go(nullptr)
{}

LevelData::~LevelData ()
{
  delete ep_g;
  delete ep_go;
}
