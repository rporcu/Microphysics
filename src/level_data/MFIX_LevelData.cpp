#include <MFIX_LevelData.H>

using namespace amrex;

LevelData::LevelData ()
  : ep_go(nullptr)
{}

LevelData::~LevelData ()
{
  delete ep_go;
}
