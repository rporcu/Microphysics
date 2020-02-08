#include <MFIX_LevelData.H>

using namespace amrex;

// Constructor
LevelData::LevelData (BoxArray const& ba,
                      DistributionMapping const& dm,
                      FabFactory<FArrayBox> const& fact,
                      const int nghost)
  : ep_g (ba, dm, 3, nghost, MFInfo(), fact)
  , ep_go (ba, dm, 3, nghost, MFInfo(), fact)
{}
