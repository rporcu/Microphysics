#ifndef MFIX_REGIONS_PARMS_H_
#define MFIX_REGIONS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

namespace REGIONS
{
  // Name and extents of regions

  void Initialize ();

  // Returns the index of region with matching name.
  // -1 is returned if the region is not matched.
  int getIndex(std::string input);
  const amrex::RealBox* getRegion(std::string name);

}

#endif
