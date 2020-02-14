#ifndef _MFIX_LEVELDATA_H_
#define _MFIX_LEVELDATA_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabFactory.H>
#include <AMReX_MultiFab.H>

class LevelData
{
  public:
    // Default Constructor
    LevelData () = default;

    // Cell-centered MultiFabs
    // Void fraction
    amrex::MultiFab* ep_g;
    amrex::MultiFab* ep_go;
};

#endif
