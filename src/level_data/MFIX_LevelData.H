#ifndef _MFIX_LEVELDATA_H_
#define _MFIX_LEVELDATA_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabFactory.H>
#include <AMReX_MultiFab.H>

class LevelData
{
  public:
    // Constructor
    LevelData ();

    // Destructor
    ~LevelData ();

    // Cell-centered MultiFabs
    // Void fraction
    amrex::MultiFab* ep_g;
    amrex::MultiFab* ep_go;
    
    // Gas pressure fraction
    amrex::MultiFab* p_g;
    amrex::MultiFab* p_go;

    // Gas density
    amrex::MultiFab* ro_g;
    amrex::MultiFab* ro_go;

    // Tracer in gas
    amrex::MultiFab* trac;
    amrex::MultiFab* trac_o;

    // Gas velocity
    amrex::MultiFab* vel_g;
    amrex::MultiFab* vel_go;

};

#endif
