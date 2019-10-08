#ifndef MFIX_MFHELPERS_H
#define MFIX_MFHELPERS_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EBMultiFabUtil.H>

namespace MFHelpers
{
    
    // Creates an exact copy ("shape" + value) of "mold"
    std::unique_ptr<amrex::MultiFab> createFrom ( amrex::MultiFab& mold );

    // Creates a copy of "mold" and initializes it to "val".
    std::unique_ptr<amrex::MultiFab> createFrom ( amrex::MultiFab& mold,
                                                  amrex::Real val );

    // Creates a copy of "mold" with nGrow ghosts and initializes it to "val"   
    std::unique_ptr<amrex::MultiFab> createFrom ( amrex::MultiFab& mold,
                                                  amrex::Real val,
                                                  int nGrow );
    
    // Copies src to dst including all components and ghost cells
    void copy (amrex::MultiFab& dst, amrex::MultiFab& src);

    // Copies src to dst including all components and nGrow ghost cells
    void copy (amrex::MultiFab& dst, amrex::MultiFab& src, int nGrow);

    // Compute norm 0 of component "comp" of EB MultiFab
    amrex::Real norm0 ( amrex::MultiFab& mf, int comp );
    
    // Compute norm 0 of component "comp" of EB MultiFab on level "lev"
    amrex::Real norm0 ( const amrex::Vector< std::unique_ptr< amrex::MultiFab> >& mf,
                        int lev,
                        int comp );
    
    // Compute max of component "comp" of EB MultiFab
    amrex::Real max ( amrex::MultiFab& mf, int comp );

    // Compute min of component "comp" of EB MultiFab
    amrex::Real min ( amrex::MultiFab& mf, int comp );

}



#endif 
