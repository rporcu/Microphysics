#ifndef _MFIX_BC_LIST_H_
#define _MFIX_BC_LIST_H_

#include <AMReX_Vector.H>
#include <AMReX_IArrayBox.H>


// Boundary conditions list auxiliary structure
//
struct BCList {

  public:
    // Types
    enum {undefined=0, ig=9, pinf=10, pout=11, minf=20, nsw=100, eb=111};

    // Constructor
    BCList (const int nlev_in);

    // Destructor
    ~BCList();

    void MakeBCArrays (int nghost,
                       bool ooo_debug,
                       amrex::Vector<amrex::Geometry>& geom);

  private:
    const int nlev;

  public:
    // Boundary conditions types
    amrex::Vector<amrex::IArrayBox*> bc_ilo;
    amrex::Vector<amrex::IArrayBox*> bc_ihi;
    amrex::Vector<amrex::IArrayBox*> bc_jlo;
    amrex::Vector<amrex::IArrayBox*> bc_jhi;
    amrex::Vector<amrex::IArrayBox*> bc_klo;
    amrex::Vector<amrex::IArrayBox*> bc_khi;
};

#endif
