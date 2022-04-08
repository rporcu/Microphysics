#ifndef _MFIX_BCLIST_H_
#define _MFIX_BCLIST_H_

#include <AMReX_Vector.H>
#include <AMReX_IArrayBox.H>

// Boundary conditions list auxiliary structure
//
class BCList {

  public:
    // Constructor
    BCList (const int nlev_in,
            const int undefined = 0,
            const int ig = 9,
            const int pinf = 10,
            const int pout = 11,
            const int minf = 20,
            const int nsw = 100,
            const int eb  = 111);

    // Destructor
    ~BCList();

    int get_undefined () const {return und_;};
    int get_ig () const {return ig_;};

    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_pinf () const {return pinf_;};
    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_pout () const {return pout_;};
    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_minf () const {return minf_;};
    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_nsw  () const {return nsw_ ;};
    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_eb   () const {return eb_  ;};

    void MakeBCArrays (int nghost,
                       bool ooo_debug,
                       amrex::Vector<amrex::Geometry>& geom);

  private:
    const int nlev;
    const int und_;
    const int ig_;
    const int pinf_;
    const int pout_;
    const int minf_;
    const int nsw_;
    const int eb_;

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
