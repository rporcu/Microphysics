#ifndef MAC_PROJECTION_H_
#define MAC_PROJECTION_H_

#include <iostream>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_BLProfiler.H>


class MacProjection
{
public:
    
   MacProjection (amrex::AmrCore* a_amrcore, int a_nghost,
                  amrex::Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> >* a_ebfactory );

   ~MacProjection ();

   void set_bcs ( amrex::IArrayBox* a_bc_ilo, amrex::IArrayBox* a_bc_ihi,
                  amrex::IArrayBox* a_bc_jlo, amrex::IArrayBox* a_bc_jhi,
                  amrex::IArrayBox* a_bc_klo, amrex::IArrayBox* a_bc_khi );

   void apply_projection ( amrex::Vector< std::unique_ptr<amrex::MultiFab> >& u, 
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& v,
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& w,
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& ep,
                           const amrex::Vector< std::unique_ptr<amrex::MultiFab> >& ro,
                           const amrex::Real c );

   void apply_projection ( amrex::Vector< std::unique_ptr<amrex::MultiFab> >& u, 
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& v,
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& w,
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& ep,
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& epo,
                           const amrex::Vector< std::unique_ptr<amrex::MultiFab> >& ro,
                           const amrex::Real c, const amrex::Real dt );

   const amrex::MultiFab& get_phi ( int lev );
    
private:

   // Unit vectors in Cartesian space
   static amrex::IntVect e_x;
   static amrex::IntVect e_y;
   static amrex::IntVect e_z; 

   int m_nghost;

   amrex::AmrCore*  m_amrcore;
   amrex::Vector< std::unique_ptr<amrex::EBFArrayBoxFactory> >*  m_ebfactory;

   // BCs
   amrex::IArrayBox* m_bc_ilo;
   amrex::IArrayBox* m_bc_ihi;
   amrex::IArrayBox* m_bc_jlo;
   amrex::IArrayBox* m_bc_jhi;
   amrex::IArrayBox* m_bc_klo;
   amrex::IArrayBox* m_bc_khi;
    
   amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_diveu;
   amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_phi;
   amrex::Vector< amrex::Vector< std::unique_ptr<amrex::MultiFab> > > m_b;
   amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,AMREX_SPACEDIM> > m_gradphi;


   // 
   // Stuff for linear solver
   // 
   // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
   std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_lobc;
   std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_hibc;

   bool verbose = true;
    
   int  m_mg_verbose = 0;
   int  m_mg_cg_verbose = 0;
   int  m_mg_max_iter = 100;
   int  m_mg_cg_maxiter = 100;
   int  m_mg_max_fmg_iter = 0;
   amrex::Real m_mg_rtol = 1.0e-11;
   amrex::Real m_mg_atol = 1.0e-14;

   void read_inputs ();
    
   void update_internals ();

   void compute_diveu ( amrex::Vector< std::unique_ptr<amrex::MultiFab> >& u,
                        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& v,
                        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& w,
                        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& ep,
                        const amrex::Real c);

   void set_velocity_bcs ( int lev,
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& u,
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& v,
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& w );
    
   void set_ccmf_bcs ( int lev, amrex::MultiFab& mf );

   void compute_b_coeff ( const amrex::Vector< std::unique_ptr<amrex::MultiFab> >& u, 
                          const amrex::Vector< std::unique_ptr<amrex::MultiFab> >& v,
                          const amrex::Vector< std::unique_ptr<amrex::MultiFab> >& w,
                          const amrex::Vector< std::unique_ptr<amrex::MultiFab> >& ep,
                          const amrex::Vector< std::unique_ptr<amrex::MultiFab> >& ro );

   void solve_for_phi ();

   void project_velocity ( amrex::Vector< std::unique_ptr<amrex::MultiFab> >& u,
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& v,
                           amrex::Vector< std::unique_ptr<amrex::MultiFab> >& w,
                           const amrex::Vector< std::unique_ptr<amrex::MultiFab> >& ro,
                           const amrex::Real c );

   //
   // Stuff to compute norms of EB-MultiFabs
   //
   amrex::Real norm0 (const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mf, int lev);
   
};


#endif
