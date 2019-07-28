#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#include <mfix_proj_F.H>
#include <mfix_F.H>
#include <mfix.H>

#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MacProjector.H>

using namespace amrex;

//
// Computes the following decomposition:
// 
//    u + c*grad(phi)/ro = u*  with  div(ep*u) = 0
//
// Inputs:
// 
//   lev    = the AMR level
//   u,v,w  = the MAC velocity field to be projected
//   ep     = the cell-centered volume fraction
//   ro     = the cell-centered density
//
// Outputs:
//
//  phi     = the projection auxiliary function
//  u,v,w   = the PROJECTED MAC velocity field 
//
// Notes:
//
//  phi is computed by solving
//
//       div(ep*grad(phi)/ro) = div(ep * u*)
//
//  WARNING: this method returns the MAC velocity with up-to-date BCs in place
// 
void 
mfix::apply_MAC_projection (Vector< std::unique_ptr<MultiFab> >& u, 
                            Vector< std::unique_ptr<MultiFab> >& v,
                            Vector< std::unique_ptr<MultiFab> >& w,
                            Vector< std::unique_ptr<MultiFab> >& ep,
                            const Vector< std::unique_ptr<MultiFab> >& ro,
                            Real time, int steady_state)
{
   BL_PROFILE("mfix::apply_MAC_projection()");

   if (m_verbose)
      Print() << "MAC Projection:\n";

   // Check that everything is consistent with amrcore
   // update_internals();

   // Setup for solve
   Vector<Array<MultiFab*,AMREX_SPACEDIM> > vel;
   vel.resize(finest_level+1);

   if (m_verbose)
      Print() << " >> Before projection\n" ; 

   // ep_face and ro_face are now temporaries, no need to keep them outside this routine
   Vector< Array< std::unique_ptr<MultiFab>, AMREX_SPACEDIM> > ep_face;
   Vector< Array< std::unique_ptr<MultiFab>, AMREX_SPACEDIM> > ro_face;

   ep_face.resize(finest_level+1);
   ro_face.resize(finest_level+1);
   
   for ( int lev=0; lev <= finest_level; ++lev )
   {
      BoxArray x_edge_ba = grids[lev];
      x_edge_ba.surroundingNodes(0);
      BoxArray y_edge_ba = grids[lev];
      y_edge_ba.surroundingNodes(1);
      BoxArray z_edge_ba = grids[lev];
      z_edge_ba.surroundingNodes(2);

      ep_face[lev][0].reset(new MultiFab(x_edge_ba,dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      ep_face[lev][1].reset(new MultiFab(y_edge_ba,dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      ep_face[lev][2].reset(new MultiFab(z_edge_ba,dmap[lev],1,0,MFInfo(),*ebfactory[lev]));

      ro_face[lev][0].reset(new MultiFab(x_edge_ba,dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      ro_face[lev][1].reset(new MultiFab(y_edge_ba,dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      ro_face[lev][2].reset(new MultiFab(z_edge_ba,dmap[lev],1,0,MFInfo(),*ebfactory[lev]));

      // Compute ep at faces
      ep[lev]->FillBoundary(geom[lev].periodicity());
     
      average_cellcenter_to_face( GetArrOfPtrs(ep_face[lev]), *ep[lev], geom[lev] );
      average_cellcenter_to_face( GetArrOfPtrs(ro_face[lev]), *ro[lev], geom[lev] );

      MultiFab::Copy( *bcoeff_cc[lev][0], *(ep_face[lev][0]), 0, 0, 1, 0 );
      MultiFab::Copy( *bcoeff_cc[lev][1], *(ep_face[lev][1]), 0, 0, 1, 0 );
      MultiFab::Copy( *bcoeff_cc[lev][2], *(ep_face[lev][2]), 0, 0, 1, 0 );

      // Set velocity bcs -- before we multiply by ep
      set_MC_velocity_bcs( lev, u, v, w, time );

      // Compute ep*u at faces and store it in u, v, w
      MultiFab::Multiply( *u[lev], *(ep_face[lev][0]), 0, 0, 1, 0 );
      MultiFab::Multiply( *v[lev], *(ep_face[lev][1]), 0, 0, 1, 0 );
      MultiFab::Multiply( *w[lev], *(ep_face[lev][2]), 0, 0, 1, 0 );

      // Compute beta coefficients for div(beta*grad(phi)) = RHS:  beta = ep / ro
      MultiFab::Divide( *bcoeff_cc[lev][0], *(ro_face[lev][0]), 0, 0, 1, 0 );
      MultiFab::Divide( *bcoeff_cc[lev][1], *(ro_face[lev][1]), 0, 0, 1, 0 );
      MultiFab::Divide( *bcoeff_cc[lev][2], *(ro_face[lev][2]), 0, 0, 1, 0 );

      // Store in temporaries
      (vel[lev])[0] = u[lev].get();
      (vel[lev])[1] = v[lev].get();
      (vel[lev])[2] = w[lev].get();

      for (int i=0; i<3; ++i)
         (vel[lev])[i]->FillBoundary( geom[lev].periodicity() );
      
      if (m_verbose)
      {
         EB_computeDivergence(*mac_rhs[lev],
                              GetArrOfConstPtrs(vel[lev]),
                              geom[lev]);

         Print() << "  * On level "<< lev
                 << " max(abs(diveu)) = " << mfix_norm0(mac_rhs,lev,0) << "\n";
      }  
   }

   //
   // If we want to set max_coarsening_level we have to send it in to the constructor
   //
   LPInfo lp_info;
   lp_info.setMaxCoarseningLevel(mac_mg_max_coarsening_level);

   //
   // Perform MAC projection
   //
   MacProjector macproj( vel, GetVecOfArrOfPtrsConst(bcoeff_cc), geom, lp_info);

   macproj.setDomainBC  ( ppe_lobc, ppe_hibc );
   macproj.setVerbose   ( mac_mg_verbose);
   macproj.setCGVerbose ( mac_mg_cg_verbose);
   macproj.setMaxIter   ( mac_mg_maxiter);
   macproj.setCGMaxIter ( mac_mg_cg_maxiter);   
   // The default bottom solver is BiCG
   // Other options include:
   ///   Hypre IJ AMG solver
   //    macproj.getMLMG().setBottomSolver(MLMG::BottomSolver::hypre);
   ///   regular smoothing
   //    macproj.getMLMG().setBottomSolver(MLMG::BottomSolver::smoother);

   if (mac_bottom_solver_type == "smoother")
   {
      macproj.setBottomSolver(MLMG::BottomSolver::smoother);
   }
   else if (mac_bottom_solver_type == "cg")
   {
      macproj.setBottomSolver(MLMG::BottomSolver::cg);
   }
   else if (mac_bottom_solver_type == "bicgcg")
   {
      macproj.setBottomSolver(MLMG::BottomSolver::bicgcg);
   }
   else if (mac_bottom_solver_type == "cgbicg")
   {
      macproj.setBottomSolver(MLMG::BottomSolver::cgbicg);
   }
   else if (mac_bottom_solver_type == "hypre")
   {
      macproj.setBottomSolver(MLMG::BottomSolver::hypre);
   }

   if (steady_state)
   {
       // Solve using mac_phi as an initial guess -- note that mac_phi is
       //       stored from iteration to iteration
       macproj.project(GetVecOfPtrs(mac_phi), mac_mg_rtol,mac_mg_atol);
   } 
   else 
   {
       // Solve with initial guess of zero
       macproj.project(mac_mg_rtol,mac_mg_atol);
   }

   // Get MAC velocities at face CENTER by dividing solution by ep at faces
   if (m_verbose)
      Print() << " >> After projection\n" ; 
   
   for ( int lev=0; lev <= finest_level ; ++lev )
   {   
      if (m_verbose)
      {
         vel[lev][0]->FillBoundary( geom[lev].periodicity() );
         vel[lev][1]->FillBoundary( geom[lev].periodicity() );
         vel[lev][2]->FillBoundary( geom[lev].periodicity() );
         
         EB_computeDivergence(*mac_rhs[lev],
                              GetArrOfConstPtrs(vel[lev]),
                              geom[lev]);

         Print() << "  * On level "<< lev
                 << " max(abs(diveu)) = " << mfix_norm0(mac_rhs,lev,0) << "\n";
      } 

      // Now convert (eps u, eps v, eps w) back to u,v,w
      MultiFab::Divide( *u[lev], *(ep_face[lev][0]), 0, 0, 1, 0 );
      MultiFab::Divide( *v[lev], *(ep_face[lev][1]), 0, 0, 1, 0 );
      MultiFab::Divide( *w[lev], *(ep_face[lev][2]), 0, 0, 1, 0 ); 

      // Set velocity bcs
      set_MC_velocity_bcs( lev, u, v, w, time );
   }
}



//
// Set the BCs for velocity only
// 
void
mfix::set_MC_velocity_bcs ( int lev,
                             Vector< std::unique_ptr<MultiFab> >& u,
                             Vector< std::unique_ptr<MultiFab> >& v,
                             Vector< std::unique_ptr<MultiFab> >& w,
                             amrex::Real time)
{
   BL_PROFILE("MacProjection::set_MAC_velocity_bcs()");

   u[lev] -> FillBoundary( geom[lev].periodicity() );
   v[lev] -> FillBoundary( geom[lev].periodicity() );
   w[lev] -> FillBoundary( geom[lev].periodicity() );
     
   Box domain(geom[lev].Domain()); 

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi((*mac_rhs[lev]), false); mfi.isValid(); ++mfi)
   {
      const Box& bx = (*mac_rhs[lev])[mfi].box();

      set_mac_velocity_bcs(&time, bx, &mfi, lev, u, v, w, domain);
   }
}

#if 0
//
// Norm 0 for EB Multifab
//
Real
MacProjection::norm0 (const Vector<std::unique_ptr<MultiFab>>& mf, int lev)
{
   MultiFab mf_tmp( mf[lev]->boxArray(), mf[lev]->DistributionMap(), mf[lev]->nComp(),
                    0,  MFInfo(), *(*m_ebfactory)[lev]);
  
   MultiFab::Copy( mf_tmp, *mf[lev], 0, 0, 1, 0 );
   EB_set_covered( mf_tmp, 0.0 );
  
   return mf_tmp.norm0(0);
}
#endif
