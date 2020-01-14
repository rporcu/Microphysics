#include <mfix.H>
#include <AMReX_MacProjector.H>
#include <MFIX_BC_Parms.H>

using namespace amrex;

//
// Computes the following decomposition:
//
//    u* = u + c*grad(phi)/ro with  div(ep u) = 0
//
// Notes:
//
//  phi is computed by solving
//
//       div(ep*grad(phi)/ro) = div(ep u*)
//
//  This method returns the MAC velocity with up-to-date BCs in place
//
void
mfix::apply_MAC_projection (Vector< std::unique_ptr<MultiFab> >& ep_u_mac,
                            Vector< std::unique_ptr<MultiFab> >& ep_v_mac,
                            Vector< std::unique_ptr<MultiFab> >& ep_w_mac,
                            Vector< std::unique_ptr<MultiFab> >& ep_in,
                            Vector< std::unique_ptr<MultiFab> >& ro_in,
                            Real time)
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

    // Set bc's on density and ep_g so ro_face and ep_face will have correct values
    mfix_set_density_bcs(time,ro_in);

   // ro_face and ep_face are temporary, no need to keep it outside this routine
   Vector< Array< std::unique_ptr<MultiFab>, AMREX_SPACEDIM> > ro_face;
   Vector< Array< std::unique_ptr<MultiFab>, AMREX_SPACEDIM> > ep_face;

   ep_face.resize(finest_level+1);
   ro_face.resize(finest_level+1);

   for ( int lev=0; lev <= finest_level; ++lev )
   {
      ep_in[lev]->FillBoundary(geom[lev].periodicity());
      ro_in[lev]->FillBoundary(geom[lev].periodicity());

      ep_face[lev][0].reset(new MultiFab(ep_u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      ep_face[lev][1].reset(new MultiFab(ep_v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      ep_face[lev][2].reset(new MultiFab(ep_w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));

      ro_face[lev][0].reset(new MultiFab(ep_u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      ro_face[lev][1].reset(new MultiFab(ep_v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));
      ro_face[lev][2].reset(new MultiFab(ep_w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]));

      // Define ep and rho on faces
      average_cellcenter_to_face(GetArrOfPtrs(ep_face[lev]), *ep_in[lev], geom[lev]);
      average_cellcenter_to_face(GetArrOfPtrs(ro_face[lev]), *ro_in[lev], geom[lev]);

      // Compute ep_face into bcoeff
      MultiFab::Copy(*bcoeff[lev][0], *(ep_face[lev][0]), 0, 0, 1, 0);
      MultiFab::Copy(*bcoeff[lev][1], *(ep_face[lev][1]), 0, 0, 1, 0);
      MultiFab::Copy(*bcoeff[lev][2], *(ep_face[lev][2]), 0, 0, 1, 0);

      // Compute beta coefficients for div(beta*grad(phi)) = RHS:  beta = ep / ro
      MultiFab::Divide(*bcoeff[lev][0], *(ro_face[lev][0]), 0, 0, 1, 0);
      MultiFab::Divide(*bcoeff[lev][1], *(ro_face[lev][1]), 0, 0, 1, 0);
      MultiFab::Divide(*bcoeff[lev][2], *(ro_face[lev][2]), 0, 0, 1, 0);

      // Store (ep * u) in temporaries
      (vel[lev])[0] = ep_u_mac[lev].get();
      (vel[lev])[1] = ep_v_mac[lev].get();
      (vel[lev])[2] = ep_w_mac[lev].get();

      for (int i=0; i<AMREX_SPACEDIM; ++i)
         (vel[lev])[i]->FillBoundary( geom[lev].periodicity() );

      if (m_verbose)
      {
         bool already_on_centroid = true;
         EB_computeDivergence(*mac_rhs[lev],
                              GetArrOfConstPtrs(vel[lev]),
                              geom[lev], already_on_centroid);

         Print() << "  * On level "<< lev
                 << " max(abs(diveu)) = "
                 << mac_rhs[lev]->norm0(0,0,false,true) << "\n";
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
   MacProjector macproj(vel, GetVecOfArrOfPtrsConst(bcoeff), geom, lp_info);

   macproj.setDomainBC(BC::ppe_lobc, BC::ppe_hibc);

   if (steady_state)
   {
       // Solve using mac_phi as an initial guess -- note that mac_phi is
       //       stored from iteration to iteration
       macproj.project(GetVecOfPtrs(mac_phi), mac_mg_rtol,mac_mg_atol,MLMG::Location::FaceCentroid);
   }
   else
   {
       // Solve with initial guess of zero
       macproj.project(mac_mg_rtol,mac_mg_atol,MLMG::Location::FaceCentroid);
   }

   // Get MAC velocities at face CENTER by dividing solution by ep at faces
   if (m_verbose)
      Print() << " >> After projection\n" ;

   for ( int lev=0; lev <= finest_level ; ++lev )
   {
      if (m_verbose)
      {
         vel[lev][0]->FillBoundary(geom[lev].periodicity());
         vel[lev][1]->FillBoundary(geom[lev].periodicity());
         vel[lev][2]->FillBoundary(geom[lev].periodicity());

         bool already_on_centroid = true;
         EB_computeDivergence(*mac_rhs[lev],
                              GetArrOfConstPtrs(vel[lev]),
                              geom[lev], already_on_centroid);

         Print() << "  * On level "<< lev
                 << " max(abs(diveu)) = "
                 << mac_rhs[lev]->norm0(0,0,false,true) << "\n";
      }

      // Set bcs on (ep * u_mac)
      set_MAC_velocity_bcs(lev, ep_u_mac, ep_v_mac, ep_w_mac, time);

      ep_u_mac[lev]->FillBoundary(geom[lev].periodicity());
      ep_v_mac[lev]->FillBoundary(geom[lev].periodicity());
      ep_w_mac[lev]->FillBoundary(geom[lev].periodicity());
   }
}
