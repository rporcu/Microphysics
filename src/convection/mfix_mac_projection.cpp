#include <mfix.H>
#include <mfix_bc_parms.H>
#include <mfix_mf_helpers.H>

#include <AMReX_MacProjector.H>

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
//       div(ep*grad(phi)/ro) = div(ep u*) - S
//
//  This method returns the MAC velocity with up-to-date BCs in place
//
void
mfix::apply_MAC_projection (const bool update_laplacians,
                            Vector< MultiFab* > const& lap_T,
                            Vector< MultiFab* > const& lap_X,
                            Vector< MultiFab* > const& ep_u_mac,
                            Vector< MultiFab* > const& ep_v_mac,
                            Vector< MultiFab* > const& ep_w_mac,
                            Vector< MultiFab* > const& ep_g_in,
                            Vector< MultiFab* > const& ro_g_in,
                            Vector< MultiFab* > const& MW_g_in,
                            Vector< MultiFab* > const& T_g_in,
                            Vector< MultiFab* > const& cp_g_in,
                            Vector< MultiFab* > const& k_g_in,
                            Vector< MultiFab* > const& T_g_on_eb_in,
                            Vector< MultiFab* > const& k_g_on_eb_in,
                            Vector< MultiFab* > const& X_gk_in,
                            Vector< MultiFab* > const& D_gk_in,
                            Vector< MultiFab* > const& h_gk_in,
                            Vector< MultiFab* > const& txfr_in,
                            Vector< MultiFab* > const& ro_gk_txfr_in,
                            Real time)
{
  BL_PROFILE("mfix::apply_MAC_projection()");

  if (m_verbose)
    Print() << "MAC Projection:\n";

  // Check that everything is consistent with amrcore
  // update_internals();

  // Setup for solve
  Vector< Array<MultiFab*,3> > vel;
  vel.resize(finest_level+1);

  if (m_verbose)
    Print() << " >> Before projection\n" ;

  // Set bc's on density and ep_g so ro_face and ep_face will have correct values
  mfix_set_density_bcs(time, ro_g_in);

  // ro_face and ep_face are temporary, no need to keep it outside this routine
  Vector< Array<MultiFab*,3> > ro_face;
  Vector< Array<MultiFab*,3> > ep_face;

  ep_face.resize(finest_level+1);
  ro_face.resize(finest_level+1);

  for ( int lev=0; lev <= finest_level; ++lev )
  {
    ep_g_in[lev]->FillBoundary(geom[lev].periodicity());
    ro_g_in[lev]->FillBoundary(geom[lev].periodicity());

    ep_face[lev][0] = new MultiFab(ep_u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);
    ep_face[lev][1] = new MultiFab(ep_v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);
    ep_face[lev][2] = new MultiFab(ep_w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);

    ro_face[lev][0] = new MultiFab(ep_u_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);
    ro_face[lev][1] = new MultiFab(ep_v_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);
    ro_face[lev][2] = new MultiFab(ep_w_mac[lev]->boxArray(),dmap[lev],1,0,MFInfo(),*ebfactory[lev]);

    // Define ep and rho on face centroids (using interpolation from cell centroids)
    // The only use of bcs in this call is to test on whether a domain boundary is ext_dir
    // average_cellcenter_to_face(ro_face[lev], *ro_g_in[lev], geom[lev]);
    // average_cellcenter_to_face(ep_face[lev], *ep_g_in[lev], geom[lev]);
    EB_interp_CellCentroid_to_FaceCentroid (*ro_g_in[lev], ro_face[lev], 0, 0, 1, geom[lev], bcs_s);
    EB_interp_CellCentroid_to_FaceCentroid (*ep_g_in[lev], ep_face[lev], 0, 0, 1, geom[lev], bcs_s);

    // Compute ep_face into bcoeff
    MultiFab::Copy(*bcoeff[lev][0], *(ep_face[lev][0]), 0, 0, 1, 0);
    MultiFab::Copy(*bcoeff[lev][1], *(ep_face[lev][1]), 0, 0, 1, 0);
    MultiFab::Copy(*bcoeff[lev][2], *(ep_face[lev][2]), 0, 0, 1, 0);

    // Compute beta coefficients for div(beta*grad(phi)) = RHS:  beta = ep / ro
    MultiFab::Divide(*bcoeff[lev][0], *(ro_face[lev][0]), 0, 0, 1, 0);
    MultiFab::Divide(*bcoeff[lev][1], *(ro_face[lev][1]), 0, 0, 1, 0);
    MultiFab::Divide(*bcoeff[lev][2], *(ro_face[lev][2]), 0, 0, 1, 0);

    // Store (ep * u) in temporaries
    (vel[lev])[0] = ep_u_mac[lev];
    (vel[lev])[1] = ep_v_mac[lev];
    (vel[lev])[2] = ep_w_mac[lev];

    for (int i=0; i < 3; ++i)
      (vel[lev])[i]->FillBoundary(geom[lev].periodicity());

    if (m_verbose)
    {
      bool already_on_centroid = true;
      EB_computeDivergence(*m_leveldata[lev]->mac_rhs, GetArrOfConstPtrs(vel[lev]),
          geom[lev], already_on_centroid);

      Print() << "  * On level "<< lev << " max(abs(diveu)) = "
              << m_leveldata[lev]->mac_rhs->norm0(0,0,false,true) << "\n";
    }
  }

  for ( int lev=0; lev <= finest_level; ++lev )
  {
    delete ep_face[lev][0];
    delete ep_face[lev][1];
    delete ep_face[lev][2];

    delete ro_face[lev][0];
    delete ro_face[lev][1];
    delete ro_face[lev][2];
  }

  // Set the incompressibility or open-system constraint
  Vector< MultiFab* > depdt(finest_level+1);
  Vector< MultiFab* > mac_rhs(finest_level+1);
  Vector< MultiFab* > S_cc(finest_level+1);

  for (int lev(0); lev <= finest_level; ++lev) {
    depdt[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();
    mac_rhs[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();
    S_cc[lev] = MFHelpers::createFrom(*m_leveldata[lev]->ep_g, 0.0, 1).release();
  }

  if (open_system_constraint) {
    mfix_open_system_rhs(mac_rhs, update_laplacians, lap_T, lap_X, ep_g_in,
        ro_g_in, MW_g_in, T_g_in, cp_g_in, k_g_in, T_g_on_eb_in, k_g_on_eb_in,
        X_gk_in, D_gk_in, h_gk_in, txfr_in, ro_gk_txfr_in);
  }

  for (int lev(0); lev <= finest_level; ++lev) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*S_cc[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      Array4< Real > const& depdt_array = depdt[lev]->array(mfi);
      Array4< Real > const& mac_rhs_array = mac_rhs[lev]->array(mfi);
      Array4< Real > const& S_cc_array = S_cc[lev]->array(mfi);

      amrex::ParallelFor ( bx, [S_cc_array,depdt_array,mac_rhs_array]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        S_cc_array(i,j,k) = mac_rhs_array(i,j,k) - depdt_array(i,j,k);
      });
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
  Vector< Array <MultiFab const*, 3>> const_bcoeff;
  const_bcoeff.reserve(bcoeff.size());
  for (const auto& x : bcoeff) const_bcoeff.push_back(GetArrOfConstPtrs(x));

  Vector< const MultiFab* > const_S_cc = GetVecOfConstPtrs(S_cc);

  MacProjector macproj(vel,           MLMG::Location::FaceCentroid,  // location of vel 
                       const_bcoeff,  MLMG::Location::FaceCentroid,  // location of beta 
                                      MLMG::Location::CellCenter,    // location of phi 
                       geom, lp_info,
                       const_S_cc   , MLMG::Location::CellCentroid); // location of MAC rhs

  macproj.setDomainBC(BC::ppe_lobc, BC::ppe_hibc);

  if (steady_state)
  {
    // Solve using mac_phi as an initial guess -- note that mac_phi is
    //       stored from iteration to iteration
    macproj.project(get_mac_phi(), mac_mg_rtol, mac_mg_atol);
  }
  else
  {
    // Solve with initial guess of zero
    macproj.project(mac_mg_rtol, mac_mg_atol);
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
      EB_computeDivergence(*m_leveldata[lev]->mac_rhs, GetArrOfConstPtrs(vel[lev]),
                           geom[lev], already_on_centroid);

      Print() << "  * On level "<< lev << " max(abs(diveu)) = "
              << m_leveldata[lev]->mac_rhs->norm0(0,0,false,true) << "\n";
    }

    // Set bcs on (ep * u_mac)
    set_MAC_velocity_bcs(lev, ep_u_mac, ep_v_mac, ep_w_mac, time);

    ep_u_mac[lev]->FillBoundary(geom[lev].periodicity());
    ep_v_mac[lev]->FillBoundary(geom[lev].periodicity());
    ep_w_mac[lev]->FillBoundary(geom[lev].periodicity());
  }

  for (int lev(0); lev <= finest_level; lev++) {
    delete depdt[lev];
    delete mac_rhs[lev];
    delete S_cc[lev];
  }
}
