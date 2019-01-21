#include <mfix.H>
#include <AMReX_EB_utils.H>

void
mfix::AllocateArrays (int lev)
{
    mfix_update_ebfactory(lev);

    // ********************************************************************************
    // Cell- or node-based arrays
    // ********************************************************************************

    // Void fraction
    ep_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    ep_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    ep_g[lev]->setVal(1.0);
    ep_go[lev]->setVal(1.0);

    // Gas density
    ro_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    ro_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    ro_g[lev]->setVal(0.);
    ro_go[lev]->setVal(0.);

    // Gas bulk density
    rop_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    rop_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    rop_g[lev]->setVal(0.);
    rop_go[lev]->setVal(0.);

    const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

    p0_g[lev].reset(new MultiFab(nd_grids,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    p0_g[lev]->setVal(0.);

    p_g[lev].reset(new MultiFab(nd_grids,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    p_g[lev]->setVal(0.);

    p_go[lev].reset(new MultiFab(nd_grids,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    p_go[lev]->setVal(0.);

    pp_g[lev].reset(new MultiFab(nd_grids,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    pp_g[lev]->setVal(0.);

    phi[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0, MFInfo(), *ebfactory[lev]));
    phi[lev]->setVal(0.);

    diveu[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0, MFInfo(), *ebfactory[lev]));
    diveu[lev]->setVal(0.);


    // Presssure gradients
    gp[lev].reset(new MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    gp0[lev].reset(new MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    gp[lev]->setVal(0.);
    gp0[lev]->setVal(0.);

    // Molecular viscosity
    mu_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    mu_g[lev]->setVal(0.);

    // Coefficient of grad(div(u)) in viscous terms
    lambda_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    lambda_g[lev]->setVal(0.);

    // Current velocity
    vel_g[lev].reset(new MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    vel_g[lev]->setVal(0.);

    // Old velocity
    vel_go[lev].reset(new  MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    vel_go[lev]->setVal(0.);

    // Div(u)
    trD_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    trD_g[lev]->setVal(0.);

    // Vorticity
    vort[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    vort[lev]->setVal(0.);

    // This is the deposition onto the grid of the beta coefficient
    // for fluid vel in the expression beta*(fluid_vel _ particle_vel)
    // Note this only needs one component since all velocity components are co-located
    f_gds[lev].reset(new  MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    f_gds[lev]->setVal(0.);

    // This is the deposition onto the grid of the drag term experienced by the particle
    drag[lev].reset(new  MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    drag[lev]->setVal(0.);

    // Arrays to store the solution and rhs for the diffusion solve
    phi_diff[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    rhs_diff[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));

    // Slopes in x-direction
    xslopes[lev].reset(new  MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    xslopes[lev] -> setVal(0.);

    // Slopes in y-direction
    yslopes[lev].reset(new  MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    yslopes[lev] -> setVal(0.);

    // Slopes in z-direction
    zslopes[lev].reset(new  MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    zslopes[lev] -> setVal(0.);

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

    // When the pressure is on nodes, bcoeff is at cell centers
    bcoeff[lev][0].reset(new  MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    bcoeff[lev][0]->setVal(0.);

    bcoeff[lev][1].reset(new  MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    bcoeff[lev][1]->setVal(0.);

    bcoeff[lev][2].reset(new  MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    bcoeff[lev][2]->setVal(0.);

    // ****************************************************************

    // Create a BoxArray on x-faces.
    BoxArray x_edge_ba = grids[lev];
    x_edge_ba.surroundingNodes(0);
    bcoeff_diff[lev][0].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    m_u_mac[lev].reset(new MultiFab(x_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));

    // Create a BoxArray on y-faces.
    BoxArray y_edge_ba = grids[lev];
    y_edge_ba.surroundingNodes(1);
    bcoeff_diff[lev][1].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    m_v_mac[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));

    // Create a BoxArray on y-faces.
    BoxArray z_edge_ba = grids[lev];
    z_edge_ba.surroundingNodes(2);
    bcoeff_diff[lev][2].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    m_w_mac[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));

    bcoeff_diff[lev][0]->setVal(0.);
    bcoeff_diff[lev][1]->setVal(0.);
    bcoeff_diff[lev][2]->setVal(0.);

    m_u_mac[lev]->setVal(0.);
    m_v_mac[lev]->setVal(0.);
    m_w_mac[lev]->setVal(0.);
}


void
mfix::RegridArrays (int lev)
{
    bool need_regrid = mfix_update_ebfactory(lev);

    // exit this function is ebfactory has not been updated because that means
    // that dm and ba haven't changed
    if (!need_regrid)
        return;

    // ********************************************************************************
    // Cell-based arrays
    // ********************************************************************************
    //
    // After calling copy() with dst_ngrow set to ng, we do not need to call
    // FillBoundary().
    //
    //

    // Void fraction
    int ng = ep_g[lev]->nGrow();
    std::unique_ptr<MultiFab> ep_g_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    ep_g_new->setVal(1.0);
    ep_g_new->copy(*ep_g[lev],0,0,1,0,ng);
    ep_g[lev] = std::move(ep_g_new);

    // Old void fraction
    ng = ep_go[lev]->nGrow();
    std::unique_ptr<MultiFab> ep_go_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev] ));
    ep_go_new->setVal(1.0);
    ep_go_new->copy(*ep_go[lev],0,0,1,0,ng);
    ep_go[lev] = std::move(ep_go_new);

    // Gas density
    ng = ro_g[lev]->nGrow();
    std::unique_ptr<MultiFab> ro_g_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    ro_g_new->setVal(0.0);
    ro_g_new->copy(*ro_g[lev],0,0,1,0,ng);
    ro_g[lev] = std::move(ro_g_new);

    // Old gas density
    ng = ro_go[lev]->nGrow();
    std::unique_ptr<MultiFab> ro_go_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    ro_go_new->setVal(0.0);
    ro_go_new->copy(*ro_go[lev],0,0,1,ng,ng);
    ro_go[lev] = std::move(ro_go_new);

    // Gas bulk density
    ng = rop_g[lev]->nGrow();
    std::unique_ptr<MultiFab> rop_g_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    rop_g_new->setVal(0.0);
    rop_g_new->copy(*rop_g[lev],0,0,1,0,ng);
    rop_g[lev] = std::move(rop_g_new);

    // Old gas bulk density
    ng = rop_go[lev]->nGrow();
    std::unique_ptr<MultiFab> rop_go_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    rop_go_new->setVal(0.0);
    rop_go_new->copy(*rop_go[lev],0,0,1,0,ng);
    rop_go[lev] = std::move(rop_go_new);

    const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

    ng = p_g[lev]->nGrow();
    std::unique_ptr<MultiFab> p_g_new(new MultiFab(nd_grids,dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    p_g_new->setVal(0.0);
    p_g_new->copy(*p_g[lev],0,0,1,0,ng);
    p_g[lev] = std::move(p_g_new);

    ng = p_go[lev]->nGrow();
    std::unique_ptr<MultiFab> p_go_new(new MultiFab(nd_grids,dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    p_go_new->setVal(0.0);
    p_go_new->copy(*p_go[lev],0,0,1,0,ng);
    p_go[lev] = std::move(p_go_new);

    ng = p0_g[lev]->nGrow();
    std::unique_ptr<MultiFab> p0_g_new(new MultiFab(nd_grids,dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    p0_g_new->setVal(0.0);
    p0_g_new->copy(*p0_g[lev],0,0,1,0,ng);
    p0_g[lev] = std::move(p0_g_new);

    ng = pp_g[lev]->nGrow();
    std::unique_ptr<MultiFab> pp_g_new(new MultiFab(nd_grids,dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    pp_g_new->FillBoundary(p0_periodicity);
    pp_g_new->copy(*pp_g[lev],0,0,1,0,ng);
    pp_g[lev] = std::move(pp_g_new);

    ng = diveu[lev]->nGrow();
    std::unique_ptr<MultiFab> diveu_new(new MultiFab(nd_grids,dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    diveu[lev] = std::move(diveu_new);
    diveu[lev]->setVal(0.);

    ng = phi[lev] -> nGrow();
    std::unique_ptr<MultiFab> phi_new(new MultiFab(nd_grids,dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    phi[lev] = std::move(phi_new);
    phi[lev]->setVal(0.);

    std::unique_ptr<MultiFab> bc0_new(new MultiFab(grids[lev],dmap[lev],1,bcoeff[lev][0]->nGrow(),
                                                  MFInfo(), *ebfactory[lev] ));
    bcoeff[lev][0] = std::move(bc0_new);
    bcoeff[lev][0]->setVal(0.);

    std::unique_ptr<MultiFab> bc1_new(new MultiFab(grids[lev],dmap[lev],1,bcoeff[lev][1]->nGrow(),
                                                  MFInfo(), *ebfactory[lev] ));
    bcoeff[lev][1] = std::move(bc1_new);
    bcoeff[lev][1]->setVal(0.);

    std::unique_ptr<MultiFab> bc2_new(new MultiFab(grids[lev],dmap[lev],1,bcoeff[lev][2]->nGrow(),
                                                  MFInfo(), *ebfactory[lev] ));
    bcoeff[lev][2] = std::move(bc2_new);
    bcoeff[lev][2]->setVal(0.);

    // Molecular viscosity
    ng = mu_g[lev]->nGrow();
    std::unique_ptr<MultiFab> mu_g_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    mu_g_new->setVal(0.0);
    mu_g_new->copy(*mu_g[lev],0,0,1,0,ng);
    mu_g[lev] = std::move(mu_g_new);

    // Lambda
    ng = lambda_g[lev]->nGrow();
    std::unique_ptr<MultiFab> lambda_g_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    lambda_g_new->setVal(0.0);
    lambda_g_new->copy(*lambda_g[lev],0,0,1,0,ng);
    lambda_g[lev] = std::move(lambda_g_new);

    // Gas velocity
    ng = vel_g[lev]->nGrow();
    std::unique_ptr<MultiFab> vel_g_new(new MultiFab(grids[lev],dmap[lev],vel_g[lev]->nComp(),ng,
                                                     MFInfo(), *ebfactory[lev]));
    vel_g_new->setVal(0.0);
    vel_g_new->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),0,ng);
    vel_g[lev] = std::move(vel_g_new);

    // Old gas velocity
    ng = vel_go[lev]->nGrow();
    std::unique_ptr<MultiFab> vel_go_new(new MultiFab(grids[lev],dmap[lev],vel_go[lev]->nComp(),ng,
                                                      MFInfo(), *ebfactory[lev]));
    vel_go_new->setVal(0.0);
    vel_go_new->copy(*vel_go[lev],0,0,vel_go[lev]->nComp(),0,ng);
    vel_go[lev] = std::move(vel_go_new);

    // Pressure gradients
    ng = gp[lev]->nGrow();
    std::unique_ptr<MultiFab> gp_new(new MultiFab(grids[lev],dmap[lev],3,ng, MFInfo(), *ebfactory[lev]));
    gp_new->setVal(0.0);
    gp_new->copy(*gp[lev],0,0,1,0,ng);
    gp[lev] = std::move(gp_new);

    // Pressure gradients
    ng = gp0[lev]->nGrow();
    std::unique_ptr<MultiFab> gp0_new(new MultiFab(grids[lev],dmap[lev],3,ng, MFInfo(), *ebfactory[lev]));
    gp0_new->setVal(0.0);
    gp0_new->copy(*gp0[lev],0,0,1,0,ng);
    gp0[lev] = std::move(gp0_new);

    // Trace(D)
    ng = trD_g[lev]->nGrow();
    std::unique_ptr<MultiFab> trD_g_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    trD_g[lev] = std::move(trD_g_new);
    trD_g[lev]->setVal(0.);

    // Vorticity
    ng = vort[lev]->nGrow();
    std::unique_ptr<MultiFab> vort_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    vort[lev] = std::move(vort_new);
    vort[lev]->setVal(0.);

    // Coefficient in drag
    ng = f_gds[lev]->nGrow();
    std::unique_ptr<MultiFab> f_gds_new(new MultiFab(grids[lev],dmap[lev],f_gds[lev]->nComp(),ng,
                                                     MFInfo(), *ebfactory[lev]));
    f_gds[lev] = std::move(f_gds_new);
    f_gds[lev]->setVal(0.);

    // Particle/fluid drag
    ng = drag[lev]->nGrow();
    std::unique_ptr<MultiFab> drag_new(new MultiFab(grids[lev],dmap[lev],drag[lev]->nComp(),ng,
                                                    MFInfo(), *ebfactory[lev]));
    drag[lev] = std::move(drag_new);
    drag[lev]->setVal(0.);

    // Arrays to store the solution and rhs for the diffusion solve
    std::unique_ptr<MultiFab> phi_diff_new(new  MultiFab(grids[lev], dmap[lev], 1, nghost,
                                                         MFInfo(), *ebfactory[lev]));
    phi_diff[lev] = std::move(phi_diff_new);
    phi_diff[lev] -> setVal(0.);

    std::unique_ptr<MultiFab> rhs_diff_new(new  MultiFab(grids[lev], dmap[lev], 1, nghost,
                                                         MFInfo(), *ebfactory[lev]));
    rhs_diff[lev] = std::move(rhs_diff_new);
    rhs_diff[lev] -> setVal(0.);

    // Slopes in x-direction
    ng = xslopes[lev] -> nGrow();
    std::unique_ptr<MultiFab> xslopes_new(new  MultiFab(grids[lev], dmap[lev],xslopes[lev]->nComp(),nghost,
                                                        MFInfo(), *ebfactory[lev]));
    xslopes[lev] = std::move(xslopes_new);
    xslopes[lev] -> setVal(0.);

    // Slopes in y-direction
    ng = yslopes[lev] -> nGrow();
    std::unique_ptr<MultiFab> yslopes_new(new  MultiFab(grids[lev], dmap[lev],yslopes[lev]->nComp(),nghost,
                                                        MFInfo(), *ebfactory[lev]));
    yslopes[lev] = std::move(yslopes_new);
    yslopes[lev] -> setVal(0.);

    // Slopes in z-direction
    ng = zslopes[lev] -> nGrow();
    std::unique_ptr<MultiFab> zslopes_new(new  MultiFab(grids[lev], dmap[lev],zslopes[lev]->nComp(),nghost,
                                                        MFInfo(), *ebfactory[lev]));
    zslopes[lev] = std::move(zslopes_new);
    zslopes[lev] -> setVal(0.);

   /****************************************************************************
    * Face-based Arrays                                                        *
    ****************************************************************************/

    BoxArray x_ba = grids[lev];
    x_ba = x_ba.surroundingNodes(0);

    // MAC velocity
    std::unique_ptr<MultiFab> u_mac_new(new MultiFab(x_ba,dmap[lev],1,nghost,MFInfo(), *ebfactory[lev]));
    m_u_mac[lev] = std::move(u_mac_new);
    m_u_mac[lev] -> setVal(0.0);

    // Diffusion coefficient on x-faces
    std::unique_ptr<MultiFab> bc0_diff_new(new MultiFab(x_ba,dmap[lev],1,nghost,MFInfo(), *ebfactory[lev]));
    bcoeff_diff[lev][0] = std::move(bc0_diff_new);
    bcoeff_diff[lev][0] -> setVal(0.0);

   //****************************************************************************

    BoxArray y_ba = grids[lev];
    y_ba = y_ba.surroundingNodes(1);

    // MAC velocity
    std::unique_ptr<MultiFab> v_mac_new(new MultiFab(y_ba,dmap[lev],1,nghost,MFInfo(), *ebfactory[lev]));
    m_v_mac[lev] = std::move(v_mac_new);
    m_v_mac[lev] -> setVal(0.0);

    // Diffusion coefficient on y-faces
    std::unique_ptr<MultiFab> bc1_diff_new(new MultiFab(y_ba,dmap[lev],1,nghost,MFInfo(), *ebfactory[lev]));
    bcoeff_diff[lev][1] = std::move(bc1_diff_new);
    bcoeff_diff[lev][1] -> setVal(0.0);

   //****************************************************************************

    BoxArray z_ba = grids[lev];
    z_ba = z_ba.surroundingNodes(2);

    // MAC velocity
    std::unique_ptr<MultiFab> w_mac_new(new MultiFab(z_ba,dmap[lev],1,nghost,MFInfo(), *ebfactory[lev]));
    m_w_mac[lev] = std::move(w_mac_new);
    m_w_mac[lev] -> setVal(0.0);

    // Diffusion coefficient on z-faces
    std::unique_ptr<MultiFab> bc2_diff_new(new MultiFab(z_ba,dmap[lev],1,nghost,MFInfo(), *ebfactory[lev]));
    bcoeff_diff[lev][2] = std::move(bc2_diff_new);
    bcoeff_diff[lev][2] -> setVal(0.0);

    // ********************************************************************************
    // Make sure we fill the ghost cells as appropriate -- this is copied from init_fluid
    // ********************************************************************************

    fill_mf_bc(lev,*ep_g[lev]);
    fill_mf_bc(lev,*ep_go[lev]);
    fill_mf_bc(lev,*ro_g[lev]);
    fill_mf_bc(lev,*ro_go[lev]);
    fill_mf_bc(lev,*rop_g[lev]);
    fill_mf_bc(lev,*rop_go[lev]);

    fill_mf_bc(lev,*mu_g[lev]);
    fill_mf_bc(lev,*lambda_g[lev]);

}


//! This function regrids the level set function and updates
//! `particle_ebfactory`. This has to be done separately from the regridding of
//! the other field variables since LS and `particle_ebfactory` "live" on the
//! particle grids. Furthermore, the LS always has an additional "refined" level
//! if not operating in multi-level mode. Hence slightly different regridding
//! rules are needed.
void
mfix::RegridLevelSetArray (int a_lev)
{
   // First check if particle_ebfactory is allocated with the proper dm and ba

   // This assert is to verify that some kind of EB geometry has already been
   // defined
   AMREX_ASSERT(not EB2::IndexSpace::empty());

   const DistributionMapping&      dm = pc->ParticleDistributionMap(a_lev);
   const BoxArray&                 ba = pc->ParticleBoxArray(a_lev);

   int changed = false;

   if ( particle_ebfactory[a_lev].get() == nullptr )
   {
      amrex::Print() << "Updating particle ebfactory 1" << std::endl;

      particle_ebfactory[a_lev].reset(
          new EBFArrayBoxFactory(* eb_levels[a_lev], geom[a_lev], ba, dm,
                                 {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                  m_eb_full_grow_cells}, m_eb_support_level)
          );

      changed = true;

   }
   else
   {
      amrex::Print() << "Updating particle ebfactory 2" << std::endl;

      const DistributionMapping&  eb_dm = particle_ebfactory[a_lev]->DistributionMap();
      const BoxArray&             eb_ba = particle_ebfactory[a_lev]->boxArray();

      if ( (dm != eb_dm) || (ba != eb_ba) )
      {
          particle_ebfactory[a_lev].reset(
              new EBFArrayBoxFactory( * eb_levels[a_lev], geom[a_lev], ba, dm,
                                      {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                       m_eb_full_grow_cells}, m_eb_support_level)
              );

         changed = true;
      }
   }

   if (changed)
   {

       Print() << "Regridding level-set on lev = " << a_lev << std::endl;

       const BoxArray nd_ba = amrex::convert(grids[a_lev], IntVect::TheNodeVector());

       std::unique_ptr<MultiFab> new_level_set(new MultiFab);
       MFUtil::regrid(* new_level_set, nd_ba, dmap[a_lev],
                      * particle_ebfactory[a_lev], * level_sets[a_lev], true);
       level_sets[a_lev] = std::move(new_level_set);

       std::unique_ptr<MultiFab> new_impfunc(new MultiFab);
       MFUtil::regrid(* new_impfunc, nd_ba, dmap[a_lev],
                      * particle_ebfactory[a_lev], * implicit_functions[a_lev], true);

       implicit_functions[a_lev] = std::move(new_impfunc);

       //________________________________________________________________________
       // If we're operating in single-level mode, the level-set has a second
       // (refined) MultiFab that also needs to be regridded.
       if ((nlev == 1) && (a_lev == 0))
       {
           Print() << "Also regridding refined level-set" << std::endl;

           BoxArray ref_nd_ba = amrex::convert(grids[a_lev], IntVect::TheNodeVector());
           ref_nd_ba.refine(levelset__refinement);

           std::unique_ptr<MultiFab> new_level_set(new MultiFab);
           MFUtil::regrid(* new_level_set, ref_nd_ba, dmap[a_lev],
                          * level_sets[a_lev + 1], true);
           level_sets[a_lev + 1] = std::move(new_level_set);

           std::unique_ptr<MultiFab> new_impfunc(new MultiFab);
           MFUtil::regrid(* new_impfunc, ref_nd_ba, dmap[a_lev],
                          * implicit_functions[a_lev + 1], true);
           implicit_functions[a_lev + 1] = std::move(new_impfunc);
       }

       // This is broken at the moment => using mfix::intersect_ls_walls () instead
       // Print() << "Modifying level set to see inflow" << std::endl;
       // mfix_set_ls_near_inflow(); //TODO move the level-set creation
       // Print() << "Done regridding level-set on lev = " << a_lev << std::endl;
   }
}


//!  This function checks if ebfactory is allocated with the proper dm and ba
bool mfix::mfix_update_ebfactory (int a_lev)
{
   // This assert is to verify that some kind of EB geometry has already been
   // defined
   AMREX_ASSERT(not EB2::IndexSpace::empty());

   const DistributionMapping & dm = DistributionMap(a_lev);
   const BoxArray &            ba = boxArray(a_lev);

   bool is_updated = false;

   if ( ebfactory[a_lev].get() == nullptr )
   {
      Print() << "Updating ebfactory" << std::endl;

      ebfactory[a_lev].reset(
          new EBFArrayBoxFactory(* eb_levels[a_lev], geom[a_lev], ba, dm,
                                 {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                  m_eb_full_grow_cells}, m_eb_support_level)
          );

      is_updated = true;
   }
   else
   {

      const DistributionMapping&  eb_dm = ebfactory[a_lev]->DistributionMap();
      const BoxArray&             eb_ba = ebfactory[a_lev]->boxArray();

      if ( (dm != eb_dm) || (ba != eb_ba) )
      {
          Print() << "Updating ebfactory" << std::endl;

          ebfactory[a_lev].reset(
              new EBFArrayBoxFactory(* eb_levels[a_lev], geom[a_lev], ba, dm,
                                     {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                      m_eb_full_grow_cells}, m_eb_support_level)
              );

          is_updated = true;
      }
   }

   return is_updated;
}
