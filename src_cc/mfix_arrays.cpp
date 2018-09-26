#include <mfix_level.H>

void
mfix_level::AllocateArrays (int lev)
{
    // ********************************************************************************
    // Cell- or node-based arrays
    // ********************************************************************************

    // Void fraction
    ep_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    ep_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    ep_g[lev]->setVal(1.0);
    ep_go[lev]->setVal(1.0);

    // Gas density
    ro_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    ro_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    ro_g[lev]->setVal(0.);
    ro_go[lev]->setVal(0.);

    // Gas bulk density
    rop_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    rop_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    rop_g[lev]->setVal(0.);
    rop_go[lev]->setVal(0.);

    if (nodal_pressure)
    {
       const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

       p0_g[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0));
        p_g[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0));
       p_go[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0));
       pp_g[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0));

    } else {

       p0_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
        p_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
       p_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
       pp_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));

    }

     p0_g[lev]->setVal(0.);
      p_g[lev]->setVal(0.);
     p_go[lev]->setVal(0.);
     pp_g[lev]->setVal(0.);

    // Presssure gradients
    gp[lev].reset(new MultiFab(grids[lev],dmap[lev],3,nghost));
    gp0[lev].reset(new MultiFab(grids[lev],dmap[lev],3,nghost));
    gp[lev]->setVal(0.);
    gp0[lev]->setVal(0.);

    // Molecular viscosity
    mu_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    mu_g[lev]->setVal(0.);

    // Coefficient of grad(div(u)) in viscous terms
    lambda_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    lambda_g[lev]->setVal(0.);

    // Current velocity
    vel_g[lev].reset(new MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    vel_g[lev]->setVal(0.);

    // Old velocity
    vel_go[lev].reset(new  MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    vel_go[lev]->setVal(0.);

    // Div(u)
    trD_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    trD_g[lev]->setVal(0.);

    // Vorticity
    vort[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    vort[lev]->setVal(0.);

    // This is the deposition onto the grid of the beta coefficient
    // for fluid vel in the expression beta*(fluid_vel _ particle_vel)
    // Note this only needs one component since all velocity components are co-located
    f_gds[lev].reset(new  MultiFab(grids[lev],dmap[lev],1,nghost));
    f_gds[lev]->setVal(0.);

    // This is the deposition onto the grid of the drag term experienced by the particle
    drag[lev].reset(new  MultiFab(grids[lev],dmap[lev],3,nghost));
    drag[lev]->setVal(0.);

    if (nodal_pressure)
    {
       const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

        phi[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0));
      diveu[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0));

    } else {

        phi[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
      diveu[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));

    }

      phi[lev]->setVal(0.);
    diveu[lev]->setVal(0.);

    // Arrays to store the solution and rhs for the diffusion solve
    phi_diff[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    rhs_diff[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));

    // Slopes in x-direction
    xslopes[lev].reset(new  MultiFab(grids[lev],dmap[lev],3,nghost));
    xslopes[lev] -> setVal(0.);

    // Slopes in y-direction
    yslopes[lev].reset(new  MultiFab(grids[lev],dmap[lev],3,nghost));
    yslopes[lev] -> setVal(0.);

    // Slopes in z-direction
    zslopes[lev].reset(new  MultiFab(grids[lev],dmap[lev],3,nghost));
    zslopes[lev] -> setVal(0.);

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

    // When the pressure is on nodes, bcoeff is at cell centers
    if (nodal_pressure)
    {
       bcoeff[lev][0].reset(new  MultiFab(grids[lev],dmap[lev],1,nghost));
       bcoeff[lev][1].reset(new  MultiFab(grids[lev],dmap[lev],1,nghost));
       bcoeff[lev][2].reset(new  MultiFab(grids[lev],dmap[lev],1,nghost));

    } else {

       // Create a BoxArray on x-faces.
       BoxArray x_edge_ba = grids[lev];
       x_edge_ba.surroundingNodes(0);
       bcoeff[lev][0].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));

       // Create a BoxArray on y-faces.
       BoxArray y_edge_ba = grids[lev];
       y_edge_ba.surroundingNodes(1);
       bcoeff[lev][1].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));

       // Create a BoxArray on y-faces.
       BoxArray z_edge_ba = grids[lev];
       z_edge_ba.surroundingNodes(2);
       bcoeff[lev][2].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    }

    bcoeff[lev][0]->setVal(0.);
    bcoeff[lev][1]->setVal(0.);
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
mfix_level::RegridArrays (int lev, BoxArray& new_grids, DistributionMapping& new_dmap)
{

    const EB2::IndexSpace & ebis = EB2::IndexSpace::top();
    const EB2::Level & ebis_level  = ebis.getLevel(geom[lev]);

    ebfactory[lev].reset(new EBFArrayBoxFactory( ebis_level,
                             geom[lev], new_grids, new_dmap,
                            {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                             m_eb_full_grow_cells}, m_eb_support_level)
                  );

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
    std::unique_ptr<MultiFab> ep_g_new(new MultiFab(new_grids,new_dmap,1,ng));
    ep_g_new->setVal(1.0);
    ep_g_new->copy(*ep_g[lev],0,0,1,0,ng);
    ep_g[lev] = std::move(ep_g_new);

    // Old void fraction
    ng = ep_go[lev]->nGrow();
    std::unique_ptr<MultiFab> ep_go_new(new MultiFab(new_grids,new_dmap,1,ep_go[lev]->nGrow()));
    ep_go_new->setVal(1.0);
    ep_go_new->copy(*ep_go[lev],0,0,1,0,ng);
    ep_go[lev] = std::move(ep_go_new);

    // Gas density
    ng = ro_g[lev]->nGrow();
    std::unique_ptr<MultiFab> ro_g_new(new MultiFab(new_grids,new_dmap,1,ro_g[lev]->nGrow()));
    ro_g_new->setVal(0.0);
    ro_g_new->copy(*ro_g[lev],0,0,1,0,ng);
    ro_g[lev] = std::move(ro_g_new);

    // Old gas density
    ng = ro_go[lev]->nGrow();
    std::unique_ptr<MultiFab> ro_go_new(new MultiFab(new_grids,new_dmap,1,ro_go[lev]->nGrow()));
    ro_go_new->setVal(0.0);
    ro_go_new->copy(*ro_go[lev],0,0,1,ng,ng);
    ro_go[lev] = std::move(ro_go_new);

    // Gas bulk density
    ng = rop_g[lev]->nGrow();
    std::unique_ptr<MultiFab> rop_g_new(new MultiFab(new_grids,new_dmap,1,rop_g[lev]->nGrow()));
    rop_g_new->setVal(0.0);
    rop_g_new->copy(*rop_g[lev],0,0,1,0,ng);
    rop_g[lev] = std::move(rop_g_new);

    // Old gas bulk density
    ng = rop_go[lev]->nGrow();
    std::unique_ptr<MultiFab> rop_go_new(new MultiFab(new_grids,new_dmap,1,rop_go[lev]->nGrow()));
    rop_go_new->setVal(0.0);
    rop_go_new->copy(*rop_go[lev],0,0,1,0,ng);
    rop_go[lev] = std::move(rop_go_new);

    if (nodal_pressure)
    {
       const BoxArray & nd_grids = amrex::convert(new_grids, IntVect{1,1,1});

       ng = p_g[lev]->nGrow();
       std::unique_ptr<MultiFab> p_g_new(new MultiFab(nd_grids,new_dmap,1,ng));
       p_g_new->setVal(0.0);
       p_g_new->copy(*p_g[lev],0,0,1,0,ng);
       p_g[lev] = std::move(p_g_new);

       ng = p_go[lev]->nGrow();
       std::unique_ptr<MultiFab> p_go_new(new MultiFab(nd_grids,new_dmap,1,ng));
       p_go_new->setVal(0.0);
       p_go_new->copy(*p_go[lev],0,0,1,0,ng);
       p_go[lev] = std::move(p_go_new);

       ng = p0_g[lev]->nGrow();
       std::unique_ptr<MultiFab> p0_g_new(new MultiFab(nd_grids,new_dmap,1,ng));
       p0_g_new->setVal(0.0);
       p0_g_new->copy(*p0_g[lev],0,0,1,0,ng);
       p0_g[lev] = std::move(p0_g_new);

       ng = pp_g[lev]->nGrow();
       std::unique_ptr<MultiFab> pp_g_new(new MultiFab(nd_grids,new_dmap,1,ng));
       pp_g_new->FillBoundary(p0_periodicity);
       pp_g_new->copy(*pp_g[lev],0,0,1,0,ng);
       pp_g[lev] = std::move(pp_g_new);

       std::unique_ptr<MultiFab> diveu_new(new MultiFab(nd_grids,new_dmap,1,diveu[lev]->nGrow()));
       diveu[lev] = std::move(diveu_new);
       diveu[lev]->setVal(0.);

       std::unique_ptr<MultiFab> phi_new(new MultiFab(nd_grids,new_dmap,1,phi[lev]->nGrow()));
       phi[lev] = std::move(phi_new);
       phi[lev]->setVal(0.);

       std::unique_ptr<MultiFab> bc0_new(new MultiFab(new_grids,new_dmap,1,bcoeff[lev][0]->nGrow(),
                                                     MFInfo(), *ebfactory[lev] ));
       bcoeff[lev][0] = std::move(bc0_new);
       bcoeff[lev][0]->setVal(0.);

       std::unique_ptr<MultiFab> bc1_new(new MultiFab(new_grids,new_dmap,1,bcoeff[lev][1]->nGrow(),
                                                     MFInfo(), *ebfactory[lev] ));
       bcoeff[lev][1] = std::move(bc1_new);
       bcoeff[lev][1]->setVal(0.);

       std::unique_ptr<MultiFab> bc2_new(new MultiFab(new_grids,new_dmap,1,bcoeff[lev][2]->nGrow(),
                                                     MFInfo(), *ebfactory[lev] ));
       bcoeff[lev][2] = std::move(bc2_new);
       bcoeff[lev][2]->setVal(0.);

    } else {

       ng = p_g[lev]->nGrow();
       std::unique_ptr<MultiFab> p_g_new(new MultiFab(new_grids,new_dmap,1,ng));
       p_g_new->setVal(0.0);
       p_g_new->copy(*p_g[lev],0,0,1,0,ng);
       p_g[lev] = std::move(p_g_new);

       ng = p_go[lev]->nGrow();
       std::unique_ptr<MultiFab> p_go_new(new MultiFab(new_grids,new_dmap,1,ng));
       p_go_new->setVal(0.0);
       p_go_new->copy(*p_go[lev],0,0,1,0,ng);
       p_go[lev] = std::move(p_go_new);

       ng = p0_g[lev]->nGrow();
       std::unique_ptr<MultiFab> p0_g_new(new MultiFab(new_grids,new_dmap,1,ng));
       p0_g_new->setVal(0.0);
       p0_g_new->copy(*p0_g[lev],0,0,1,0,ng);
       p0_g[lev] = std::move(p0_g_new);

       ng = pp_g[lev]->nGrow();
       std::unique_ptr<MultiFab> pp_g_new(new MultiFab(new_grids,new_dmap,1,ng));
       pp_g[lev]->setVal(0.0);
       pp_g_new->copy(*pp_g[lev],0,0,1,0,ng);
       pp_g[lev] = std::move(pp_g_new);

       std::unique_ptr<MultiFab> phi_new(new MultiFab(new_grids,new_dmap,1,phi[lev]->nGrow(),
                                                      MFInfo(), *ebfactory[lev]));
       phi[lev] = std::move(phi_new);
       phi[lev]->setVal(0.);

       std::unique_ptr<MultiFab> diveu_new(new MultiFab(new_grids,new_dmap,1,diveu[lev]->nGrow(),
                                                      MFInfo(), *ebfactory[lev]));
       diveu[lev] = std::move(diveu_new);
       diveu[lev]->setVal(0.);

       // Cell-centered pressure uses face-based coefficients
       BoxArray x_ba = new_grids;
       x_ba = x_ba.surroundingNodes(0);
       std::unique_ptr<MultiFab> bc0_new(new MultiFab(x_ba,new_dmap,1,nghost));
       bcoeff[lev][0] = std::move(bc0_new);
       bcoeff[lev][0] -> setVal(0.0);

       BoxArray y_ba = new_grids;
       y_ba = y_ba.surroundingNodes(1);
       std::unique_ptr<MultiFab> bc1_new(new MultiFab(y_ba,new_dmap,1,nghost));
       bcoeff[lev][1] = std::move(bc1_new);
       bcoeff[lev][1] -> setVal(0.0);

       BoxArray z_ba = new_grids;
       z_ba = z_ba.surroundingNodes(2);
       std::unique_ptr<MultiFab> bc2_new(new MultiFab(z_ba,new_dmap,1,nghost));
       bcoeff[lev][2] = std::move(bc2_new);
       bcoeff[lev][2] -> setVal(0.0);
    }

    // Molecular viscosity
    ng = mu_g[lev]->nGrow();
    std::unique_ptr<MultiFab> mu_g_new(new MultiFab(new_grids,new_dmap,1,ng));
    mu_g_new->setVal(0.0);
    mu_g_new->copy(*mu_g[lev],0,0,1,0,ng);
    mu_g[lev] = std::move(mu_g_new);

    // Lambda
    ng = lambda_g[lev]->nGrow();
    std::unique_ptr<MultiFab> lambda_g_new(new MultiFab(new_grids,new_dmap,1,ng));
    lambda_g_new->setVal(0.0);
    lambda_g_new->copy(*lambda_g[lev],0,0,1,0,ng);
    lambda_g[lev] = std::move(lambda_g_new);

    // Gas velocity
    ng = vel_g[lev]->nGrow();
    std::unique_ptr<MultiFab> vel_g_new(new MultiFab(new_grids,new_dmap,vel_g[lev]->nComp(),ng));
    vel_g_new->setVal(0.0);
    vel_g_new->copy(*vel_g[lev],0,0,vel_g[lev]->nComp(),0,ng);
    vel_g[lev] = std::move(vel_g_new);

    // Old gas velocity
    ng = vel_go[lev]->nGrow();
    std::unique_ptr<MultiFab> vel_go_new(new MultiFab(new_grids,new_dmap,vel_go[lev]->nComp(),ng));
    vel_go_new->setVal(0.0);
    vel_go_new->copy(*vel_go[lev],0,0,vel_go[lev]->nComp(),0,ng);
    vel_go[lev] = std::move(vel_go_new);

    // Pressure gradients
    ng = gp[lev]->nGrow();
    std::unique_ptr<MultiFab> gp_new(new MultiFab(new_grids,new_dmap,3,ng));
    gp_new->setVal(0.0);
    gp_new->copy(*gp[lev],0,0,1,0,ng);
    gp[lev] = std::move(gp_new);

    // Pressure gradients
    ng = gp0[lev]->nGrow();
    std::unique_ptr<MultiFab> gp0_new(new MultiFab(new_grids,new_dmap,3,ng));
    gp0_new->setVal(0.0);
    gp0_new->copy(*gp0[lev],0,0,1,0,ng);
    gp0[lev] = std::move(gp0_new);

    // Trace(D)
    ng = trD_g[lev]->nGrow();
    std::unique_ptr<MultiFab> trD_g_new(new MultiFab(new_grids,new_dmap,1,ng));
    trD_g[lev] = std::move(trD_g_new);
    trD_g[lev]->setVal(0.);

    // Vorticity
    ng = vort[lev]->nGrow();
    std::unique_ptr<MultiFab> vort_new(new MultiFab(new_grids,new_dmap,1,ng));
    vort[lev] = std::move(vort_new);
    vort[lev]->setVal(0.);

    // Coefficient in drag
    ng = f_gds[lev]->nGrow();
    std::unique_ptr<MultiFab> f_gds_new(new MultiFab(new_grids,new_dmap,f_gds[lev]->nComp(),ng));
    f_gds[lev] = std::move(f_gds_new);
    f_gds[lev]->setVal(0.);

    // Particle/fluid drag
    ng = drag[lev]->nGrow();
    std::unique_ptr<MultiFab> drag_new(new MultiFab(new_grids,new_dmap,drag[lev]->nComp(),ng));
    drag[lev] = std::move(drag_new);
    drag[lev]->setVal(0.);

   /****************************************************************************
    * Face-based Arrays                                                        *
    ****************************************************************************/

    BoxArray x_ba = new_grids;
    x_ba = x_ba.surroundingNodes(0);

    // MAC velocity
    std::unique_ptr<MultiFab> u_mac_new(new MultiFab(x_ba,new_dmap,1,nghost));
    m_u_mac[lev] = std::move(u_mac_new);
    m_u_mac[lev] -> setVal(0.0);

    // Diffusion coefficient on x-faces
    std::unique_ptr<MultiFab> bc0_new(new MultiFab(x_ba,new_dmap,1,nghost));
    bcoeff_diff[lev][0] = std::move(bc0_new);
    bcoeff_diff[lev][0] -> setVal(0.0);

   //****************************************************************************

    BoxArray y_ba = new_grids;
    y_ba = y_ba.surroundingNodes(1);

    // MAC velocity
    std::unique_ptr<MultiFab> v_mac_new(new MultiFab(y_ba,new_dmap,1,nghost));
    m_v_mac[lev] = std::move(v_mac_new);
    m_v_mac[lev] -> setVal(0.0);

    // Diffusion coefficient on y-faces
    std::unique_ptr<MultiFab> bc1_new(new MultiFab(y_ba,new_dmap,1,nghost));
    bcoeff_diff[lev][1] = std::move(bc1_new);
    bcoeff_diff[lev][1] -> setVal(0.0);

   //****************************************************************************

    BoxArray z_ba = new_grids;
    z_ba = z_ba.surroundingNodes(2);

    // MAC velocity
    std::unique_ptr<MultiFab> w_mac_new(new MultiFab(z_ba,new_dmap,1,nghost));
    m_w_mac[lev] = std::move(w_mac_new);
    m_w_mac[lev] -> setVal(0.0);

    // Diffusion coefficient on z-faces
    std::unique_ptr<MultiFab> bc2_new(new MultiFab(z_ba,new_dmap,1,nghost));
    bcoeff[lev][2] = std::move(bc2_new);
    bcoeff[lev][2] -> setVal(0.0);

   /****************************************************************************
    * Nodal Arrays                                                             *
    ****************************************************************************/

    // Create a nodal BoxArray
    const BoxArray & new_nodal_grids = amrex::convert(new_grids, IntVect{1,1,1});

    // Level-set
    // For the level set we set src_nghost to ng as well, even if it's potentially
    // not totally safe. This way we have the domain ghost cells filled properly.
    ng = ls[lev]->nGrow();
    std::unique_ptr<MultiFab> ls_new(new MultiFab(new_nodal_grids, new_dmap, 1, ng ));
    ls_new->copy(*ls[lev],0,0,1,ng,ng);
    ls[lev] = std::move(ls_new);


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

    if (!nodal_pressure){
       fill_mf_bc(lev,*p_g[lev]);
       fill_mf_bc(lev,*p_go[lev]);
    }
}
