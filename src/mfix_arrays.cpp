#include <mfix_level.H>

void
mfix_level::AllocateArrays (int lev)
{
    // ********************************************************************************
    // Cell-based arrays
    // ********************************************************************************

    // Void fraction
    ep_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    ep_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    ep_g[lev]->setVal(1.0);
    ep_go[lev]->setVal(1.0);

    // Gas pressure fraction
    p_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    p_go[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    p_g[lev]->setVal(0.);
    p_go[lev]->setVal(0.);

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

    // Base pressure that captures delp and/or p_in and p_out
    p0_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    p0_g[lev]->setVal(0.);

    // Pressure correction equation
    pp_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    pp_g[lev]->setVal(0.);

    // Molecular viscosity
    mu_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    mu_g[lev]->setVal(0.);

    //
    lambda_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    lambda_g[lev]->setVal(0.);
    //
    trD_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    trD_g[lev]->setVal(0.);

    // Vorticity
    vort[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    vort[lev]->setVal(0.);

    // Pressure increment
    phi[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    phi[lev]->setVal(0.);

    // div(e u)
    diveu[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost));
    diveu[lev]->setVal(0.);

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

    // Create a BoxArray on x-faces.
    BoxArray x_edge_ba = grids[lev];
    x_edge_ba.surroundingNodes(0);

    // X-axis gas velocity
    u_g[lev].reset(new MultiFab(x_edge_ba,dmap[lev],1,nghost));
    u_go[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    u_gt[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    fp_x[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    bcoeff[lev][0].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    u_g[lev]->setVal(0.);
    u_go[lev]->setVal(0.);
    u_gt[lev]->setVal(0.);
    fp_x[lev]->setVal(0.);
    bcoeff[lev][0]->setVal(0.);

    d_e[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    d_e[lev]->setVal(0.);

    tau_u_g[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    tau_u_g[lev]->setVal(0.);

    fluxX[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    fluxX[lev]->setVal(0.);

    ropX[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    ropX[lev]->setVal(0.);

    f_gds_u[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,1));
    f_gds_u[lev]->setVal(0.);

    drag_u[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,1));
    drag_u[lev]->setVal(0.);

    // u-velocity slopes. Note that the number of components is not 1, but 3!
    slopes_u[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],3,nghost));
    slopes_u[lev] -> setVal(0.);

    // u acceleration terms
    uacc[lev].reset(new  MultiFab(x_edge_ba,dmap[lev],1,nghost));
    uacc[lev] -> setVal(0.);

    // ********************************************************************************
    // Y-face-based arrays
    // ********************************************************************************

    // Create a BoxArray on y-faces.
    BoxArray y_edge_ba = grids[lev];
    y_edge_ba.surroundingNodes(1);

    // Y-axis gas velocity
    v_g[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    v_go[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    v_gt[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    fp_y[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    bcoeff[lev][1].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    v_g[lev]->setVal(0.);
    v_go[lev]->setVal(0.);
    v_gt[lev]->setVal(0.);
    fp_y[lev]->setVal(0.);
    bcoeff[lev][1]->setVal(0.);

    d_n[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost));
    d_n[lev]->setVal(0.);

    tau_v_g[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost));
    tau_v_g[lev]->setVal(0.);

    fluxY[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost));
    fluxY[lev]->setVal(0.);

    ropY[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost));
    ropY[lev]->setVal(0.);

    f_gds_v[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,1));
    f_gds_v[lev]->setVal(0.);

    drag_v[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,1));
    drag_v[lev]->setVal(0.);

    // v-velocity slopes. Note that the number of components is not 1, but 3!
    slopes_v[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],3,nghost));
    slopes_v[lev] -> setVal(0.);

    // v acceleration terms
    vacc[lev].reset(new  MultiFab(y_edge_ba,dmap[lev],1,nghost));
    vacc[lev] -> setVal(0.);



    // ********************************************************************************
    // Z-face-based arrays
    // ********************************************************************************

    // Create a BoxArray on y-faces.
    BoxArray z_edge_ba = grids[lev];
    z_edge_ba.surroundingNodes(2);

    // Z-axis gas velocity
    w_g[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    w_go[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    w_gt[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    fp_z[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    bcoeff[lev][2].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    w_g[lev]->setVal(0.);
    w_go[lev]->setVal(0.);
    w_gt[lev]->setVal(0.);
    fp_z[lev]->setVal(0.);
    bcoeff[lev][2]->setVal(0.);

    d_t[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    d_t[lev]->setVal(0.);

    tau_w_g[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    tau_w_g[lev]->setVal(0.);

    fluxZ[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    fluxZ[lev]->setVal(0.);

    ropZ[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,nghost));
    ropZ[lev]->setVal(0.);

    f_gds_w[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,1));
    f_gds_w[lev]->setVal(0.);

    drag_w[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,1));
    drag_w[lev]->setVal(0.);

    // w-velocity slopes. Note that the number of components is not 1, but 3!
    slopes_w[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],3,nghost));
    slopes_w[lev] -> setVal(0.);

    // w acceleration terms
    wacc[lev].reset(new  MultiFab(z_edge_ba,dmap[lev],1,nghost));
    wacc[lev] -> setVal(0.);

}

void
mfix_level::RegridArrays (int lev, BoxArray& new_grids, DistributionMapping& new_dmap)
{
   /****************************************************************************
    * Cell-based arrays                                                        *
    ****************************************************************************/

    // Void fraction
    int ng = ep_g[lev]->nGrow();
    std::unique_ptr<MultiFab> ep_g_new(new MultiFab(new_grids,new_dmap,1,ng));
    ep_g_new->copy(*ep_g[lev],0,0,1,ng,ng);
    ep_g_new->FillBoundary(geom[lev].periodicity());
    ep_g[lev] = std::move(ep_g_new);

    // Old void fraction
    ng = ep_go[lev]->nGrow();
    std::unique_ptr<MultiFab> ep_go_new(new MultiFab(new_grids,new_dmap,1,ep_go[lev]->nGrow()));
    ep_go_new->copy(*ep_go[lev],0,0,1,ng,ng);
    ep_go_new->FillBoundary(geom[lev].periodicity());
    ep_go[lev] = std::move(ep_go_new);

    // Gas pressure fraction
    ng = p_g[lev]->nGrow();
    std::unique_ptr<MultiFab> p_g_new(new MultiFab(new_grids,new_dmap,1,p_g[lev]->nGrow()));
    p_g_new->copy(*p_g[lev],0,0,1,ng,ng);
    p_g_new->FillBoundary(geom[lev].periodicity());
    p_g[lev] = std::move(p_g_new);

    // Old gas pressure fraction
    ng = p_go[lev]->nGrow();
    std::unique_ptr<MultiFab> p_go_new(new MultiFab(new_grids,new_dmap,1,p_go[lev]->nGrow()));
    p_go_new->copy(*p_go[lev],0,0,1,ng,ng);
    p_go_new->FillBoundary(geom[lev].periodicity());
    p_go[lev] = std::move(p_go_new);

    // Gas density
    ng = ro_g[lev]->nGrow();
    std::unique_ptr<MultiFab> ro_g_new(new MultiFab(new_grids,new_dmap,1,ro_g[lev]->nGrow()));
    ro_g_new->copy(*ro_g[lev],0,0,1,ng,ng);
    ro_g_new->FillBoundary(geom[lev].periodicity());
    ro_g[lev] = std::move(ro_g_new);

    // Old gas density
    ng = ro_go[lev]->nGrow();
    std::unique_ptr<MultiFab> ro_go_new(new MultiFab(new_grids,new_dmap,1,ro_go[lev]->nGrow()));
    ro_go_new->copy(*ro_go[lev],0,0,1,ng,ng);
    ro_go_new->FillBoundary(geom[lev].periodicity());
    ro_go[lev] = std::move(ro_go_new);

    // Gas bulk density
    ng = rop_g[lev]->nGrow();
    std::unique_ptr<MultiFab> rop_g_new(new MultiFab(new_grids,new_dmap,1,rop_g[lev]->nGrow()));
    rop_g_new->copy(*rop_g[lev],0,0,1,ng,ng);
    rop_g_new->FillBoundary(geom[lev].periodicity());
    rop_g[lev] = std::move(rop_g_new);

    // Old gas bulk density
    ng = rop_go[lev]->nGrow();
    std::unique_ptr<MultiFab> rop_go_new(new MultiFab(new_grids,new_dmap,1,rop_go[lev]->nGrow()));
    rop_go_new->copy(*rop_go[lev],0,0,1,ng,ng);
    rop_go_new->FillBoundary(geom[lev].periodicity());
    rop_go[lev] = std::move(rop_go_new);

    // Base pressure
    ng = p0_g[lev]->nGrow();
    std::unique_ptr<MultiFab> p0_g_new(new MultiFab(new_grids,new_dmap,1,p0_g[lev]->nGrow()));
    p0_g_new->copy(*p0_g[lev],0,0,1,ng,ng);
    p0_g_new->FillBoundary(p0_periodicity);
    p0_g[lev] = std::move(p0_g_new);

    // Pressure correction
    ng = pp_g[lev]->nGrow();
    std::unique_ptr<MultiFab> pp_g_new(new MultiFab(new_grids,new_dmap,1,pp_g[lev]->nGrow()));
    pp_g_new->copy(*pp_g[lev],0,0,1,ng,ng);
    pp_g_new->FillBoundary(geom[lev].periodicity());
    pp_g[lev] = std::move(pp_g_new);

    // Molecular viscosity
    ng = mu_g[lev]->nGrow();
    std::unique_ptr<MultiFab> mu_g_new(new MultiFab(new_grids,new_dmap,1,mu_g[lev]->nGrow()));
    mu_g_new->copy(*mu_g[lev],0,0,1,ng,ng);
    mu_g_new->FillBoundary(geom[lev].periodicity());
    mu_g[lev] = std::move(mu_g_new);

    // Lambda
    ng = lambda_g[lev]->nGrow();
    std::unique_ptr<MultiFab> lambda_g_new(new MultiFab(new_grids,new_dmap,1,lambda_g[lev]->nGrow()));
    lambda_g_new->copy(*lambda_g[lev],0,0,1,ng,ng);
    lambda_g_new->FillBoundary(geom[lev].periodicity());
    lambda_g[lev] = std::move(lambda_g_new);

    // Trace(D)
    ng = trD_g[lev]->nGrow();
    std::unique_ptr<MultiFab> trD_g_new(new MultiFab(new_grids,new_dmap,1,trD_g[lev]->nGrow()));
    trD_g_new->copy(*trD_g[lev],0,0,1,ng,ng);
    trD_g_new->FillBoundary(geom[lev].periodicity());
    trD_g[lev] = std::move(trD_g_new);

    // Vorticity
    ng = vort[lev]->nGrow();
    std::unique_ptr<MultiFab> vort_new(new MultiFab(new_grids,new_dmap,1,vort[lev]->nGrow()));
    vort_new->copy(*vort[lev],0,0,1,ng,ng);
    vort_new->FillBoundary(geom[lev].periodicity());
    vort[lev] = std::move(vort_new);

    // Pressure increment (phi)
    ng = phi[lev]->nGrow();
    std::unique_ptr<MultiFab> phi_new(new MultiFab(new_grids,new_dmap,1,diveu[lev]->nGrow()));
    phi_new->copy(*phi[lev],0,0,1,ng,ng);
    phi_new->FillBoundary(geom[lev].periodicity());
    phi[lev] = std::move(phi_new);

    // Diveu
    ng = diveu[lev]->nGrow();
    std::unique_ptr<MultiFab> diveu_new(new MultiFab(new_grids,new_dmap,1,diveu[lev]->nGrow()));
    diveu_new->copy(*diveu[lev],0,0,1,ng,ng);
    diveu_new->FillBoundary(geom[lev].periodicity());
    diveu[lev] = std::move(diveu_new);

   /****************************************************************************
    * X-face-based arrays                                                      *
    ****************************************************************************/

    // Create a BoxArray on x-faces.
    BoxArray x_edge_ba = new_grids;
    x_edge_ba.surroundingNodes(0);

    // X-axis gas velocity
    ng = u_g[lev]->nGrow();
    std::unique_ptr<MultiFab> u_g_new(new MultiFab(x_edge_ba,new_dmap,1,ng));
    u_g_new->copy(*u_g[lev],0,0,1,ng,ng);
    u_g_new->FillBoundary(geom[lev].periodicity());
    u_g[lev] = std::move(u_g_new);

    ng = u_go[lev]->nGrow();
    std::unique_ptr<MultiFab> u_go_new(new MultiFab(x_edge_ba,new_dmap,1,ng));
    u_go_new->copy(*u_go[lev],0,0,1,ng,ng);
    u_go_new->FillBoundary(geom[lev].periodicity());
    u_go[lev] = std::move(u_go_new);

    ng = u_gt[lev]->nGrow();
    std::unique_ptr<MultiFab> u_gt_new(new MultiFab(x_edge_ba,new_dmap,1,ng));
    u_gt_new->copy(*u_gt[lev],0,0,1,ng,ng);
    u_gt_new->FillBoundary(geom[lev].periodicity());
    u_gt[lev] = std::move(u_gt_new);

    ng = ropX[lev]->nGrow();
    std::unique_ptr<MultiFab> ropX_new(new MultiFab(x_edge_ba,new_dmap,1,ng));
    ropX_new->copy(*ropX[lev],0,0,1,ng,ng);
    ropX_new->FillBoundary(geom[lev].periodicity());
    ropX[lev] = std::move(ropX_new);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> d_e_new(new MultiFab(x_edge_ba,new_dmap,1,d_e[lev]->nGrow()));
    d_e[lev] = std::move(d_e_new);
    d_e[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> tau_u_g_new(new MultiFab(x_edge_ba,new_dmap,1,tau_u_g[lev]->nGrow()));
    tau_u_g[lev] = std::move(tau_u_g_new);
    tau_u_g[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> fluxX_new(new MultiFab(x_edge_ba,new_dmap,1,fluxX[lev]->nGrow()));
    fluxX[lev] = std::move(fluxX_new);
    fluxX[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> f_gds_u_new(new MultiFab(x_edge_ba,new_dmap,1,f_gds_u[lev]->nGrow()));
    f_gds_u[lev] = std::move(f_gds_u_new);
    f_gds_u[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> drag_u_new(new MultiFab(x_edge_ba,new_dmap,1,drag_u[lev]->nGrow()));
    drag_u[lev] = std::move(drag_u_new);
    drag_u[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> bcoeffx(new MultiFab(x_edge_ba,new_dmap,1,bcoeff[lev][0]->nGrow()));
    bcoeff[lev][0] = std::move(bcoeffx);
    bcoeff[lev][0]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> slopes_u_new(new MultiFab(x_edge_ba,new_dmap,3,slopes_u[lev]->nGrow()));
    slopes_u[lev] = std::move(slopes_u_new);
    slopes_u[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> uacc_new(new MultiFab(x_edge_ba,new_dmap,1,uacc[lev]->nGrow()));
    uacc[lev] = std::move(uacc_new);
    uacc[lev]->setVal(0.);


   /****************************************************************************
    * Y-face-based arrays                                                      *
    ****************************************************************************/

    // Create a BoxArray on y-faces.
    BoxArray y_edge_ba = new_grids;
    y_edge_ba.surroundingNodes(1);

    // Y-axis gas velocity
    ng = v_g[lev]->nGrow();
    std::unique_ptr<MultiFab> v_g_new(new MultiFab(y_edge_ba,new_dmap,1,ng));
    v_g_new->copy(*v_g[lev],0,0,1,ng,ng);
    v_g_new->FillBoundary(geom[lev].periodicity());
    v_g[lev] = std::move(v_g_new);

    ng = v_go[lev]->nGrow();
    std::unique_ptr<MultiFab> v_go_new(new MultiFab(y_edge_ba,new_dmap,1,ng));
    v_go_new->copy(*v_go[lev],0,0,1,ng,ng);
    v_go_new->FillBoundary(geom[lev].periodicity());
    v_go[lev] = std::move(v_go_new);

    ng = v_gt[lev]->nGrow();
    std::unique_ptr<MultiFab> v_gt_new(new MultiFab(y_edge_ba,new_dmap,1,ng));
    v_gt_new->copy(*v_gt[lev],0,0,1,ng,ng);
    v_gt_new->FillBoundary(geom[lev].periodicity());
    v_gt[lev] = std::move(v_gt_new);

    ng = ropY[lev]->nGrow();
    std::unique_ptr<MultiFab> ropY_new(new MultiFab(y_edge_ba,new_dmap,1,ng));
    ropY_new->copy(*ropY[lev],0,0,1,ng,ng);
    ropY_new->FillBoundary(geom[lev].periodicity());
    ropY[lev] = std::move(ropY_new);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> d_n_new(new MultiFab(y_edge_ba,new_dmap,1,d_n[lev]->nGrow()));
    d_n[lev] = std::move(d_n_new);
    d_n[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> tau_v_g_new(new MultiFab(y_edge_ba,new_dmap,1,tau_v_g[lev]->nGrow()));
    tau_v_g[lev] = std::move(tau_v_g_new);
    tau_v_g[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> fluxY_new(new MultiFab(y_edge_ba,new_dmap,1,fluxY[lev]->nGrow()));
    fluxY[lev] = std::move(fluxY_new);
    fluxY[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> f_gds_v_new(new MultiFab(y_edge_ba,new_dmap,1,f_gds_v[lev]->nGrow()));
    f_gds_v[lev] = std::move(f_gds_v_new);
    f_gds_v[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> drag_v_new(new MultiFab(y_edge_ba,new_dmap,1,drag_v[lev]->nGrow()));
    drag_v[lev] = std::move(drag_v_new);
    drag_v[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> bcoeffy(new MultiFab(y_edge_ba,new_dmap,1,bcoeff[lev][1]->nGrow()));
    bcoeff[lev][1] = std::move(bcoeffy);
    bcoeff[lev][1]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    ng = slopes_v[lev]->nGrow();
    std::unique_ptr<MultiFab> slopes_v_new(new MultiFab(y_edge_ba,new_dmap,3,ng));
    slopes_v[lev] = std::move(slopes_v_new);
    slopes_v[lev]->setVal(0.);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> vacc_new(new MultiFab(y_edge_ba,new_dmap,1,vacc[lev]->nGrow()));
    vacc[lev] = std::move(vacc_new);
    vacc[lev]->setVal(0.);

   /****************************************************************************
    * Z-face-based arrays                                                      *
    ****************************************************************************/

    // Create a BoxArray on z-faces.
    BoxArray z_edge_ba = new_grids;
    z_edge_ba.surroundingNodes(2);

    // Z-axis gas velocity
    ng = w_g[lev]->nGrow();
    std::unique_ptr<MultiFab> w_g_new(new MultiFab(z_edge_ba,new_dmap,1,ng));
    w_g_new->copy(*w_g[lev],0,0,1,ng,ng);
    w_g_new->FillBoundary(geom[lev].periodicity());
    w_g[lev] = std::move(w_g_new);

    ng = w_go[lev]->nGrow();
    std::unique_ptr<MultiFab> w_go_new(new MultiFab(z_edge_ba,new_dmap,1,ng));
    w_go_new->copy(*w_go[lev],0,0,1,ng,ng);
    w_go_new->FillBoundary(geom[lev].periodicity());
    w_go[lev] = std::move(w_go_new);

    ng = w_gt[lev]->nGrow();
    std::unique_ptr<MultiFab> w_gt_new(new MultiFab(z_edge_ba,new_dmap,1,ng));
    w_gt_new->copy(*w_gt[lev],0,0,1,ng,ng);
    w_gt_new->FillBoundary(geom[lev].periodicity());
    w_gt[lev] = std::move(w_gt_new);

    ng = ropZ[lev]->nGrow();
    std::unique_ptr<MultiFab> ropZ_new(new MultiFab(z_edge_ba,new_dmap,1,ng));
    ropZ_new->copy(*ropZ[lev],0,0,1,ng,ng);
    ropZ_new->FillBoundary(geom[lev].periodicity());
    ropZ[lev] = std::move(ropZ_new);

    // These are temporaries so we don't need to copy in old data
    ng = d_t[lev]->nGrow();
    std::unique_ptr<MultiFab> d_t_new(new MultiFab(z_edge_ba,new_dmap,1,ng));
    d_t[lev] = std::move(d_t_new);
    d_t[lev]->setVal(0.0);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> tau_w_g_new(new MultiFab(z_edge_ba,new_dmap,1,tau_w_g[lev]->nGrow()));
    tau_w_g[lev] = std::move(tau_w_g_new);
    tau_w_g[lev]->setVal(0.0);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> fluxZ_new(new MultiFab(z_edge_ba,new_dmap,1,fluxZ[lev]->nGrow()));
    fluxZ[lev] = std::move(fluxZ_new);
    fluxZ[lev]->setVal(0.0);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> f_gds_w_new(new MultiFab(z_edge_ba,new_dmap,1,f_gds_w[lev]->nGrow()));
    f_gds_w[lev] = std::move(f_gds_w_new);
    f_gds_w[lev]->setVal(0.0);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> drag_w_new(new MultiFab(z_edge_ba,new_dmap,1,drag_w[lev]->nGrow()));
    drag_w[lev] = std::move(drag_w_new);
    drag_w[lev]->setVal(0.0);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> bcoeffz(new MultiFab(z_edge_ba,new_dmap,1,bcoeff[lev][2]->nGrow()));
    bcoeff[lev][2] = std::move(bcoeffz);
    bcoeff[lev][2]->setVal(0.0);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> slopes_w_new(new MultiFab(z_edge_ba,new_dmap,3,slopes_w[lev]->nGrow()));
    slopes_w[lev] = std::move(slopes_w_new);
    slopes_w[lev]->setVal(0.0);

    // These are temporaries so we don't need to copy in old data
    std::unique_ptr<MultiFab> wacc_new(new MultiFab(z_edge_ba,new_dmap,1,wacc[lev]->nGrow()));
    wacc[lev] = std::move(wacc_new);
    wacc[lev]->setVal(0.);

   /****************************************************************************
    * Nodal Arrays                                                             *
    ****************************************************************************/

    // Create a nodal BoxArray
    const BoxArray & new_nodal_grids = amrex::convert(new_grids, IntVect{1,1,1});

    // Level-set
    ng = ls[lev]->nGrow();
    std::unique_ptr<MultiFab> ls_new(new MultiFab(new_nodal_grids, new_dmap, 1, ls[lev]->nGrow()));
    ls_new->copy(*ls[lev],0,0,1,ng,ng);
    ls_new->FillBoundary(geom[lev].periodicity());
    ls[lev] = std::move(ls_new);


   /****************************************************************************
    * Make sure we fill the ghost cells as appropriate -- this is copied from  *
    * init_fluid                                                               *
    ****************************************************************************/

    fill_mf_bc(lev,*ep_g[lev]);
    fill_mf_bc(lev,*ep_go[lev]);
    fill_mf_bc(lev,*p_g[lev]);
    fill_mf_bc(lev,*p_go[lev]);
    fill_mf_bc(lev,*ro_g[lev]);
    fill_mf_bc(lev,*ro_go[lev]);
    fill_mf_bc(lev,*rop_g[lev]);
    fill_mf_bc(lev,*rop_go[lev]);

    fill_mf_bc(lev,*mu_g[lev]);
    fill_mf_bc(lev,*lambda_g[lev]);
}
