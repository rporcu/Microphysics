#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

void
mfix_level::Regrid (int lev, BoxArray& new_grids, DistributionMapping& new_dmap)
{
    // ********************************************************************************
    // Cell-based arrays
    // ********************************************************************************

    // Void fraction
    std::unique_ptr<MultiFab> ep_g_new(new MultiFab(new_grids,new_dmap,1,ep_g[lev]->nGrow()));
    ep_g_new->copy(*ep_g[lev]);
    ep_g[lev] = std::move(ep_g_new);

    // Gas pressure fraction
    std::unique_ptr<MultiFab> p_g_new(new MultiFab(new_grids,new_dmap,1,p_g[lev]->nGrow()));
    p_g_new->copy(*p_g[lev]);
    p_g[lev] = std::move(p_g_new);

    // Old gas pressure fraction
    std::unique_ptr<MultiFab> p_go_new(new MultiFab(new_grids,new_dmap,1,p_go[lev]->nGrow()));
    p_go_new->copy(*p_go[lev]);
    p_go[lev] = std::move(p_go_new);

    // Gas density
    std::unique_ptr<MultiFab> ro_g_new(new MultiFab(new_grids,new_dmap,1,ro_g[lev]->nGrow()));
    ro_g_new->copy(*ro_g[lev]);
    ro_g[lev] = std::move(ro_g_new);

    // Old gas density
    std::unique_ptr<MultiFab> ro_go_new(new MultiFab(new_grids,new_dmap,1,ro_go[lev]->nGrow()));
    ro_go_new->copy(*ro_go[lev]);
    ro_go[lev] = std::move(ro_go_new);

    // Gas bulk density
    std::unique_ptr<MultiFab> rop_g_new(new MultiFab(new_grids,new_dmap,1,rop_g[lev]->nGrow()));
    rop_g_new->copy(*rop_g[lev]);
    rop_g[lev] = std::move(rop_g_new);

    // Old gas bulk density
    std::unique_ptr<MultiFab> rop_go_new(new MultiFab(new_grids,new_dmap,1,rop_go[lev]->nGrow()));
    rop_go_new->copy(*rop_go[lev]);
    rop_go[lev] = std::move(rop_go_new);

    // Pressure correction equation
    std::unique_ptr<MultiFab> pp_g_new(new MultiFab(new_grids,new_dmap,1,pp_g[lev]->nGrow()));
    pp_g_new->copy(*pp_g[lev]);
    pp_g[lev] = std::move(pp_g_new);

    // Molecular viscosity 
    std::unique_ptr<MultiFab> mu_g_new(new MultiFab(new_grids,new_dmap,1,mu_g[lev]->nGrow()));
    mu_g_new->copy(*mu_g[lev]);
    mu_g[lev] = std::move(mu_g_new);

    // Lambda
    std::unique_ptr<MultiFab> lambda_g_new(new MultiFab(new_grids,new_dmap,1,lambda_g[lev]->nGrow()));
    lambda_g_new->copy(*lambda_g[lev]);
    lambda_g[lev] = std::move(lambda_g_new);

    // Trace(D)
    std::unique_ptr<MultiFab> trD_g_new(new MultiFab(new_grids,new_dmap,1,trD_g[lev]->nGrow()));
    trD_g_new->copy(*trD_g[lev]);
    trD_g[lev] = std::move(trD_g_new);

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

    // Create a BoxArray on x-faces.
    BoxArray x_edge_ba = new_grids;
    x_edge_ba.surroundingNodes(0);

    // X-axis gas velocity
    std::unique_ptr<MultiFab> u_g_new(new MultiFab(new_grids,new_dmap,1,u_g[lev]->nGrow()));
    u_g_new->copy(*u_g[lev]);
    u_g[lev] = std::move(u_g_new);

    std::unique_ptr<MultiFab> u_go_new(new MultiFab(new_grids,new_dmap,1,u_go[lev]->nGrow()));
    u_go_new->copy(*u_go[lev]);
    u_go[lev] = std::move(u_go_new);

    std::unique_ptr<MultiFab> u_gt_new(new MultiFab(new_grids,new_dmap,1,u_gt[lev]->nGrow()));
    u_gt_new->copy(*u_gt[lev]);
    u_gt[lev] = std::move(u_gt_new);

    std::unique_ptr<MultiFab> d_e_new(new MultiFab(new_grids,new_dmap,1,d_e[lev]->nGrow()));
    d_e_new->copy(*d_e[lev]);
    d_e[lev] = std::move(d_e_new);

    std::unique_ptr<MultiFab> tau_u_g_new(new MultiFab(new_grids,new_dmap,1,tau_u_g[lev]->nGrow()));
    tau_u_g_new->copy(*tau_u_g[lev]);
    tau_u_g[lev] = std::move(tau_u_g_new);

    std::unique_ptr<MultiFab> fluxX_new(new MultiFab(new_grids,new_dmap,1,fluxX[lev]->nGrow()));
    fluxX_new->copy(*fluxX[lev]);
    fluxX[lev] = std::move(fluxX_new);

    std::unique_ptr<MultiFab> ropX_new(new MultiFab(new_grids,new_dmap,1,ropX[lev]->nGrow()));
    ropX_new->copy(*ropX[lev]);
    ropX[lev] = std::move(ropX_new);

    std::unique_ptr<MultiFab> f_gds_u_new(new MultiFab(new_grids,new_dmap,1,f_gds_u[lev]->nGrow()));
    f_gds_u_new->copy(*f_gds_u[lev]);
    f_gds_u[lev] = std::move(f_gds_u_new);

    std::unique_ptr<MultiFab> drag_u_new(new MultiFab(new_grids,new_dmap,1,drag_u[lev]->nGrow()));
    drag_u_new->copy(*drag_u[lev]);
    drag_u[lev] = std::move(drag_u_new);

    // ********************************************************************************
    // Y-face-based arrays
    // ********************************************************************************

    // Create a BoxArray on y-faces.
    BoxArray y_edge_ba = new_grids;
    y_edge_ba.surroundingNodes(1);

    // Y-axis gas velocity
    std::unique_ptr<MultiFab> v_g_new(new MultiFab(new_grids,new_dmap,1,v_g[lev]->nGrow()));
    v_g_new->copy(*v_g[lev]);
    v_g[lev] = std::move(v_g_new);

    std::unique_ptr<MultiFab> v_go_new(new MultiFab(new_grids,new_dmap,1,v_go[lev]->nGrow()));
    v_go_new->copy(*v_go[lev]);
    v_go[lev] = std::move(v_go_new);

    std::unique_ptr<MultiFab> v_gt_new(new MultiFab(new_grids,new_dmap,1,v_gt[lev]->nGrow()));
    v_gt_new->copy(*v_gt[lev]);
    v_gt[lev] = std::move(v_gt_new);

    std::unique_ptr<MultiFab> d_n_new(new MultiFab(new_grids,new_dmap,1,d_n[lev]->nGrow()));
    d_n_new->copy(*d_n[lev]);
    d_n[lev] = std::move(d_n_new);

    std::unique_ptr<MultiFab> tau_v_g_new(new MultiFab(new_grids,new_dmap,1,tau_v_g[lev]->nGrow()));
    tau_v_g_new->copy(*tau_v_g[lev]);
    tau_v_g[lev] = std::move(tau_v_g_new);

    std::unique_ptr<MultiFab> fluxY_new(new MultiFab(new_grids,new_dmap,1,fluxY[lev]->nGrow()));
    fluxY_new->copy(*fluxY[lev]);
    fluxY[lev] = std::move(fluxY_new);

    std::unique_ptr<MultiFab> ropY_new(new MultiFab(new_grids,new_dmap,1,ropY[lev]->nGrow()));
    ropY_new->copy(*ropY[lev]);
    ropY[lev] = std::move(ropY_new);

    std::unique_ptr<MultiFab> f_gds_v_new(new MultiFab(new_grids,new_dmap,1,f_gds_v[lev]->nGrow()));
    f_gds_v_new->copy(*f_gds_v[lev]);
    f_gds_v[lev] = std::move(f_gds_v_new);

    std::unique_ptr<MultiFab> drag_v_new(new MultiFab(new_grids,new_dmap,1,drag_v[lev]->nGrow()));
    drag_v_new->copy(*drag_v[lev]);
    drag_v[lev] = std::move(drag_v_new);

    // ********************************************************************************
    // Z-face-based arrays
    // ********************************************************************************

    // Create a BoxArray on z-faces.
    BoxArray z_edge_ba = new_grids;
    z_edge_ba.surroundingNodes(2);

    // Z-axis gas velocity
    std::unique_ptr<MultiFab> w_g_new(new MultiFab(new_grids,new_dmap,1,w_g[lev]->nGrow()));
    w_g_new->copy(*w_g[lev]);
    w_g[lev] = std::move(w_g_new);

    std::unique_ptr<MultiFab> w_go_new(new MultiFab(new_grids,new_dmap,1,w_go[lev]->nGrow()));
    w_go_new->copy(*w_go[lev]);
    w_go[lev] = std::move(w_go_new);

    std::unique_ptr<MultiFab> w_gt_new(new MultiFab(new_grids,new_dmap,1,w_gt[lev]->nGrow()));
    w_gt_new->copy(*w_gt[lev]);
    w_gt[lev] = std::move(w_gt_new);

    std::unique_ptr<MultiFab> d_t_new(new MultiFab(new_grids,new_dmap,1,d_t[lev]->nGrow()));
    d_t_new->copy(*d_t[lev]);
    d_t[lev] = std::move(d_t_new);

    std::unique_ptr<MultiFab> tau_w_g_new(new MultiFab(new_grids,new_dmap,1,tau_w_g[lev]->nGrow()));
    tau_w_g_new->copy(*tau_w_g[lev]);
    tau_w_g[lev] = std::move(tau_w_g_new);

    std::unique_ptr<MultiFab> fluxZ_new(new MultiFab(new_grids,new_dmap,1,fluxZ[lev]->nGrow()));
    fluxZ_new->copy(*fluxZ[lev]);
    fluxZ[lev] = std::move(fluxZ_new);

    std::unique_ptr<MultiFab> ropZ_new(new MultiFab(new_grids,new_dmap,1,ropZ[lev]->nGrow()));
    ropZ_new->copy(*ropZ[lev]);
    ropZ[lev] = std::move(ropZ_new);

    std::unique_ptr<MultiFab> f_gds_w_new(new MultiFab(new_grids,new_dmap,1,f_gds_w[lev]->nGrow()));
    f_gds_w_new->copy(*f_gds_w[lev]);
    f_gds_w[lev] = std::move(f_gds_w_new);

    std::unique_ptr<MultiFab> drag_w_new(new MultiFab(new_grids,new_dmap,1,drag_w[lev]->nGrow()));
    drag_w_new->copy(*drag_w[lev]);
    drag_w[lev] = std::move(drag_w_new);
}
