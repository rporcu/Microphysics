#include <AMReX_ParmParse.H>

#include <mfix_F.H>
#include <mfix_level.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

void
mfix_level::Regrid (int lev, int nstep, int dual_grid)
{
    amrex::Print() << "In Regrid at step " << nstep << std::endl;

    if (load_balance_type == "KDTree")
    {
       // This creates a new BA and new DM, re-defines the particle BA and DM to be these new ones,
       //      and calls Redistribute.  This doesn't touch the fluid grids.
       pc -> BalanceParticleLoad_KDTree ();

       if (dual_grid == 0)
       {
          SetBoxArray(lev, pc->ParticleBoxArray(lev));
          SetDistributionMap(lev, pc->ParticleDistributionMap(lev));

          // Since we have already allocated the fluid data we need to re-define those arrays
          //   and copy from the old BoxArray to the new one.  Note that the SetBoxArray and
          //   SetDistributionMap calls above have re-defined grids and dmap to be the new ones.
          if (solve_fluid) 
             RegridArrays(lev,grids[lev],dmap[lev]);
       }

       mfix_set_bc0(lev);
    }
    else if (load_balance_type == "KnapSack") {

        amrex::Print() << "Load balancing using KnapSack " << std::endl;
        
        if (ParallelDescriptor::NProcs() == 1) return;
        
        AMREX_ALWAYS_ASSERT(costs[0] != nullptr);

        for (int lev = 0; lev <= finestLevel(); ++lev)
        {
            DistributionMapping newdm = DistributionMapping::makeKnapSack(*costs[lev]);
            RegridArrays(lev, grids[lev], newdm);
            SetDistributionMap(lev, newdm);            
        }
        
        pc->Regrid(dmap[lev], grids[lev]);

        mfix_set_bc0(lev);
    }
}

void
mfix_level::RegridOnRestart (int lev)
{
    amrex::Print() << "In RegridOnRestart " << std::endl;

    if (load_balance_type == "FixedSize" || load_balance_type == "KnapSack")
    {
       // We hold on to the old_ba so that we can test the new BoxArray against it
       //   to see if the grids have changed
       BoxArray old_ba(grids[lev]);

       // This creates a new BoxArray (based on the new max_grid_size)
       const BoxArray& ba = MakeBaseGrids();

       // This creates the associated Distribution Mapping
       DistributionMapping dm(ba, ParallelDescriptor::NProcs());

       // This sets grids[lev] = ba 
       SetBoxArray(lev, ba);

       // This sets dmap[lev] = dm
       SetDistributionMap(lev, dm);

       // If the grids have changed, we need to re-define those arrays 
       // and copy from the old BoxArray to the new one since we have already read in the old fluid data
       if ( (grids[0] != old_ba) && solve_fluid) 
          RegridArrays(lev,grids[lev],dmap[lev]);

       mfix_set_bc0(lev);

       pc->BuildLevelMask(lev, geom[lev], dm, ba); 
       pc->Redistribute();

    }
}

void
mfix_level::RegridArrays (int lev, BoxArray& new_grids, DistributionMapping& new_dmap)
{
    // ********************************************************************************
    // Cell-based arrays
    // ********************************************************************************

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

    // Pressure correction equation
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

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

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

    ng = d_e[lev]->nGrow();
    std::unique_ptr<MultiFab> d_e_new(new MultiFab(x_edge_ba,new_dmap,1,ng));
    d_e_new->copy(*d_e[lev],0,0,1,ng,ng);
    d_e_new->FillBoundary(geom[lev].periodicity());
    d_e[lev] = std::move(d_e_new);

    ng = tau_u_g[lev]->nGrow();
    std::unique_ptr<MultiFab> tau_u_g_new(new MultiFab(x_edge_ba,new_dmap,1,ng));
    tau_u_g_new->copy(*tau_u_g[lev],0,0,1,ng,ng);
    tau_u_g_new->FillBoundary(geom[lev].periodicity());
    tau_u_g[lev] = std::move(tau_u_g_new);

    ng = fluxX[lev]->nGrow();
    std::unique_ptr<MultiFab> fluxX_new(new MultiFab(x_edge_ba,new_dmap,1,ng));
    fluxX_new->copy(*fluxX[lev],0,0,1,ng,ng);
    fluxX_new->FillBoundary(geom[lev].periodicity());
    fluxX[lev] = std::move(fluxX_new);

    ng = ropX[lev]->nGrow();
    std::unique_ptr<MultiFab> ropX_new(new MultiFab(x_edge_ba,new_dmap,1,ng));
    ropX_new->copy(*ropX[lev],0,0,1,ng,ng);
    ropX_new->FillBoundary(geom[lev].periodicity());
    ropX[lev] = std::move(ropX_new);

    ng = f_gds_u[lev]->nGrow();
    std::unique_ptr<MultiFab> f_gds_u_new(new MultiFab(x_edge_ba,new_dmap,1,ng));
    f_gds_u_new->copy(*f_gds_u[lev],0,0,1,ng,ng);
    f_gds_u_new->FillBoundary(geom[lev].periodicity());
    f_gds_u[lev] = std::move(f_gds_u_new);

    ng = drag_u[lev]->nGrow();
    std::unique_ptr<MultiFab> drag_u_new(new MultiFab(x_edge_ba,new_dmap,1,ng));
    drag_u_new->copy(*drag_u[lev],0,0,1,ng,ng);
    drag_u_new->FillBoundary(geom[lev].periodicity());
    drag_u[lev] = std::move(drag_u_new);

    // ********************************************************************************
    // Y-face-based arrays
    // ********************************************************************************

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

    ng = d_n[lev]->nGrow();
    std::unique_ptr<MultiFab> d_n_new(new MultiFab(y_edge_ba,new_dmap,1,ng));
    d_n_new->copy(*d_n[lev],0,0,1,ng,ng);
    d_n_new->FillBoundary(geom[lev].periodicity());
    d_n[lev] = std::move(d_n_new);

    ng = tau_v_g[lev]->nGrow();
    std::unique_ptr<MultiFab> tau_v_g_new(new MultiFab(y_edge_ba,new_dmap,1,ng));
    tau_v_g_new->copy(*tau_v_g[lev],0,0,1,ng,ng);
    tau_v_g_new->FillBoundary(geom[lev].periodicity());
    tau_v_g[lev] = std::move(tau_v_g_new);

    ng = fluxY[lev]->nGrow();
    std::unique_ptr<MultiFab> fluxY_new(new MultiFab(y_edge_ba,new_dmap,1,ng));
    fluxY_new->copy(*fluxY[lev],0,0,1,ng,ng);
    fluxY_new->FillBoundary(geom[lev].periodicity());
    fluxY[lev] = std::move(fluxY_new);

    ng = ropY[lev]->nGrow();
    std::unique_ptr<MultiFab> ropY_new(new MultiFab(y_edge_ba,new_dmap,1,ng));
    ropY_new->copy(*ropY[lev],0,0,1,ng,ng);
    ropY_new->FillBoundary(geom[lev].periodicity());
    ropY[lev] = std::move(ropY_new);

    ng = f_gds_v[lev]->nGrow();
    std::unique_ptr<MultiFab> f_gds_v_new(new MultiFab(y_edge_ba,new_dmap,1,ng));
    f_gds_v_new->copy(*f_gds_v[lev],0,0,1,ng,ng);
    f_gds_v_new->FillBoundary(geom[lev].periodicity());
    f_gds_v[lev] = std::move(f_gds_v_new);

    ng = drag_v[lev]->nGrow();
    std::unique_ptr<MultiFab> drag_v_new(new MultiFab(y_edge_ba,new_dmap,1,ng));
    drag_v_new->copy(*drag_v[lev],0,0,1,ng,ng);
    drag_v_new->FillBoundary(geom[lev].periodicity());
    drag_v[lev] = std::move(drag_v_new);

    // ********************************************************************************
    // Z-face-based arrays
    // ********************************************************************************

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

    ng = d_t[lev]->nGrow();
    std::unique_ptr<MultiFab> d_t_new(new MultiFab(z_edge_ba,new_dmap,1,ng));
    d_t_new->copy(*d_t[lev],0,0,1,ng,ng);
    d_t_new->FillBoundary(geom[lev].periodicity());
    d_t[lev] = std::move(d_t_new);

    ng = tau_w_g[lev]->nGrow();
    std::unique_ptr<MultiFab> tau_w_g_new(new MultiFab(z_edge_ba,new_dmap,1,ng));
    tau_w_g_new->copy(*tau_w_g[lev],0,0,1,ng,ng);
    tau_w_g_new->FillBoundary(geom[lev].periodicity());
    tau_w_g[lev] = std::move(tau_w_g_new);

    ng = fluxZ[lev]->nGrow();
    std::unique_ptr<MultiFab> fluxZ_new(new MultiFab(z_edge_ba,new_dmap,1,ng));
    fluxZ_new->copy(*fluxZ[lev],0,0,1,ng,ng);
    fluxZ_new->FillBoundary(geom[lev].periodicity());
    fluxZ[lev] = std::move(fluxZ_new);

    ng = ropZ[lev]->nGrow();
    std::unique_ptr<MultiFab> ropZ_new(new MultiFab(z_edge_ba,new_dmap,1,ng));
    ropZ_new->copy(*ropZ[lev],0,0,1,ng,ng);
    ropZ_new->FillBoundary(geom[lev].periodicity());
    ropZ[lev] = std::move(ropZ_new);

    ng = f_gds_w[lev]->nGrow();
    std::unique_ptr<MultiFab> f_gds_w_new(new MultiFab(z_edge_ba,new_dmap,1,ng));
    f_gds_w_new->copy(*f_gds_w[lev],0,0,1,ng,ng);
    f_gds_w_new->FillBoundary(geom[lev].periodicity());
    f_gds_w[lev] = std::move(f_gds_w_new);

    ng = drag_w[lev]->nGrow();
    std::unique_ptr<MultiFab> drag_w_new(new MultiFab(z_edge_ba,new_dmap,1,ng));
    drag_w_new->copy(*drag_w[lev],0,0,1,ng,ng);
    drag_w_new->FillBoundary(geom[lev].periodicity());
    drag_w[lev] = std::move(drag_w_new);

    // ********************************************************************************
    // Make sure we fill the ghost cells as appropriate -- this is copied from init_fluid
    // ********************************************************************************

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

    // ********************************************************************************
    // EB stuff
    // ********************************************************************************
    int m_eb_basic_grow_cells = 2;
    int m_eb_volume_grow_cells = 2;
    int m_eb_full_grow_cells = 2;
    EBSupport m_eb_support_level = EBSupport::full;

    if (ebfactory)
       ebfactory.reset(new EBFArrayBoxFactory(geom[lev], new_grids, new_dmap,
                    {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells}, m_eb_support_level));

    if (costs[lev] != nullptr) {
        costs[lev].reset(new MultiFab(new_grids, new_dmap, 1, 0));
        costs[lev]->setVal(0.0);
    }
}
