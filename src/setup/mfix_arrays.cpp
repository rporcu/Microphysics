#include <mfix.H>
#include <AMReX_EB_utils.H>

void
mfix::ResizeArrays ()
{
    int nlevs_max = maxLevel() + 1;

    ep_g.resize(nlevs_max);
    ep_go.resize(nlevs_max);

    p_g.resize(nlevs_max);
    p_go.resize(nlevs_max);

    p0_g.resize(nlevs_max);

    ro_g.resize(nlevs_max);
    ro_go.resize(nlevs_max);

    trac.resize(nlevs_max);
    trac_o.resize(nlevs_max);

    phi_nd.resize(nlevs_max);
    diveu.resize(nlevs_max);

    // RHS arrays for cell-centered solves
    diff_rhs.resize(nlevs_max);
    diff_rhs1.resize(nlevs_max);
    diff_rhs4.resize(nlevs_max);

    // Solution array for diffusion solves
    diff_phi.resize(nlevs_max);
    diff_phi1.resize(nlevs_max);
    diff_phi4.resize(nlevs_max);

    // MAC velocities at faces
    u_mac.resize(nlevs_max);
    v_mac.resize(nlevs_max);
    w_mac.resize(nlevs_max);

    // RHS array for MAC projection
    mac_rhs.resize(nlevs_max);

    // Solution array for MAC projection
    mac_phi.resize(nlevs_max);

    // Current (vel_g) and old (vel_go) velocities
    vel_g.resize(nlevs_max);
    vel_go.resize(nlevs_max);

    // Pressure gradients
    gp.resize(nlevs_max);

    drag.resize(nlevs_max);

    mu_g.resize(nlevs_max);

    // Vorticity
    vort.resize(nlevs_max);

    xslopes_u.resize(nlevs_max);
    yslopes_u.resize(nlevs_max);
    zslopes_u.resize(nlevs_max);

    xslopes_s.resize(nlevs_max);
    yslopes_s.resize(nlevs_max);
    zslopes_s.resize(nlevs_max);

    bcoeff.resize(nlevs_max);

    // Fuid cost (load balancing)
    fluid_cost.resize(nlevs_max);

    // Fluid grid EB factory
    ebfactory.resize(nlevs_max);

    /****************************************************************************
     *                                                                          *
     * Initialize particle data (including level-set data)                      *
     * NOTE: the level-set data (as well as implicit functions) live on at      *
     *       least two levels                                                   *
     *                                                                          *
     ***************************************************************************/

    // Particle costs (load balancing)
    particle_cost.resize(nlevs_max);

    // Particle grid EB factory
    particle_ebfactory.resize(nlevs_max);

    eb_levels.resize(std::max(2, nlevs_max));
    particle_eb_levels.resize(std::max(2, nlevs_max));

    level_sets.resize(std::max(2, nlevs_max));
}

void
mfix::AllocateArrays (int lev)
{
    if (ooo_debug) amrex::Print() << "AllocateArrays" << std::endl;
    mfix_update_ebfactory(lev);

    // ********************************************************************************
    // Cell- or node-based arrays
    // ********************************************************************************

    // Void fraction
    ep_g[lev] = new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]);
    ep_go[lev] = new MultiFab(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]);
    ep_g[lev]->setVal(1.);
    ep_go[lev]->setVal(1.);

    // Gas density
    ro_g[lev] = new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]);
    ro_go[lev] = new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]);
    ro_g[lev]->setVal(0.);
    ro_go[lev]->setVal(0.);

    // Tracer in gas
    trac[lev] = new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]);
    trac_o[lev] = new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]);
    trac[lev]->setVal(0.);
    trac_o[lev]->setVal(0.);

    const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

    p0_g[lev] = new MultiFab(nd_grids,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]);
    p0_g[lev]->setVal(0.);

    p_g[lev] = new MultiFab(nd_grids,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]);
    p_g[lev]->setVal(0.);

    p_go[lev] = new MultiFab(nd_grids,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]);
    p_go[lev]->setVal(0.);

    phi_nd[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0, MFInfo(), *ebfactory[lev]));
    phi_nd[lev]->setVal(0.);

    diveu[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0, MFInfo(), *ebfactory[lev]));
    diveu[lev]->setVal(0.);

    // Presssure gradients
    gp[lev] = new MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]);
    gp[lev]->setVal(0.);

    // Molecular viscosity
    mu_g[lev] = new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]);
    mu_g[lev]->setVal(0.);

    // Current velocity
    vel_g[lev] = new MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]);
    vel_g[lev]->setVal(0.);

    // Old velocity
    vel_go[lev] = new  MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]);
    vel_go[lev]->setVal(0.);

    // Vorticity
    vort[lev] = new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]);
    vort[lev]->setVal(0.);

    // This is the deposition of the drag force onto the grid
    // 0,1,2 is (drag coefficient * particle velocity)
    // 4 is drag coefficient
    drag[lev] = new  MultiFab(grids[lev],dmap[lev],4,nghost, MFInfo(), *ebfactory[lev]);
    drag[lev]->setVal(0.);

    // Array to store the rhs for diffusion solves -- no ghost cells needed
    diff_rhs[lev].reset(new MultiFab(grids[lev],dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
    diff_rhs[lev]->setVal(0.);

    diff_rhs1[lev].reset(new MultiFab(grids[lev],dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]));
    diff_rhs1[lev]->setVal(0.);

    diff_rhs4[lev].reset(new MultiFab(grids[lev],dmap[lev], 4, 0, MFInfo(), *ebfactory[lev]));
    diff_rhs4[lev]->setVal(0.);

    // Array to store the solution for diffusion solves -- only 1 ghost cell needed
    diff_phi[lev].reset(new MultiFab(grids[lev],dmap[lev], 3, 1, MFInfo(), *ebfactory[lev]));
    diff_phi[lev]->setVal(0.);

    diff_phi1[lev].reset(new MultiFab(grids[lev],dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]));
    diff_phi1[lev]->setVal(0.);

    diff_phi4[lev].reset(new MultiFab(grids[lev],dmap[lev], 4, 1, MFInfo(), *ebfactory[lev]));
    diff_phi4[lev]->setVal(0.);

    // Array to store the rhs for MAC projection
    mac_rhs[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    mac_rhs[lev]->setVal(0.);

    // Array to store the solution for MAC projections
    mac_phi[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    mac_phi[lev]->setVal(0.);

    // Slopes in x-direction
    xslopes_u[lev] = new MultiFab(grids[lev],dmap[lev], 3, nghost, MFInfo(), *ebfactory[lev]);
    xslopes_u[lev]->setVal(0.);
    xslopes_s[lev] = new MultiFab(grids[lev],dmap[lev], 2, nghost, MFInfo(), *ebfactory[lev]);
    xslopes_s[lev]->setVal(0.);

    // Slopes in y-direction
    yslopes_u[lev] = new MultiFab(grids[lev],dmap[lev], 3, nghost, MFInfo(), *ebfactory[lev]);
    yslopes_u[lev]->setVal(0.);
    yslopes_s[lev] = new MultiFab(grids[lev],dmap[lev], 2, nghost, MFInfo(), *ebfactory[lev]);
    yslopes_s[lev]->setVal(0.);

    // Slopes in z-direction
    zslopes_u[lev] = new MultiFab(grids[lev],dmap[lev], 3, nghost, MFInfo(), *ebfactory[lev]);
    zslopes_u[lev]->setVal(0.);
    zslopes_s[lev] = new MultiFab(grids[lev],dmap[lev], 2, nghost, MFInfo(), *ebfactory[lev]);
    zslopes_s[lev]->setVal(0.);

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

    // ****************************************************************

    // Create a BoxArray on x-faces.
    BoxArray x_edge_ba = grids[lev];
    x_edge_ba.surroundingNodes(0);

    // x-face-based coefficient for MAC and diffusive solves
    bcoeff[lev][0].reset(new MultiFab(x_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));
    bcoeff[lev][0]->setVal(0.);

    // U velocity at x-faces (MAC)
    u_mac[lev].reset(new MultiFab(x_edge_ba,dmap[lev],1,2,MFInfo(),*ebfactory[lev]));
    u_mac[lev]->setVal(covered_val); // Covered val as initial value

    // Create a BoxArray on y-faces.
    BoxArray y_edge_ba = grids[lev];
    y_edge_ba.surroundingNodes(1);

    // y-face-based coefficient for MAC and diffusive solves
    bcoeff[lev][1].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));
    bcoeff[lev][1]->setVal(0.);

    // V velocity at y-faces (MAC)
    v_mac[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,2,MFInfo(),*ebfactory[lev]));
    v_mac[lev]->setVal(covered_val);

    // Create a BoxArray on z-faces.
    BoxArray z_edge_ba = grids[lev];
    z_edge_ba.surroundingNodes(2);

    // z-face-based coefficient for MAC and diffusive solves
    bcoeff[lev][2].reset(new MultiFab(z_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));
    bcoeff[lev][2]->setVal(0.);

    // W velocity at z-faces (MAC)
    w_mac[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,2,MFInfo(),*ebfactory[lev]));
    w_mac[lev]->setVal(covered_val);

}

void
mfix::RegridArrays (int lev)
{
    if (ooo_debug) amrex::Print() << "RegridArrays" << std::endl;
    bool need_regrid = mfix_update_ebfactory(lev);

    // exit this function if ebfactory has not been updated because that means
    // that dm and ba haven't changed
    if (!need_regrid)
        return;

    // ********************************************************************************
    // Cell-based arrays
    // ********************************************************************************
    //
    // Note: after calling copy() using dst_ngrow, we do not need to call FillBoundary().
    //
    //

    // Void fraction
    MultiFab* ep_g_new = new MultiFab(grids[lev],dmap[lev], ep_g[lev]->nComp(),
                                      ep_g[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    ep_g_new->setVal(1.);
    ep_g_new->copy(*ep_g[lev], 0, 0, ep_g[lev]->nComp(), ep_g[lev]->nGrow(), ep_g[lev]->nGrow());
    std::swap(ep_g[lev], ep_g_new);
    delete ep_g_new;

    // Old void fraction
    MultiFab* ep_go_new = new MultiFab(grids[lev],dmap[lev], ep_go[lev]->nComp(),
                                       ep_go[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    ep_go_new->setVal(1.);
    ep_go_new->copy(*ep_go[lev], 0, 0, ep_go[lev]->nComp(), ep_go[lev]->nGrow(), ep_go[lev]->nGrow());
    std::swap(ep_go[lev], ep_go_new);
    delete ep_go_new;

    // Gas density
    MultiFab* ro_g_new = new MultiFab(grids[lev],dmap[lev], ro_g[lev]->nComp(),
                                       ro_g[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    ro_g_new->setVal(0.);
    ro_g_new->copy(*ro_g[lev], 0, 0, ro_g[lev]->nComp(), ro_g[lev]->nGrow(), ro_g[lev]->nGrow());
    std::swap(ro_g[lev], ro_g_new);
    delete ro_g_new;

    // Old gas density
    MultiFab* ro_go_new = new MultiFab(grids[lev],dmap[lev], ro_go[lev]->nComp(),
                                       ro_go[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    ro_go_new->setVal(0.);
    ro_go_new->copy(*ro_go[lev], 0, 0, ro_go[lev]->nComp(), ro_go[lev]->nGrow(), ro_go[lev]->nGrow());
    std::swap(ro_go[lev], ro_go_new);
    delete ro_go_new;

    // Tracer in gas
    MultiFab* trac_new = new MultiFab(grids[lev],dmap[lev], trac[lev]->nComp(),
                                       trac[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    trac_new->setVal(0.);
    trac_new->copy(*trac[lev], 0, 0, trac[lev]->nComp(), trac[lev]->nGrow(), trac[lev]->nGrow());
    std::swap(trac[lev], trac_new);
    delete trac_new;

    // Old tracer in gas
    MultiFab* trac_o_new = new MultiFab(grids[lev],dmap[lev], trac_o[lev]->nComp(),
                                       trac_o[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    trac_o_new->setVal(0.);
    trac_o_new->copy(*trac_o[lev], 0, 0, trac_o[lev]->nComp(), trac_o[lev]->nGrow(), trac_o[lev]->nGrow());
    std::swap(trac_o[lev], trac_o_new);
    delete trac_o_new;

    const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

    MultiFab* p_g_new = new MultiFab(grids[lev],dmap[lev], p_g[lev]->nComp(),
                                       p_g[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    p_g_new->setVal(0.);
    p_g_new->copy(*p_g[lev], 0, 0, p_g[lev]->nComp(), p_g[lev]->nGrow(), p_g[lev]->nGrow());
    std::swap(p_g[lev], p_g_new);
    delete p_g_new;

    MultiFab* p_go_new = new MultiFab(grids[lev],dmap[lev], p_go[lev]->nComp(),
                                       p_go[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    p_go_new->setVal(0.);
    p_go_new->copy(*p_go[lev], 0, 0, p_go[lev]->nComp(), p_go[lev]->nGrow(), p_go[lev]->nGrow());
    std::swap(p_go[lev], p_go_new);
    delete p_go_new;

    MultiFab* p0_g_new = new MultiFab(grids[lev],dmap[lev], p0_g[lev]->nComp(),
                                       p0_g[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    p0_g_new->setVal(0.);
    p0_g_new->copy(*p0_g[lev], 0, 0, p0_g[lev]->nComp(), p0_g[lev]->nGrow(), p0_g[lev]->nGrow());
    std::swap(p0_g[lev], p0_g_new);
    delete p0_g_new;

    std::unique_ptr<MultiFab> diveu_new(new MultiFab(nd_grids,dmap[lev],
                                        diveu[lev]->nComp(),diveu[lev]->nGrow(),MFInfo(),*ebfactory[lev]));
    diveu[lev] = std::move(diveu_new);
    diveu[lev]->setVal(0.);

    std::unique_ptr<MultiFab> phi_new(new MultiFab(nd_grids,dmap[lev],
                                      phi_nd[lev]->nComp(),phi_nd[lev]->nGrow(),MFInfo(),*ebfactory[lev]));
    phi_nd[lev] = std::move(phi_new);
    phi_nd[lev]->setVal(0.);

    // Molecular viscosity
    MultiFab* mu_g_new = new MultiFab(grids[lev],dmap[lev], mu_g[lev]->nComp(),
                                       mu_g[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    mu_g_new->setVal(0.);
    mu_g_new->copy(*mu_g[lev], 0, 0, mu_g[lev]->nComp(), mu_g[lev]->nGrow(), mu_g[lev]->nGrow());
    std::swap(mu_g[lev], mu_g_new);
    delete mu_g_new;

    // Gas velocity
    MultiFab* vel_g_new = new MultiFab(grids[lev],dmap[lev], vel_g[lev]->nComp(),
                                       vel_g[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    vel_g_new->setVal(0.);
    vel_g_new->copy(*vel_g[lev], 0, 0, vel_g[lev]->nComp(), vel_g[lev]->nGrow(), vel_g[lev]->nGrow());
    std::swap(vel_g[lev], vel_g_new);
    delete vel_g_new;

    // Old gas velocity
    MultiFab* vel_go_new = new MultiFab(grids[lev],dmap[lev], vel_go[lev]->nComp(),
                                       vel_go[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    vel_go_new->setVal(0.);
    vel_go_new->copy(*vel_go[lev], 0, 0, vel_go[lev]->nComp(), vel_go[lev]->nGrow(), vel_go[lev]->nGrow());
    std::swap(vel_go[lev], vel_go_new);
    delete vel_go_new;

    // Pressure gradients
    MultiFab* gp_new = new MultiFab(grids[lev],dmap[lev], gp[lev]->nComp(),
                                       gp[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    gp_new->setVal(0.);
    gp_new->copy(*gp[lev], 0, 0, gp[lev]->nComp(), gp[lev]->nGrow(), gp[lev]->nGrow());
    std::swap(gp[lev], gp_new);
    delete gp_new;

    // Vorticity
    MultiFab* vort_new = new MultiFab(grids[lev],dmap[lev], vort[lev]->nComp(),
                                       vort[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    vort_new->setVal(0.);
    vort_new->copy(*vort[lev], 0, 0, vort[lev]->nComp(), vort[lev]->nGrow(), vort[lev]->nGrow());
    std::swap(vort[lev], vort_new);
    delete vort_new;

    // Particle/fluid drag -- note it is important to copy from previous step in order to use in dt calculation
    MultiFab* drag_new = new MultiFab(grids[lev],dmap[lev], drag[lev]->nComp(),
                                       drag[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    drag_new->setVal(0.);
    drag_new->copy(*drag[lev], 0, 0, drag[lev]->nComp(), drag[lev]->nGrow(), drag[lev]->nGrow());
    std::swap(drag[lev], drag_new);
    delete drag_new;

    // Array to store the rhs for tensor diffusion solve
    std::unique_ptr<MultiFab> diff_rhs_new(new  MultiFab(grids[lev], dmap[lev],
                                           diff_rhs[lev]->nComp(),diff_rhs[lev]->nGrow(),MFInfo(),*ebfactory[lev]));
    diff_rhs[lev] = std::move(diff_rhs_new);
    diff_rhs[lev] -> setVal(0.);

    // Array to store the rhs for tensor diffusion solve
    std::unique_ptr<MultiFab> diff_rhs1_new(new  MultiFab(grids[lev], dmap[lev],
                                                          diff_rhs1[lev]->nComp(), diff_rhs1[lev]->nGrow(),
                                                          MFInfo(), *ebfactory[lev]));
    diff_rhs1[lev] = std::move(diff_rhs1_new);
    diff_rhs1[lev] -> setVal(0.);

    // // Array to store the rhs for tensor diffusion solve
    std::unique_ptr<MultiFab> diff_rhs4_new(new  MultiFab(grids[lev], dmap[lev],
                                                          diff_rhs4[lev]->nComp(), diff_rhs4[lev]->nGrow(),
                                                          MFInfo(), *ebfactory[lev]));
    diff_rhs4[lev] = std::move(diff_rhs4_new);
    diff_rhs4[lev] -> setVal(0.);

    // Arrays to store the solution for diffusion solves
    std::unique_ptr<MultiFab> diff_phi_new(new  MultiFab(grids[lev], dmap[lev],
                                           diff_phi[lev]->nComp(),diff_phi[lev]->nGrow(),MFInfo(),*ebfactory[lev]));
    diff_phi[lev] = std::move(diff_phi_new);
    diff_phi[lev] -> setVal(0.);

    // Arrays to store the solution for diffusion solves
    std::unique_ptr<MultiFab> diff_phi1_new(new  MultiFab(grids[lev], dmap[lev],
                                                          diff_phi1[lev]->nComp(), diff_phi1[lev]->nGrow(),
                                                          MFInfo(), *ebfactory[lev]));
    diff_phi1[lev] = std::move(diff_phi1_new);
    diff_phi1[lev] -> setVal(0.);

    // Arrays to store the solution for diffusion solves
    std::unique_ptr<MultiFab> diff_phi4_new(new  MultiFab(grids[lev], dmap[lev],
                                                          diff_phi4[lev]->nComp(), diff_phi4[lev]->nGrow(),
                                                          MFInfo(), *ebfactory[lev]));
    diff_phi4[lev] = std::move(diff_phi4_new);
    diff_phi4[lev] -> setVal(0.);

    // Array to store the rhs for cell-centered solves
    std::unique_ptr<MultiFab> mac_rhs_new(new  MultiFab(grids[lev], dmap[lev],
                                          mac_rhs[lev]->nComp(),mac_rhs[lev]->nGrow(),MFInfo(),*ebfactory[lev]));
    mac_rhs[lev] = std::move(mac_rhs_new);
    mac_rhs[lev] -> setVal(0.);

    // Arrays to store the solution for the MAC projection
    std::unique_ptr<MultiFab> mac_phi_new(new  MultiFab(grids[lev], dmap[lev],
                                          mac_phi[lev]->nComp(),mac_phi[lev]->nGrow(),MFInfo(),*ebfactory[lev]));
    mac_phi[lev] = std::move(mac_phi_new);
    mac_phi[lev] -> setVal(0.);

    // Slopes in x-direction
    MultiFab* xslopes_u_new = new MultiFab(grids[lev], dmap[lev],
                                           xslopes_u[lev]->nComp(),
                                           xslopes_u[lev]->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    xslopes_u_new->setVal(0.);
    std::swap(xslopes_u[lev], xslopes_u_new);
    delete xslopes_u_new;

    MultiFab* xslopes_s_new = new MultiFab(grids[lev], dmap[lev],
                                           xslopes_s[lev]->nComp(),
                                           xslopes_s[lev]->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    xslopes_s_new->setVal(0.);
    std::swap(xslopes_s[lev], xslopes_s_new);
    delete xslopes_s_new;

    // Slopes in y-direction
    MultiFab* yslopes_u_new = new MultiFab(grids[lev], dmap[lev],
                                           yslopes_u[lev]->nComp(),
                                           yslopes_u[lev]->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    yslopes_u_new->setVal(0.);
    std::swap(yslopes_u[lev], yslopes_u_new);
    delete yslopes_u_new;

    MultiFab* yslopes_s_new = new MultiFab(grids[lev], dmap[lev],
                                           yslopes_s[lev]->nComp(),
                                           yslopes_s[lev]->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    yslopes_s_new->setVal(0.);
    std::swap(yslopes_s[lev], yslopes_s_new);
    delete yslopes_s_new;

    // Slopes in z-direction
    MultiFab* zslopes_u_new = new MultiFab(grids[lev], dmap[lev],
                                           zslopes_u[lev]->nComp(),
                                           zslopes_u[lev]->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    zslopes_u_new->setVal(0.);
    std::swap(zslopes_u[lev], zslopes_u_new);
    delete zslopes_u_new;

    MultiFab* zslopes_s_new = new MultiFab(grids[lev], dmap[lev],
                                           zslopes_s[lev]->nComp(),
                                           zslopes_s[lev]->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    zslopes_s_new->setVal(0.);
    std::swap(zslopes_s[lev], zslopes_s_new);
    delete zslopes_s_new;

   /****************************************************************************
    * x-face-based arrays                                                        *
    ****************************************************************************/

    BoxArray x_ba = grids[lev];
    x_ba = x_ba.surroundingNodes(0);

    // MAC/diffusion coefficient on x-faces
    std::unique_ptr<MultiFab> bcx_mac_new(new MultiFab(x_ba,dmap[lev],
                                          bcoeff[lev][0]->nComp(),bcoeff[lev][0]->nGrow(),MFInfo(),*ebfactory[lev]));
    bcoeff[lev][0] = std::move(bcx_mac_new);
    bcoeff[lev][0] -> setVal(0.0);

    // x component of MAC velocity
    std::unique_ptr<MultiFab> u_mac_new(new MultiFab(x_ba,dmap[lev],
                                        u_mac[lev]->nComp(),u_mac[lev]->nGrow(),MFInfo(),*ebfactory[lev]));
    u_mac[lev] = std::move(u_mac_new);
    u_mac[lev] -> setVal(0.0);

   /****************************************************************************
    * y-face-based arrays                                                        *
    ****************************************************************************/

    BoxArray y_ba = grids[lev];
    y_ba = y_ba.surroundingNodes(1);

    // MAC/diffusion coefficient on y-faces
    std::unique_ptr<MultiFab> bcy_mac_new(new MultiFab(y_ba,dmap[lev],
                                          bcoeff[lev][1]->nComp(),bcoeff[lev][1]->nGrow(),MFInfo(),*ebfactory[lev]));
    bcoeff[lev][1] = std::move(bcy_mac_new);
    bcoeff[lev][1] -> setVal(0.0);

    // y component of MAC velocity
    std::unique_ptr<MultiFab> v_mac_new(new MultiFab(y_ba,dmap[lev],
                                        v_mac[lev]->nComp(),v_mac[lev]->nGrow(),MFInfo(),*ebfactory[lev]));
    v_mac[lev] = std::move(v_mac_new);
    v_mac[lev] -> setVal(0.0);

   /****************************************************************************
    * z-face-based arrays                                                        *
    ****************************************************************************/

    BoxArray z_ba = grids[lev];
    z_ba = z_ba.surroundingNodes(2);

    // MAC/diffusion coefficient on z-faces
    std::unique_ptr<MultiFab> bcz_mac_new(new MultiFab(z_ba,dmap[lev],
                                          bcoeff[lev][2]->nComp(),bcoeff[lev][2]->nGrow(),MFInfo(),*ebfactory[lev]));
    bcoeff[lev][2] = std::move(bcz_mac_new);
    bcoeff[lev][2] -> setVal(0.0);

    // z component of MAC velocity
    std::unique_ptr<MultiFab> w_mac_new(new MultiFab(z_ba,dmap[lev],
                                        w_mac[lev]->nComp(),w_mac[lev]->nGrow(),MFInfo(),*ebfactory[lev]));
    w_mac[lev] = std::move(w_mac_new);
    w_mac[lev] -> setVal(0.0);
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
    if (ooo_debug) amrex::Print() << "RegridLevelSetArray" << std::endl;
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
          new EBFArrayBoxFactory(* particle_eb_levels[a_lev], geom[a_lev], ba, dm,
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
              new EBFArrayBoxFactory( * particle_eb_levels[a_lev], geom[a_lev], ba, dm,
                                      {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                       m_eb_full_grow_cells}, m_eb_support_level)
              );

         changed = true;
      }
   }

   if (changed)
   {

       Print() << "Regridding level-set on lev = " << a_lev << std::endl;

       const BoxArray nd_ba = amrex::convert(ba, IntVect::TheNodeVector());

       MultiFab* new_level_set = new MultiFab();

       if (level_sets[a_lev]->boxArray() == nd_ba)
       {
           MFUtil::regrid(*new_level_set, nd_ba, dm, *level_sets[a_lev], true);
       }
       else
       {
           int nc = level_sets[a_lev]->nComp();
           int ng = level_sets[a_lev]->nGrow();
           const Periodicity& period = geom[a_lev].periodicity();
           new_level_set->define(nd_ba, dm, nc, ng);
           new_level_set->setVal(0.);
           new_level_set->copy(*level_sets[a_lev], 0, 0, nc, 0, ng, period);
       }

       std::swap(level_sets[a_lev], new_level_set);
       delete new_level_set;

       //________________________________________________________________________
       // If we're operating in single-level mode, the level-set has a second
       // (refined) MultiFab that also needs to be regridded.
       if ((nlev == 1) && (a_lev == 0))
       {
           Print() << "Also regridding refined level-set" << std::endl;

           BoxArray ref_nd_ba = amrex::convert(ba, IntVect::TheNodeVector());
           ref_nd_ba.refine(levelset_refinement);

           MultiFab* new_level_set = new MultiFab();

           if (level_sets[a_lev+1]->boxArray() == ref_nd_ba)
           {
               MFUtil::regrid(*new_level_set, ref_nd_ba, dm, *level_sets[a_lev+1], true);
           }
           else
           {
               int nc = level_sets[a_lev+1]->nComp();
               int ng = level_sets[a_lev+1]->nGrow();
               const Periodicity& period = geom[a_lev].periodicity();
               new_level_set->define(ref_nd_ba, dm, nc, ng);
               new_level_set->setVal(0.0);
               new_level_set->copy(*level_sets[a_lev+1], 0, 0, nc, 0, ng, period);
           }

           std::swap(level_sets[a_lev+1], new_level_set);
           delete new_level_set;
       }
   }
}


//!  This function checks if ebfactory is allocated with the proper dm and ba
bool mfix::mfix_update_ebfactory (int a_lev)
{
    if (ooo_debug) amrex::Print() << "mfix_update_ebfactory" << std::endl;
   // This assert is to verify that some kind of EB geometry has already been
   // defined
   AMREX_ASSERT(not EB2::IndexSpace::empty());

   const DistributionMapping & dm = DistributionMap(a_lev);
   const BoxArray &            ba = boxArray(a_lev);

   bool is_updated = false;

   if ( ebfactory[a_lev].get() == nullptr )
   {
      Print() << "Updating ebfactory from nullptr" << std::endl;

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
          Print() << "Updating ebfactory from existing" << std::endl;

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
