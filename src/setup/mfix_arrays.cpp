#include <mfix.H>
#include <AMReX_EB_utils.H>

void
mfix::ResizeArrays ()
{
    int nlevs_max = maxLevel() + 1;

    m_leveldata.resize(nlevs_max);
    for (int lev(0); lev < nlevs_max; ++lev)
      m_leveldata[lev].reset(new LevelData());

    phi_nd.resize(nlevs_max);

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
    m_leveldata[lev]->ep_g = new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                          MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->ep_go = new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                           MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->ep_g->setVal(1);
    m_leveldata[lev]->ep_go->setVal(1);

    // Gas density
    m_leveldata[lev]->ro_g = new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                          MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->ro_go = new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                           MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->ro_g->setVal(0);
    m_leveldata[lev]->ro_go->setVal(0);

    // Tracer in gas
    m_leveldata[lev]->trac = new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                          MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->trac_o = new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                            MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->trac->setVal(0);
    m_leveldata[lev]->trac_o->setVal(0);

    const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

    m_leveldata[lev]->p0_g = new MultiFab(nd_grids, dmap[lev], 1, nghost,
                                          MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->p0_g->setVal(0.);

    m_leveldata[lev]->p_g = new MultiFab(nd_grids, dmap[lev], 1, nghost,
                                         MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->p_g->setVal(0.);

    m_leveldata[lev]->p_go = new MultiFab(nd_grids, dmap[lev], 1, nghost,
                                          MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->p_go->setVal(0.);

    phi_nd[lev] = new MultiFab(nd_grids, dmap[lev], 1, 0,
                               MFInfo(), *ebfactory[lev]);
    phi_nd[lev]->setVal(0.);

    m_leveldata[lev]->diveu = new MultiFab(nd_grids, dmap[lev], 1, 0,
                                           MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->diveu->setVal(0);

    // Presssure gradients
    m_leveldata[lev]->gp = new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                        MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->gp->setVal(0.);

    // Molecular viscosity
    m_leveldata[lev]->mu_g = new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                          MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->mu_g->setVal(0.);

    // Current velocity
    m_leveldata[lev]->vel_g = new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                           MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->vel_g->setVal(0);

    // Old velocity
    m_leveldata[lev]->vel_go = new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                            MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->vel_go->setVal(0.);

    // Vorticity
    m_leveldata[lev]->vort = new MultiFab(grids[lev], dmap[lev], 1, nghost,
                                          MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->vort->setVal(0.);

    // This is the deposition of the drag force onto the grid
    // 0,1,2 is (drag coefficient * particle velocity)
    // 4 is drag coefficient
    m_leveldata[lev]->drag = new MultiFab(grids[lev], dmap[lev], 4, nghost,
                                          MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->drag->setVal(0.);

    // Array to store the rhs for diffusion solves -- no ghost cells needed
    diff_rhs[lev] = new MultiFab(grids[lev],dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
    diff_rhs[lev]->setVal(0.);

    diff_rhs1[lev] = new MultiFab(grids[lev],dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
    diff_rhs1[lev]->setVal(0.);

    diff_rhs4[lev] = new MultiFab(grids[lev],dmap[lev], 4, 0, MFInfo(), *ebfactory[lev]);
    diff_rhs4[lev]->setVal(0.);

    // Array to store the solution for diffusion solves -- only 1 ghost cell needed
    diff_phi[lev] = new MultiFab(grids[lev],dmap[lev], 3, 1, MFInfo(), *ebfactory[lev]);
    diff_phi[lev]->setVal(0.);

    diff_phi1[lev] = new MultiFab(grids[lev],dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
    diff_phi1[lev]->setVal(0.);

    diff_phi4[lev] = new MultiFab(grids[lev],dmap[lev], 4, 1, MFInfo(), *ebfactory[lev]);
    diff_phi4[lev]->setVal(0.);

    // Array to store the rhs for MAC projection
    mac_rhs[lev] = new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]);
    mac_rhs[lev]->setVal(0.);

    // Array to store the solution for MAC projections
    mac_phi[lev] = new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]);
    mac_phi[lev]->setVal(0.);

    // Slopes in x-direction
    m_leveldata[lev]->xslopes_u = new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                               MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->xslopes_u->setVal(0.);

    m_leveldata[lev]->xslopes_s = new MultiFab(grids[lev], dmap[lev], 2, nghost,
                                               MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->xslopes_s->setVal(0);

    // Slopes in y-direction
    m_leveldata[lev]->yslopes_u = new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                               MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->yslopes_u->setVal(0.);

    m_leveldata[lev]->yslopes_s = new MultiFab(grids[lev], dmap[lev], 2, nghost,
                                               MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->yslopes_s->setVal(0.);

    // Slopes in z-direction
    m_leveldata[lev]->zslopes_u = new MultiFab(grids[lev], dmap[lev], 3, nghost,
                                  MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->zslopes_u->setVal(0.);

    m_leveldata[lev]->zslopes_s = new MultiFab(grids[lev], dmap[lev], 2, nghost,
                                  MFInfo(), *ebfactory[lev]);
    m_leveldata[lev]->zslopes_s->setVal(0.);

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

    // ****************************************************************

    // Create a BoxArray on x-faces.
    BoxArray x_edge_ba = grids[lev];
    x_edge_ba.surroundingNodes(0);

    // x-face-based coefficient for MAC and diffusive solves
    bcoeff[lev][0] = new MultiFab(x_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]);
    bcoeff[lev][0]->setVal(0.);

    // U velocity at x-faces (MAC)
    u_mac[lev] = new MultiFab(x_edge_ba,dmap[lev],1,2,MFInfo(),*ebfactory[lev]);
    u_mac[lev]->setVal(covered_val); // Covered val as initial value

    // Create a BoxArray on y-faces.
    BoxArray y_edge_ba = grids[lev];
    y_edge_ba.surroundingNodes(1);

    // y-face-based coefficient for MAC and diffusive solves
    bcoeff[lev][1] = new MultiFab(y_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]);
    bcoeff[lev][1]->setVal(0.);

    // V velocity at y-faces (MAC)
    v_mac[lev] = new MultiFab(y_edge_ba,dmap[lev],1,2,MFInfo(),*ebfactory[lev]);
    v_mac[lev]->setVal(covered_val);

    // Create a BoxArray on z-faces.
    BoxArray z_edge_ba = grids[lev];
    z_edge_ba.surroundingNodes(2);

    // z-face-based coefficient for MAC and diffusive solves
    bcoeff[lev][2] = new MultiFab(z_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]);
    bcoeff[lev][2]->setVal(0.);

    // W velocity at z-faces (MAC)
    w_mac[lev] = new MultiFab(z_edge_ba,dmap[lev],1,2,MFInfo(),*ebfactory[lev]);
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
    MultiFab& ep_g = *(m_leveldata[lev]->ep_g);
    MultiFab* ep_g_new = new MultiFab(grids[lev], dmap[lev], ep_g.nComp(),
                                      ep_g.nGrow(), MFInfo(), *ebfactory[lev]);
    ep_g_new->setVal(1);
    ep_g_new->copy(ep_g, 0, 0, ep_g.nComp(), ep_g.nGrow(), ep_g.nGrow());
    std::swap((m_leveldata[lev]->ep_g), ep_g_new);
    delete ep_g_new;

    // Old void fraction
    MultiFab& ep_go = *(m_leveldata[lev]->ep_go);
    MultiFab* ep_go_new = new MultiFab(grids[lev], dmap[lev], ep_go.nComp(),
                                      ep_go.nGrow(), MFInfo(), *ebfactory[lev]);
    ep_go_new->setVal(1);
    ep_go_new->copy(ep_go, 0, 0, ep_go.nComp(), ep_go.nGrow(), ep_go.nGrow());
    std::swap((m_leveldata[lev]->ep_go), ep_go_new);
    delete ep_go_new;

    // Gas density
    MultiFab* ro_g_new = new MultiFab(grids[lev], dmap[lev],
                                      m_leveldata[lev]->ro_g->nComp(),
                                      m_leveldata[lev]->ro_g->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    ro_g_new->setVal(0);
    ro_g_new->copy(*m_leveldata[lev]->ro_g, 0, 0, m_leveldata[lev]->ro_g->nComp(),
                   m_leveldata[lev]->ro_g->nGrow(), m_leveldata[lev]->ro_g->nGrow());
    std::swap(m_leveldata[lev]->ro_g, ro_g_new);
    delete ro_g_new;

    // Old gas density
    MultiFab* ro_go_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->ro_go->nComp(),
                                       m_leveldata[lev]->ro_go->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
    ro_go_new->setVal(0);
    ro_go_new->copy(*m_leveldata[lev]->ro_go, 0, 0, m_leveldata[lev]->ro_go->nComp(),
                    m_leveldata[lev]->ro_go->nGrow(), m_leveldata[lev]->ro_go->nGrow());
    std::swap(m_leveldata[lev]->ro_go, ro_go_new);
    delete ro_go_new;

    // Tracer in gas
    MultiFab* trac_new = new MultiFab(grids[lev], dmap[lev],
                                      m_leveldata[lev]->trac->nComp(),
                                      m_leveldata[lev]->trac->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    trac_new->setVal(0);
    trac_new->copy(*m_leveldata[lev]->trac, 0, 0, m_leveldata[lev]->trac->nComp(),
                   m_leveldata[lev]->trac->nGrow(), m_leveldata[lev]->trac->nGrow());
    std::swap(m_leveldata[lev]->trac, trac_new);
    delete trac_new;

    // Old tracer in gas
    MultiFab* trac_o_new = new MultiFab(grids[lev], dmap[lev],
                                        m_leveldata[lev]->trac_o->nComp(),
                                        m_leveldata[lev]->trac_o->nGrow(),
                                        MFInfo(), *ebfactory[lev]);
    trac_o_new->setVal(0);
    trac_o_new->copy(*m_leveldata[lev]->trac_o, 0, 0, m_leveldata[lev]->trac_o->nComp(),
                     m_leveldata[lev]->trac_o->nGrow(), m_leveldata[lev]->trac_o->nGrow());
    std::swap(m_leveldata[lev]->trac_o, trac_o_new);
    delete trac_o_new;

    const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

    MultiFab* p_g_new = new MultiFab(nd_grids, dmap[lev],
                                     m_leveldata[lev]->p_g->nComp(),
                                     m_leveldata[lev]->p_g->nGrow(),
                                     MFInfo(), *ebfactory[lev]);
    p_g_new->setVal(0);
    p_g_new->copy(*m_leveldata[lev]->p_g, 0, 0, m_leveldata[lev]->p_g->nComp(),
                  m_leveldata[lev]->p_g->nGrow(), m_leveldata[lev]->p_g->nGrow());
    std::swap(m_leveldata[lev]->p_g, p_g_new);
    delete p_g_new;

    MultiFab* p_go_new = new MultiFab(nd_grids, dmap[lev],
                                      m_leveldata[lev]->p_go->nComp(),
                                      m_leveldata[lev]->p_go->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    p_go_new->setVal(0);
    p_go_new->copy(*m_leveldata[lev]->p_go, 0, 0, m_leveldata[lev]->p_go->nComp(),
                   m_leveldata[lev]->p_go->nGrow(), m_leveldata[lev]->p_go->nGrow());
    std::swap(m_leveldata[lev]->p_go, p_go_new);
    delete p_go_new;

    MultiFab* p0_g_new = new MultiFab(nd_grids, dmap[lev],
                                      m_leveldata[lev]->p0_g->nComp(),
                                      m_leveldata[lev]->p0_g->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    p0_g_new->setVal(0);
    p0_g_new->copy(*m_leveldata[lev]->p0_g, 0, 0,
                   m_leveldata[lev]->p0_g->nComp(),
                   m_leveldata[lev]->p0_g->nGrow(),
                   m_leveldata[lev]->p0_g->nGrow());
    std::swap(m_leveldata[lev]->p0_g, p0_g_new);
    delete p0_g_new;

    MultiFab* diveu_new = new MultiFab(nd_grids, dmap[lev], m_leveldata[lev]->diveu->nComp(),
                                       m_leveldata[lev]->diveu->nGrow(), MFInfo(), *ebfactory[lev]);
    diveu_new->setVal(0);
    diveu_new->copy(*m_leveldata[lev]->diveu, 0, 0, m_leveldata[lev]->diveu->nComp(), m_leveldata[lev]->diveu->nGrow(), m_leveldata[lev]->diveu->nGrow());
    std::swap(m_leveldata[lev]->diveu, diveu_new);
    delete diveu_new;

    MultiFab* phi_new = new MultiFab(nd_grids, dmap[lev], phi_nd[lev]->nComp(),
                                     phi_nd[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    phi_new->setVal(0);
    std::swap(phi_nd[lev], phi_new);
    delete phi_new;

    // Molecular viscosity
    MultiFab* mu_g_new = new MultiFab(grids[lev], dmap[lev],
                                      m_leveldata[lev]->mu_g->nComp(),
                                      m_leveldata[lev]->mu_g->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    mu_g_new->setVal(0);
    mu_g_new->copy(*m_leveldata[lev]->mu_g, 0, 0, m_leveldata[lev]->mu_g->nComp(),
                   m_leveldata[lev]->mu_g->nGrow(), m_leveldata[lev]->mu_g->nGrow());
    std::swap(m_leveldata[lev]->mu_g, mu_g_new);
    delete mu_g_new;

    // Gas velocity
    MultiFab* vel_g_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->vel_g->nComp(),
                                       m_leveldata[lev]->vel_g->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
    vel_g_new->setVal(0);
    vel_g_new->copy(*m_leveldata[lev]->vel_g, 0, 0, m_leveldata[lev]->vel_g->nComp(),
                    m_leveldata[lev]->vel_g->nGrow(), m_leveldata[lev]->vel_g->nGrow());
    std::swap(m_leveldata[lev]->vel_g, vel_g_new);
    delete vel_g_new;

    // Old gas velocity
    MultiFab* vel_go_new = new MultiFab(grids[lev], dmap[lev],
                                        m_leveldata[lev]->vel_go->nComp(),
                                        m_leveldata[lev]->vel_go->nGrow(),
                                        MFInfo(), *ebfactory[lev]);
    vel_go_new->setVal(0);
    vel_go_new->copy(*m_leveldata[lev]->vel_go, 0, 0, m_leveldata[lev]->vel_go->nComp(),
                     m_leveldata[lev]->vel_go->nGrow(), m_leveldata[lev]->vel_go->nGrow());
    std::swap(m_leveldata[lev]->vel_go, vel_go_new);
    delete vel_go_new;

    // Pressure gradients
    MultiFab* gp_new = new MultiFab(grids[lev], dmap[lev],
                                    m_leveldata[lev]->gp->nComp(),
                                    m_leveldata[lev]->gp->nGrow(),
                                    MFInfo(), *ebfactory[lev]);
    gp_new->setVal(0);
    gp_new->copy(*m_leveldata[lev]->gp, 0, 0, m_leveldata[lev]->gp->nComp(),
                 m_leveldata[lev]->gp->nGrow(), m_leveldata[lev]->gp->nGrow());
    std::swap(m_leveldata[lev]->gp, gp_new);
    delete gp_new;

    // Vorticity
    MultiFab* vort_new = new MultiFab(grids[lev], dmap[lev],
                                      m_leveldata[lev]->vort->nComp(),
                                      m_leveldata[lev]->vort->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    vort_new->setVal(0);
    vort_new->copy(*m_leveldata[lev]->vort, 0, 0, m_leveldata[lev]->vort->nComp(),
                   m_leveldata[lev]->vort->nGrow(), m_leveldata[lev]->vort->nGrow());
    std::swap(m_leveldata[lev]->vort, vort_new);
    delete vort_new;

    // Particle/fluid drag -- note it is important to copy from previous step in order to use in dt calculation
    MultiFab* drag_new = new MultiFab(grids[lev], dmap[lev],
                                      m_leveldata[lev]->drag->nComp(),
                                      m_leveldata[lev]->drag->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    drag_new->setVal(0.);
    drag_new->copy(*m_leveldata[lev]->drag, 0, 0, m_leveldata[lev]->drag->nComp(),
                   m_leveldata[lev]->drag->nGrow(), m_leveldata[lev]->drag->nGrow());
    std::swap(m_leveldata[lev]->drag, drag_new);
    delete drag_new;

    // Array to store the rhs for tensor diffusion solve
    MultiFab* diff_rhs_new = new MultiFab(grids[lev], dmap[lev], diff_rhs[lev]->nComp(),
                                          diff_rhs[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    diff_rhs_new->setVal(0);
    std::swap(diff_rhs[lev], diff_rhs_new);
    delete diff_rhs_new;

    // Array to store the rhs for tensor diffusion solve
    MultiFab* diff_rhs1_new = new MultiFab(grids[lev], dmap[lev], diff_rhs1[lev]->nComp(),
                                          diff_rhs1[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    diff_rhs1_new->setVal(0);
    std::swap(diff_rhs1[lev], diff_rhs1_new);
    delete diff_rhs1_new;

    // // Array to store the rhs for tensor diffusion solve
    MultiFab* diff_rhs4_new = new MultiFab(grids[lev], dmap[lev], diff_rhs4[lev]->nComp(),
                                          diff_rhs4[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    diff_rhs4_new->setVal(0);
    std::swap(diff_rhs4[lev], diff_rhs4_new);
    delete diff_rhs4_new;

    // Arrays to store the solution for diffusion solves
    MultiFab* diff_phi_new = new MultiFab(grids[lev], dmap[lev], diff_phi[lev]->nComp(),
                                          diff_phi[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    diff_phi_new->setVal(0);
    std::swap(diff_phi[lev], diff_phi_new);
    delete diff_phi_new;

    // Arrays to store the solution for diffusion solves
    MultiFab* diff_phi1_new = new MultiFab(grids[lev], dmap[lev], diff_phi1[lev]->nComp(),
                                          diff_phi1[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    diff_phi1_new->setVal(0);
    std::swap(diff_phi1[lev], diff_phi1_new);
    delete diff_phi1_new;

    // Arrays to store the solution for diffusion solves
    MultiFab* diff_phi4_new = new MultiFab(grids[lev], dmap[lev], diff_phi4[lev]->nComp(),
                                          diff_phi4[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    diff_phi4_new->setVal(0);
    std::swap(diff_phi4[lev], diff_phi4_new);
    delete diff_phi4_new;

    // Array to store the rhs for cell-centered solves
    MultiFab* mac_rhs_new = new MultiFab(grids[lev], dmap[lev], mac_rhs[lev]->nComp(),
                                          mac_rhs[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    mac_rhs_new->setVal(0);
    std::swap(mac_rhs[lev], mac_rhs_new);
    delete mac_rhs_new;

    // Arrays to store the solution for the MAC projection
    MultiFab* mac_phi_new = new MultiFab(grids[lev], dmap[lev], mac_phi[lev]->nComp(),
                                          mac_phi[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    mac_phi_new->setVal(0);
    std::swap(mac_phi[lev], mac_phi_new);
    delete mac_phi_new;

    // Slopes in x-direction
    MultiFab* xslopes_u_new = new MultiFab(grids[lev], dmap[lev],
                                           m_leveldata[lev]->xslopes_u->nComp(),
                                           m_leveldata[lev]->xslopes_u->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    xslopes_u_new->setVal(0);
    std::swap(m_leveldata[lev]->xslopes_u, xslopes_u_new);
    delete xslopes_u_new;

    MultiFab* xslopes_s_new = new MultiFab(grids[lev], dmap[lev],
                                           m_leveldata[lev]->xslopes_s->nComp(),
                                           m_leveldata[lev]->xslopes_s->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    xslopes_s_new->setVal(0);
    std::swap(m_leveldata[lev]->xslopes_s, xslopes_s_new);
    delete xslopes_s_new;

    // Slopes in y-direction
    MultiFab* yslopes_u_new = new MultiFab(grids[lev], dmap[lev],
                                           m_leveldata[lev]->yslopes_u->nComp(),
                                           m_leveldata[lev]->yslopes_u->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    yslopes_u_new->setVal(0);
    std::swap(m_leveldata[lev]->yslopes_u, yslopes_u_new);
    delete yslopes_u_new;

    MultiFab* yslopes_s_new = new MultiFab(grids[lev], dmap[lev],
                                           m_leveldata[lev]->yslopes_s->nComp(),
                                           m_leveldata[lev]->yslopes_s->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    yslopes_s_new->setVal(0.);
    std::swap(m_leveldata[lev]->yslopes_s, yslopes_s_new);
    delete yslopes_s_new;

    // Slopes in z-direction
    MultiFab* zslopes_u_new = new MultiFab(grids[lev], dmap[lev],
                                           m_leveldata[lev]->zslopes_u->nComp(),
                                           m_leveldata[lev]->zslopes_u->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    zslopes_u_new->setVal(0);
    std::swap(m_leveldata[lev]->zslopes_u, zslopes_u_new);
    delete zslopes_u_new;

    MultiFab* zslopes_s_new = new MultiFab(grids[lev], dmap[lev],
                                           m_leveldata[lev]->zslopes_s->nComp(),
                                           m_leveldata[lev]->zslopes_s->nGrow(),
                                           MFInfo(), *ebfactory[lev]);
    zslopes_s_new->setVal(0);
    std::swap(m_leveldata[lev]->zslopes_s, zslopes_s_new);
    delete zslopes_s_new;

   /****************************************************************************
    * x-face-based arrays                                                        *
    ****************************************************************************/

    BoxArray x_ba = grids[lev];
    x_ba = x_ba.surroundingNodes(0);

    // MAC/diffusion coefficient on x-faces
    MultiFab* bcx_mac_new = new MultiFab(x_ba, dmap[lev], bcoeff[lev][0]->nComp(),
                                         bcoeff[lev][0]->nGrow(), MFInfo(), *ebfactory[lev]);
    bcx_mac_new->setVal(0);
    std::swap(bcoeff[lev][0], bcx_mac_new);
    delete bcx_mac_new;

    // x component of MAC velocity
    MultiFab* u_mac_new = new MultiFab(x_ba, dmap[lev], u_mac[lev]->nComp(),
                                       u_mac[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    u_mac_new->setVal(0);
    std::swap(u_mac[lev], u_mac_new);
    delete u_mac_new;

   /****************************************************************************
    * y-face-based arrays                                                        *
    ****************************************************************************/

    BoxArray y_ba = grids[lev];
    y_ba = y_ba.surroundingNodes(1);

    // MAC/diffusion coefficient on y-faces
    MultiFab* bcy_mac_new = new MultiFab(y_ba, dmap[lev], bcoeff[lev][1]->nComp(),
                                         bcoeff[lev][1]->nGrow(), MFInfo(), *ebfactory[lev]);
    bcy_mac_new->setVal(0);
    std::swap(bcoeff[lev][1], bcy_mac_new);
    delete bcy_mac_new;

    // y component of MAC velocity
    MultiFab* v_mac_new = new MultiFab(y_ba, dmap[lev], v_mac[lev]->nComp(),
                                       v_mac[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    v_mac_new->setVal(0);
    std::swap(v_mac[lev], v_mac_new);
    delete v_mac_new;

   /****************************************************************************
    * z-face-based arrays                                                        *
    ****************************************************************************/

    BoxArray z_ba = grids[lev];
    z_ba = z_ba.surroundingNodes(2);

    // MAC/diffusion coefficient on z-faces
    MultiFab* bcz_mac_new = new MultiFab(z_ba, dmap[lev], bcoeff[lev][2]->nComp(),
                                         bcoeff[lev][2]->nGrow(), MFInfo(), *ebfactory[lev]);
    bcz_mac_new->setVal(0);
    std::swap(bcoeff[lev][2], bcz_mac_new);
    delete bcz_mac_new;

    // z component of MAC velocity
    MultiFab* w_mac_new = new MultiFab(z_ba, dmap[lev], w_mac[lev]->nComp(),
                                       w_mac[lev]->nGrow(), MFInfo(), *ebfactory[lev]);
    w_mac_new->setVal(0);
    std::swap(w_mac[lev], w_mac_new);
    delete w_mac_new;
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

   if ( particle_ebfactory[a_lev] == nullptr )
   {
      amrex::Print() << "Updating particle ebfactory 1" << std::endl;

      particle_ebfactory[a_lev] =
        new EBFArrayBoxFactory(*particle_eb_levels[a_lev], geom[a_lev], ba, dm,
                               {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                m_eb_full_grow_cells}, m_eb_support_level);

      changed = true;

   }
   else
   {
      amrex::Print() << "Updating particle ebfactory 2" << std::endl;

      const DistributionMapping&  eb_dm = particle_ebfactory[a_lev]->DistributionMap();
      const BoxArray&             eb_ba = particle_ebfactory[a_lev]->boxArray();

      if ( (dm != eb_dm) || (ba != eb_ba) )
      {
          particle_ebfactory[a_lev] =
              new EBFArrayBoxFactory(*particle_eb_levels[a_lev], geom[a_lev], ba, dm,
                                     {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                      m_eb_full_grow_cells}, m_eb_support_level);

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

   if ( ebfactory[a_lev] == nullptr )
   {
      Print() << "Updating ebfactory from nullptr" << std::endl;

      ebfactory[a_lev] =
          new EBFArrayBoxFactory(*eb_levels[a_lev], geom[a_lev], ba, dm,
                                 {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                  m_eb_full_grow_cells}, m_eb_support_level);

      is_updated = true;
   }
   else
   {
      const DistributionMapping&  eb_dm = ebfactory[a_lev]->DistributionMap();
      const BoxArray&             eb_ba = ebfactory[a_lev]->boxArray();

      if ( (dm != eb_dm) || (ba != eb_ba) )
      {
          Print() << "Updating ebfactory from existing" << std::endl;

          ebfactory[a_lev] =
              new EBFArrayBoxFactory(*eb_levels[a_lev], geom[a_lev], ba, dm,
                                     {m_eb_basic_grow_cells, m_eb_volume_grow_cells,
                                      m_eb_full_grow_cells}, m_eb_support_level);

          is_updated = true;
      }
   }

   return is_updated;
}
