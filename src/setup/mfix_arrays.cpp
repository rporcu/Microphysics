#include <mfix.H>
#include <AMReX_EB_utils.H>
#include <mfix_eb_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>

void
mfix::ResizeArrays ()
{
    int nlevs_max = maxLevel() + 1;

    m_leveldata.resize(nlevs_max);
    for (int lev(0); lev < nlevs_max; ++lev)
      m_leveldata[lev].reset(new LevelData());

    bcoeff.resize(nlevs_max);

    // Fluid grid EB factory
    ebfactory.resize(nlevs_max);

    /****************************************************************************
     *                                                                          *
     * Initialize particle data (including level-set data)                      *
     * NOTE: the level-set data (as well as implicit functions) live on at      *
     *       least two levels                                                   *
     *                                                                          *
     ***************************************************************************/

    // Particle grid EB factory
    particle_ebfactory.resize(nlevs_max);

    eb_levels.resize(amrex::max(2, nlevs_max));
    particle_eb_levels.resize(amrex::max(2, nlevs_max));

    level_sets.resize(amrex::max(2, nlevs_max));
}

void
mfix::AllocateArrays (int lev)
{
    if (ooo_debug) amrex::Print() << "AllocateArrays" << std::endl;
    mfix_update_ebfactory(lev);

    // ********************************************************************************
    // Cell- or node-based arrays
    // ********************************************************************************

    m_leveldata[lev].reset(new LevelData(grids[lev], dmap[lev], nghost_state(),
                                         *ebfactory[lev]));
    m_leveldata[lev]->resetValues(covered_val);

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

    // ****************************************************************

    BoxArray ba = grids[lev];
    // Create a BoxArray on x-faces.
    // x-face-based coefficient for MAC and diffusive solves
    if (bcoeff[lev][0] != nullptr)
      delete bcoeff[lev][0];

    bcoeff[lev][0] = new MultiFab(BoxArray(ba).surroundingNodes(0), dmap[lev], 1,
                                  nghost_state(), MFInfo(), *ebfactory[lev]);
    bcoeff[lev][0]->setVal(0.);

    // Create a BoxArray on y-faces.
    // y-face-based coefficient for MAC and diffusive solves
    if (bcoeff[lev][1] != nullptr)
      delete bcoeff[lev][1];

    bcoeff[lev][1] = new MultiFab(BoxArray(ba).surroundingNodes(1), dmap[lev], 1,
                                  nghost_state(), MFInfo(), *ebfactory[lev]);
    bcoeff[lev][1]->setVal(0.);

    // Create a BoxArray on z-faces.
    // z-face-based coefficient for MAC and diffusive solves
    if (bcoeff[lev][2] != nullptr)
      delete bcoeff[lev][2];

    bcoeff[lev][2] = new MultiFab(BoxArray(ba).surroundingNodes(2), dmap[lev], 1,
                                  nghost_state(), MFInfo(), *ebfactory[lev]);
    bcoeff[lev][2]->setVal(0.);
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
    //       However, we want to be sure to only use valid regions of the src, so we use
    //       src_ngrow = 0 and dst_ngrow = all
    //
    //
    int src_ngrow = 0;

    // Void fraction
    MultiFab& ep_g = *(m_leveldata[lev]->ep_g);
    MultiFab* ep_g_new = new MultiFab(grids[lev], dmap[lev], ep_g.nComp(),
                                      ep_g.nGrow(), MFInfo(), *ebfactory[lev]);
    ep_g_new->setVal(1.);
    ep_g_new->ParallelCopy(ep_g, 0, 0, ep_g.nComp(), src_ngrow,  ep_g.nGrow(), geom[lev].periodicity());
    std::swap((m_leveldata[lev]->ep_g), ep_g_new);
    delete ep_g_new;

    // Gas density
    MultiFab* ro_g_new = new MultiFab(grids[lev], dmap[lev],
                                      m_leveldata[lev]->ro_g->nComp(),
                                      m_leveldata[lev]->ro_g->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    ro_g_new->setVal(0);
    ro_g_new->ParallelCopy(*m_leveldata[lev]->ro_g, 0, 0, m_leveldata[lev]->ro_g->nComp(),
                   src_ngrow, m_leveldata[lev]->ro_g->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->ro_g, ro_g_new);
    delete ro_g_new;

    // Old gas density
    MultiFab* ro_go_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->ro_go->nComp(),
                                       m_leveldata[lev]->ro_go->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
    ro_go_new->setVal(0);
    ro_go_new->ParallelCopy(*m_leveldata[lev]->ro_go, 0, 0, m_leveldata[lev]->ro_go->nComp(),
                    src_ngrow, m_leveldata[lev]->ro_go->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->ro_go, ro_go_new);
    delete ro_go_new;

    // Gas molecular weight
    MultiFab* MW_g_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->MW_g->nComp(),
                                       m_leveldata[lev]->MW_g->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
    MW_g_new->setVal(0);
    MW_g_new->ParallelCopy(*m_leveldata[lev]->MW_g, 0, 0, m_leveldata[lev]->MW_g->nComp(),
                    src_ngrow, m_leveldata[lev]->MW_g->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->MW_g, MW_g_new);
    delete MW_g_new;

    // Tracer in gas
    MultiFab* trac_new = new MultiFab(grids[lev], dmap[lev],
                                      m_leveldata[lev]->trac->nComp(),
                                      m_leveldata[lev]->trac->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    trac_new->setVal(0);
    trac_new->ParallelCopy(*m_leveldata[lev]->trac, 0, 0, m_leveldata[lev]->trac->nComp(),
                   src_ngrow, m_leveldata[lev]->trac->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->trac, trac_new);
    delete trac_new;

    // Old tracer in gas
    MultiFab* trac_o_new = new MultiFab(grids[lev], dmap[lev],
                                        m_leveldata[lev]->trac_o->nComp(),
                                        m_leveldata[lev]->trac_o->nGrow(),
                                        MFInfo(), *ebfactory[lev]);
    trac_o_new->setVal(0);
    trac_o_new->ParallelCopy(*m_leveldata[lev]->trac_o, 0, 0, m_leveldata[lev]->trac_o->nComp(),
                     src_ngrow, m_leveldata[lev]->trac_o->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->trac_o, trac_o_new);
    delete trac_o_new;

    const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

    // Gas perturbational pressure
    MultiFab* p_g_new = new MultiFab(nd_grids, dmap[lev],
                                     m_leveldata[lev]->p_g->nComp(),
                                     m_leveldata[lev]->p_g->nGrow(),
                                     MFInfo(), *ebfactory[lev]);
    p_g_new->setVal(0);
    p_g_new->ParallelCopy(*m_leveldata[lev]->p_g, 0, 0, m_leveldata[lev]->p_g->nComp(),
                  src_ngrow, m_leveldata[lev]->p_g->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->p_g, p_g_new);
    delete p_g_new;

    // Old gas perturbational pressure
    MultiFab* p_go_new = new MultiFab(nd_grids, dmap[lev],
                                      m_leveldata[lev]->p_go->nComp(),
                                      m_leveldata[lev]->p_go->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    p_go_new->setVal(0);
    p_go_new->ParallelCopy(*m_leveldata[lev]->p_go, 0, 0, m_leveldata[lev]->p_go->nComp(),
                   src_ngrow, m_leveldata[lev]->p_go->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->p_go, p_go_new);
    delete p_go_new;

    // 
    MultiFab* p0_g_new = new MultiFab(nd_grids, dmap[lev],
                                      m_leveldata[lev]->p0_g->nComp(),
                                      m_leveldata[lev]->p0_g->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    p0_g_new->setVal(0);
    p0_g_new->ParallelCopy(*m_leveldata[lev]->p0_g, 0, 0, m_leveldata[lev]->p0_g->nComp(),
                   src_ngrow, m_leveldata[lev]->p0_g->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->p0_g, p0_g_new);
    delete p0_g_new;

    MultiFab* diveu_new = new MultiFab(nd_grids, dmap[lev], m_leveldata[lev]->diveu->nComp(),
                                       m_leveldata[lev]->diveu->nGrow(), MFInfo(), *ebfactory[lev]);
    diveu_new->setVal(0);
    diveu_new->ParallelCopy(*m_leveldata[lev]->diveu, 0, 0, m_leveldata[lev]->diveu->nComp(),
                    src_ngrow, m_leveldata[lev]->diveu->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->diveu, diveu_new);
    delete diveu_new;

    // Molecular viscosity
    MultiFab* mu_g_new = new MultiFab(grids[lev], dmap[lev],
                                      m_leveldata[lev]->mu_g->nComp(),
                                      m_leveldata[lev]->mu_g->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    mu_g_new->setVal(0);
    mu_g_new->ParallelCopy(*m_leveldata[lev]->mu_g, 0, 0, m_leveldata[lev]->mu_g->nComp(),
                   src_ngrow, m_leveldata[lev]->mu_g->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->mu_g, mu_g_new);
    delete mu_g_new;

    // rank of fluid grids
    MultiFab* ba_proc_new = new MultiFab(grids[lev], dmap[lev],
                                         m_leveldata[lev]->ba_proc->nComp(),
                                         m_leveldata[lev]->ba_proc->nGrow(),
                                         MFInfo(), *ebfactory[lev]);
    const Real proc = Real(ParallelDescriptor::MyProc());
    for (MFIter mfi(*ba_proc_new, false); mfi.isValid(); ++mfi)
    {
      amrex::Array4<Real> const& bx_proc = ba_proc_new->array(mfi);
      ParallelFor(mfi.validbox(), [bx_proc, proc] 
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                  { bx_proc(i,j,k) = proc; });
    }
    std::swap(m_leveldata[lev]->ba_proc, ba_proc_new);
    delete ba_proc_new;

    if (advect_enthalpy) {
      // Gas thermodynamic pressure
      MultiFab* pressure_g_new = new MultiFab(grids[lev], dmap[lev],
                                              m_leveldata[lev]->pressure_g->nComp(),
                                              m_leveldata[lev]->pressure_g->nGrow(),
                                              MFInfo(), *ebfactory[lev]);
      pressure_g_new->setVal(0);
      pressure_g_new->ParallelCopy(*m_leveldata[lev]->pressure_g, 0, 0,
                                   m_leveldata[lev]->pressure_g->nComp(),
                                   src_ngrow, m_leveldata[lev]->pressure_g->nGrow(),
                                   geom[lev].periodicity());
      std::swap(m_leveldata[lev]->pressure_g, pressure_g_new);
      delete pressure_g_new;

      // Old gas thermodynamic pressure
      MultiFab* pressure_go_new = new MultiFab(grids[lev], dmap[lev],
                                               m_leveldata[lev]->pressure_go->nComp(),
                                               m_leveldata[lev]->pressure_go->nGrow(),
                                               MFInfo(), *ebfactory[lev]);
      pressure_go_new->setVal(0);
      pressure_go_new->ParallelCopy(*m_leveldata[lev]->pressure_go, 0, 0,
                                    m_leveldata[lev]->pressure_go->nComp(),
                                    src_ngrow, m_leveldata[lev]->pressure_go->nGrow(),
                                    geom[lev].periodicity());
      std::swap(m_leveldata[lev]->pressure_go, pressure_go_new);
      delete pressure_go_new;

      // Gas temperature
      MultiFab* T_g_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->T_g->nComp(),
                                       m_leveldata[lev]->T_g->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
      T_g_new->setVal(0);
      T_g_new->ParallelCopy(*m_leveldata[lev]->T_g, 0, 0, m_leveldata[lev]->T_g->nComp(),
                     src_ngrow, m_leveldata[lev]->T_g->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->T_g, T_g_new);
      delete T_g_new;

      // Old gas temperature
      MultiFab* T_go_new = new MultiFab(grids[lev], dmap[lev],
                                        m_leveldata[lev]->T_go->nComp(),
                                        m_leveldata[lev]->T_go->nGrow(),
                                        MFInfo(), *ebfactory[lev]);
      T_go_new->setVal(0);
      T_go_new->ParallelCopy(*m_leveldata[lev]->T_go, 0, 0, m_leveldata[lev]->T_go->nComp(),
                      src_ngrow, m_leveldata[lev]->T_go->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->T_go, T_go_new);
      delete T_go_new;

      // Specific heat
      MultiFab* cp_g_new = new MultiFab(grids[lev], dmap[lev],
                                        m_leveldata[lev]->cp_g->nComp(),
                                        m_leveldata[lev]->cp_g->nGrow(),
                                        MFInfo(), *ebfactory[lev]);
      cp_g_new->setVal(0);
      cp_g_new->ParallelCopy(*m_leveldata[lev]->cp_g, 0, 0, m_leveldata[lev]->cp_g->nComp(),
                     src_ngrow, m_leveldata[lev]->cp_g->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->cp_g, cp_g_new);
      delete cp_g_new;

      // Thermal conductivity
      MultiFab* k_g_new = new MultiFab(grids[lev], dmap[lev],
                                        m_leveldata[lev]->k_g->nComp(),
                                        m_leveldata[lev]->k_g->nGrow(),
                                        MFInfo(), *ebfactory[lev]);
      k_g_new->setVal(0);
      k_g_new->ParallelCopy(*m_leveldata[lev]->k_g, 0, 0, m_leveldata[lev]->k_g->nComp(),
                     src_ngrow, m_leveldata[lev]->k_g->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->k_g, k_g_new);
      delete k_g_new;

      // Gas enthalpy
      MultiFab* h_g_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->h_g->nComp(),
                                       m_leveldata[lev]->h_g->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
      h_g_new->setVal(0);
      h_g_new->ParallelCopy(*m_leveldata[lev]->h_g, 0, 0, m_leveldata[lev]->h_g->nComp(),
                     src_ngrow, m_leveldata[lev]->h_g->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->h_g, h_g_new);
      delete h_g_new;

      // Old gas enthalpy
      MultiFab* h_go_new = new MultiFab(grids[lev], dmap[lev],
                                        m_leveldata[lev]->h_go->nComp(),
                                        m_leveldata[lev]->h_go->nGrow(),
                                        MFInfo(), *ebfactory[lev]);
      h_go_new->setVal(0);
      h_go_new->ParallelCopy(*m_leveldata[lev]->h_go, 0, 0, m_leveldata[lev]->h_go->nComp(),
                      src_ngrow, m_leveldata[lev]->h_go->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->h_go, h_go_new);
      delete h_go_new;

      if (EB::fix_temperature) {
        // Dirichlet temperature values on the EB
        MultiFab* T_g_on_eb_new = new MultiFab(grids[lev], dmap[lev],
            m_leveldata[lev]->T_g_on_eb->nComp(),
            m_leveldata[lev]->T_g_on_eb->nGrow(), MFInfo(), *ebfactory[lev]);

        T_g_on_eb_new->setVal(0);

        T_g_on_eb_new->ParallelCopy(*m_leveldata[lev]->T_g_on_eb, 0, 0,
            m_leveldata[lev]->T_g_on_eb->nComp(), src_ngrow,
            m_leveldata[lev]->T_g_on_eb->nGrow(), geom[lev].periodicity());

        std::swap(m_leveldata[lev]->T_g_on_eb, T_g_on_eb_new);
        delete T_g_on_eb_new;

        // Dirichlet thermal diffusivity values on the EB
        MultiFab* k_g_on_eb_new = new MultiFab(grids[lev], dmap[lev],
            m_leveldata[lev]->k_g_on_eb->nComp(),
            m_leveldata[lev]->k_g_on_eb->nGrow(), MFInfo(), *ebfactory[lev]);

        k_g_on_eb_new->setVal(0);

        k_g_on_eb_new->ParallelCopy(*m_leveldata[lev]->k_g_on_eb, 0, 0,
            m_leveldata[lev]->k_g_on_eb->nComp(), src_ngrow,
            m_leveldata[lev]->k_g_on_eb->nGrow(), geom[lev].periodicity());

        std::swap(m_leveldata[lev]->k_g_on_eb, k_g_on_eb_new);
        delete k_g_on_eb_new;
      }
    }

    if (advect_fluid_species) {
      // Gas species mass fraction
      MultiFab* X_gk_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->X_gk->nComp(),
                                       m_leveldata[lev]->X_gk->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
      X_gk_new->setVal(0);
      X_gk_new->ParallelCopy(*m_leveldata[lev]->X_gk, 0, 0,
          m_leveldata[lev]->X_gk->nComp(), src_ngrow,
          m_leveldata[lev]->X_gk->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->X_gk, X_gk_new);
      delete X_gk_new;

      // Old gas species mass fraction
      MultiFab* X_gko_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->X_gk->nComp(),
                                       m_leveldata[lev]->X_gk->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
      X_gko_new->setVal(0);
      X_gko_new->ParallelCopy(*m_leveldata[lev]->X_gko, 0, 0,
          m_leveldata[lev]->X_gko->nComp(), src_ngrow,
          m_leveldata[lev]->X_gko->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->X_gko, X_gko_new);
      delete X_gko_new;

      // Species diffusion coefficients
      MultiFab* D_gk_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->D_gk->nComp(),
                                       m_leveldata[lev]->D_gk->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
      D_gk_new->setVal(0);
      D_gk_new->ParallelCopy(*m_leveldata[lev]->D_gk, 0, 0, m_leveldata[lev]->D_gk->nComp(),
                     src_ngrow, m_leveldata[lev]->D_gk->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->D_gk, D_gk_new);
      delete D_gk_new;
    }

    if (advect_enthalpy and advect_fluid_species) {
      // Species specific heat
      MultiFab* cp_gk_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->cp_gk->nComp(),
                                       m_leveldata[lev]->cp_gk->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
      cp_gk_new->setVal(0);
      cp_gk_new->ParallelCopy(*m_leveldata[lev]->cp_gk, 0, 0, m_leveldata[lev]->cp_gk->nComp(),
                     src_ngrow, m_leveldata[lev]->cp_gk->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->cp_gk, cp_gk_new);
      delete cp_gk_new;

      // Species enthalpy
      MultiFab* h_gk_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->h_gk->nComp(),
                                       m_leveldata[lev]->h_gk->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
      h_gk_new->setVal(0);
      h_gk_new->ParallelCopy(*m_leveldata[lev]->h_gk, 0, 0, m_leveldata[lev]->h_gk->nComp(),
                     src_ngrow, m_leveldata[lev]->h_gk->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->h_gk, h_gk_new);
      delete h_gk_new;
    }

    // Gas velocity
    MultiFab* vel_g_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->vel_g->nComp(),
                                       m_leveldata[lev]->vel_g->nGrow(),
                                       MFInfo(), *ebfactory[lev]);
    vel_g_new->setVal(0);
    vel_g_new->ParallelCopy(*m_leveldata[lev]->vel_g, 0, 0, m_leveldata[lev]->vel_g->nComp(),
                    src_ngrow, m_leveldata[lev]->vel_g->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->vel_g, vel_g_new);
    delete vel_g_new;

    // Old gas velocity
    MultiFab* vel_go_new = new MultiFab(grids[lev], dmap[lev],
                                        m_leveldata[lev]->vel_go->nComp(),
                                        m_leveldata[lev]->vel_go->nGrow(),
                                        MFInfo(), *ebfactory[lev]);
    vel_go_new->setVal(0);
    vel_go_new->ParallelCopy(*m_leveldata[lev]->vel_go, 0, 0, m_leveldata[lev]->vel_go->nComp(),
                     src_ngrow, m_leveldata[lev]->vel_go->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->vel_go, vel_go_new);
    delete vel_go_new;

    // Pressure gradients
    MultiFab* gp_new = new MultiFab(grids[lev], dmap[lev],
                                    m_leveldata[lev]->gp->nComp(),
                                    m_leveldata[lev]->gp->nGrow(),
                                    MFInfo(), *ebfactory[lev]);
    gp_new->setVal(0);
    gp_new->ParallelCopy(*m_leveldata[lev]->gp, 0, 0, m_leveldata[lev]->gp->nComp(),
                 src_ngrow, m_leveldata[lev]->gp->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->gp, gp_new);
    delete gp_new;

    // Vorticity
    MultiFab* vort_new = new MultiFab(grids[lev], dmap[lev],
                                      m_leveldata[lev]->vort->nComp(),
                                      m_leveldata[lev]->vort->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    vort_new->setVal(0);
    vort_new->ParallelCopy(*m_leveldata[lev]->vort, 0, 0, m_leveldata[lev]->vort->nComp(),
                   src_ngrow, m_leveldata[lev]->vort->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->vort, vort_new);
    delete vort_new;

    // Particle/fluid interphase transfer -- note it is important to copy from previous step in order to use in dt calculation
    MultiFab* txfr_new = new MultiFab(grids[lev], dmap[lev],
                                      m_leveldata[lev]->txfr->nComp(),
                                      m_leveldata[lev]->txfr->nGrow(),
                                      MFInfo(), *ebfactory[lev]);
    txfr_new->setVal(0.);
    txfr_new->ParallelCopy(*m_leveldata[lev]->txfr, 0, 0, m_leveldata[lev]->txfr->nComp(),
                   src_ngrow, m_leveldata[lev]->txfr->nGrow(), geom[lev].periodicity());
    std::swap(m_leveldata[lev]->txfr, txfr_new);
    delete txfr_new;

    if (advect_fluid_species and solve_reactions) {
      // Species mass transfer rates
      MultiFab* chem_txfr_new = new MultiFab(grids[lev], dmap[lev],
                                        m_leveldata[lev]->chem_txfr->nComp(),
                                        m_leveldata[lev]->chem_txfr->nGrow(),
                                        MFInfo(), *ebfactory[lev]);
      chem_txfr_new->setVal(0.);
      chem_txfr_new->ParallelCopy(*m_leveldata[lev]->chem_txfr, 0, 0, m_leveldata[lev]->chem_txfr->nComp(),
                     src_ngrow, m_leveldata[lev]->chem_txfr->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->chem_txfr, chem_txfr_new);
      delete chem_txfr_new;
    }

    // Arrays to store the solution for the MAC projection
    MultiFab* mac_phi_new = new MultiFab(grids[lev], dmap[lev],
                                         m_leveldata[lev]->mac_phi->nComp(),
                                         m_leveldata[lev]->mac_phi->nGrow(),
                                         MFInfo(), *ebfactory[lev]);
    mac_phi_new->setVal(0);
    std::swap(m_leveldata[lev]->mac_phi, mac_phi_new);
    delete mac_phi_new;


    // Arrays to store the diffusive term div(tau)
    MultiFab* divtau_new = new MultiFab(grids[lev], dmap[lev],
                                         m_leveldata[lev]->divtau_o->nComp(),
                                         m_leveldata[lev]->divtau_o->nGrow(),
                                         MFInfo(), *ebfactory[lev]);
    divtau_new->setVal(0);
    std::swap(m_leveldata[lev]->divtau_o, divtau_new);
    delete divtau_new;


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

      particle_ebfactory[a_lev].reset(
        new EBFArrayBoxFactory(*particle_eb_levels[a_lev], geom[a_lev], ba, dm,
                               {nghost_eb_basic(), nghost_eb_volume(),
                                nghost_eb_full()}, m_eb_support_level));

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
              new EBFArrayBoxFactory(*particle_eb_levels[a_lev], geom[a_lev], ba, dm,
                                     {nghost_eb_basic(), nghost_eb_volume(),
                                      nghost_eb_full()}, m_eb_support_level));

         changed = true;
      }
   }

   if (changed)
   {

       Print() << "Regridding level-set on lev = " << a_lev << std::endl;

       const BoxArray nd_ba = amrex::convert(ba, IntVect::TheNodeVector());

       std::unique_ptr<MultiFab> new_level_set{new MultiFab()};

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

           new_level_set->copy(*level_sets[a_lev], 0, 0, nc, ng,
               ng, period);
       }

       std::swap(level_sets[a_lev], new_level_set);

       //________________________________________________________________________
       // If we're operating in single-level mode, the level-set has a second
       // (refined) MultiFab that also needs to be regridded.
       if ((nlev == 1) && (a_lev == 0))
       {
           Print() << "Also regridding refined level-set" << std::endl;

           BoxArray ref_nd_ba = amrex::convert(ba, IntVect::TheNodeVector());
           ref_nd_ba.refine(levelset_refinement);

           std::unique_ptr<MultiFab> new_level_set_lev {new MultiFab()};

           if (level_sets[a_lev+1]->boxArray() == ref_nd_ba)
           {
               MFUtil::regrid(*new_level_set_lev, ref_nd_ba, dm, *level_sets[a_lev+1], true);
           }
           else
           {
               int nc = level_sets[a_lev+1]->nComp();
               int ng = level_sets[a_lev+1]->nGrow();
               const Periodicity& period = geom[a_lev].periodicity();
               new_level_set_lev->define(ref_nd_ba, dm, nc, ng);
               new_level_set_lev->setVal(0.0);

               new_level_set_lev->copy(*level_sets[a_lev+1], 0, 0, nc,
                   ng, ng, period);
           }

           std::swap(level_sets[a_lev+1], new_level_set_lev);
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

      ebfactory[a_lev].reset(
          new EBFArrayBoxFactory(*eb_levels[a_lev], geom[a_lev], ba, dm,
                                 {nghost_eb_basic(), nghost_eb_volume(),
                                  nghost_eb_full()}, m_eb_support_level));

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
              new EBFArrayBoxFactory(*eb_levels[a_lev], geom[a_lev], ba, dm,
                                     {nghost_eb_basic(), nghost_eb_volume(),
                                      nghost_eb_full()}, m_eb_support_level));

          is_updated = true;
      }
   }

   return is_updated;
}
