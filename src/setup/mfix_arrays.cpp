#include <mfix.H>
#include <AMReX_EB_utils.H>

void
mfix::AllocateArrays (int lev)
{
    if (ooo_debug) amrex::Print() << "AllocateArrays" << std::endl;
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

    // Tracer in gas
    trac[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    trac_o[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    trac[lev]->setVal(0.);
    trac_o[lev]->setVal(0.);

    const BoxArray & nd_grids = amrex::convert(grids[lev], IntVect{1,1,1});

    p0_g[lev].reset(new MultiFab(nd_grids,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    p0_g[lev]->setVal(0.);

    p_g[lev].reset(new MultiFab(nd_grids,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    p_g[lev]->setVal(0.);

    p_go[lev].reset(new MultiFab(nd_grids,dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    p_go[lev]->setVal(0.);

    phi_nd[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0, MFInfo(), *ebfactory[lev]));
    phi_nd[lev]->setVal(0.);

    diveu[lev].reset(new MultiFab(nd_grids,dmap[lev],1,0, MFInfo(), *ebfactory[lev]));
    diveu[lev]->setVal(0.);

    // Presssure gradients
    gp[lev].reset(new MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    gp[lev]->setVal(0.);

    // Molecular viscosity
    mu_g[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    mu_g[lev]->setVal(0.);

    // Current velocity
    vel_g[lev].reset(new MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    vel_g[lev]->setVal(0.);

    // Old velocity
    vel_go[lev].reset(new  MultiFab(grids[lev],dmap[lev],3,nghost, MFInfo(), *ebfactory[lev]));
    vel_go[lev]->setVal(0.);

    // Vorticity
    vort[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    vort[lev]->setVal(0.);

    // This is the deposition of the drag force onto the grid
    // 0,1,2 is (drag coefficient * particle velocity)
    // 4 is drag coefficient
    drag[lev].reset(new  MultiFab(grids[lev],dmap[lev],4,nghost, MFInfo(), *ebfactory[lev]));
    drag[lev]->setVal(0.);

    // Array to store the rhs for diffusion solves -- no ghost cells needed
    diff_rhs[lev].reset(new MultiFab(grids[lev],dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
    diff_rhs[lev] -> setVal(0.);

    // Array to store the solution for diffusion solves -- only 1 ghost cell needed
    diff_phi[lev].reset(new MultiFab(grids[lev],dmap[lev], 3, 1, MFInfo(), *ebfactory[lev]));
    diff_phi[lev] -> setVal(0.);

    // Array to store the rhs for MAC projection
    mac_rhs[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    mac_rhs[lev] -> setVal(0.);

    // Array to store the solution for MAC projections
    mac_phi[lev].reset(new MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    mac_phi[lev] -> setVal(0.);

    // Slopes in x-direction
    xslopes_u[lev].reset(new  MultiFab(grids[lev],dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), *ebfactory[lev]));
    xslopes_u[lev] -> setVal(0.);
    xslopes_s[lev].reset(new  MultiFab(grids[lev],dmap[lev], 2, nghost, MFInfo(), *ebfactory[lev]));
    xslopes_s[lev] -> setVal(0.);

    // Slopes in y-direction
    yslopes_u[lev].reset(new  MultiFab(grids[lev],dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), *ebfactory[lev]));
    yslopes_u[lev] -> setVal(0.);
    yslopes_s[lev].reset(new  MultiFab(grids[lev],dmap[lev], 2, nghost, MFInfo(), *ebfactory[lev]));
    yslopes_s[lev] -> setVal(0.);

    // Slopes in z-direction
    zslopes_u[lev].reset(new  MultiFab(grids[lev],dmap[lev], AMREX_SPACEDIM, nghost, MFInfo(), *ebfactory[lev]));
    zslopes_u[lev] -> setVal(0.);
    zslopes_s[lev].reset(new  MultiFab(grids[lev],dmap[lev], 2, nghost, MFInfo(), *ebfactory[lev]));
    zslopes_s[lev] -> setVal(0.);

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

    // When the pressure is on nodes, bcoeff_nd is at cell centers
    bcoeff_nd[lev].reset(new  MultiFab(grids[lev],dmap[lev],1,nghost, MFInfo(), *ebfactory[lev]));
    bcoeff_nd[lev]->setVal(0.);

    // ****************************************************************

    // Create a BoxArray on x-faces.
    BoxArray x_edge_ba = grids[lev];
    x_edge_ba.surroundingNodes(0);
    bcoeff[lev][0].reset(new MultiFab(x_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));

    // Create a BoxArray on y-faces.
    BoxArray y_edge_ba = grids[lev];
    y_edge_ba.surroundingNodes(1);
    bcoeff[lev][1].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));

    // Create a BoxArray on z-faces.
    BoxArray z_edge_ba = grids[lev];
    z_edge_ba.surroundingNodes(2);
    bcoeff[lev][2].reset(new MultiFab(z_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));

    bcoeff[lev][0]->setVal(0.);
    bcoeff[lev][1]->setVal(0.);
    bcoeff[lev][2]->setVal(0.);
}


void
mfix::RegridArrays (int lev)
{
    if (ooo_debug) amrex::Print() << "RegridArrays" << std::endl;
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

    // Tracer in gas
    ng = trac[lev]->nGrow();
    std::unique_ptr<MultiFab> trac_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    trac_new->setVal(0.0);
    trac_new->copy(*trac[lev],0,0,1,0,ng);
    trac[lev] = std::move(trac_new);

    // Old tracer in gas
    ng = trac_o[lev]->nGrow();
    std::unique_ptr<MultiFab> trac_o_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    trac_o_new->setVal(0.0);
    trac_o_new->copy(*trac_o[lev],0,0,1,0,ng);
    trac_o[lev] = std::move(trac_o_new);

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

    ng = diveu[lev]->nGrow();
    std::unique_ptr<MultiFab> diveu_new(new MultiFab(nd_grids,dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    diveu[lev] = std::move(diveu_new);
    diveu[lev]->setVal(0.);

    ng = phi_nd[lev] -> nGrow();
    std::unique_ptr<MultiFab> phi_new(new MultiFab(nd_grids,dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    phi_nd[lev] = std::move(phi_new);
    phi_nd[lev]->setVal(0.);

    std::unique_ptr<MultiFab> bc0_new(new MultiFab(grids[lev],dmap[lev],1,bcoeff_nd[lev]->nGrow(),
                                                   MFInfo(), *ebfactory[lev] ));
    bcoeff_nd[lev] = std::move(bc0_new);
    bcoeff_nd[lev]->setVal(0.);

    // Molecular viscosity
    ng = mu_g[lev]->nGrow();
    std::unique_ptr<MultiFab> mu_g_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    mu_g_new->setVal(0.0);
    mu_g_new->copy(*mu_g[lev],0,0,1,0,ng);
    mu_g[lev] = std::move(mu_g_new);

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
    gp_new->copy(*gp[lev],0,0,gp[lev]->nComp(),0,ng);
    gp[lev] = std::move(gp_new);

    // Vorticity
    ng = vort[lev]->nGrow();
    std::unique_ptr<MultiFab> vort_new(new MultiFab(grids[lev],dmap[lev],1,ng, MFInfo(), *ebfactory[lev]));
    vort[lev] = std::move(vort_new);
    vort[lev]->setVal(0.);

    // Particle/fluid drag
    ng = drag[lev]->nGrow();
    std::unique_ptr<MultiFab> drag_new(new MultiFab(grids[lev],dmap[lev],drag[lev]->nComp(),ng,
                                                    MFInfo(), *ebfactory[lev]));
    drag[lev] = std::move(drag_new);
    drag[lev]->setVal(0.);

    // Array to store the rhs for tensor diffusion solve
    std::unique_ptr<MultiFab> diff_rhs_new(new  MultiFab(grids[lev], dmap[lev], 3, diff_rhs[lev]->nComp(), 
                                                         MFInfo(), *ebfactory[lev]));
    diff_rhs[lev] = std::move(diff_rhs_new);
    diff_rhs[lev] -> setVal(0.);

    // Arrays to store the solution for diffusion solves
    std::unique_ptr<MultiFab> diff_phi_new(new  MultiFab(grids[lev], dmap[lev], 3, diff_phi[lev]->nComp(),
                                                         MFInfo(), *ebfactory[lev]));
    diff_phi[lev] = std::move(diff_phi_new);
    diff_phi[lev] -> setVal(0.);

    // Array to store the rhs for cell-centered solves
    std::unique_ptr<MultiFab> mac_rhs_new(new  MultiFab(grids[lev], dmap[lev], 1, nghost,
                                                        MFInfo(), *ebfactory[lev]));
    mac_rhs[lev] = std::move(mac_rhs_new);
    mac_rhs[lev] -> setVal(0.);

    // Arrays to store the solution for the MAC projection
    std::unique_ptr<MultiFab> mac_phi_new(new  MultiFab(grids[lev], dmap[lev], 1, nghost,
                                                        MFInfo(), *ebfactory[lev]));
    mac_phi[lev] = std::move(mac_phi_new);
    mac_phi[lev] -> setVal(0.);

    // Slopes in x-direction
    ng = xslopes_u[lev] -> nGrow();
    std::unique_ptr<MultiFab> xslopes_u_new(new  MultiFab(grids[lev], dmap[lev],xslopes_u[lev]->nComp(),nghost,
                                                          MFInfo(), *ebfactory[lev]));
    xslopes_u[lev] = std::move(xslopes_u_new);
    xslopes_u[lev] -> setVal(0.);

    ng = xslopes_s[lev] -> nGrow();
    std::unique_ptr<MultiFab> xslopes_s_new(new  MultiFab(grids[lev], dmap[lev],xslopes_s[lev]->nComp(),nghost,
                                                          MFInfo(), *ebfactory[lev]));
    xslopes_s[lev] = std::move(xslopes_s_new);

    // Slopes in y-direction
    ng = yslopes_u[lev] -> nGrow();
    std::unique_ptr<MultiFab> yslopes_u_new(new  MultiFab(grids[lev], dmap[lev],yslopes_u[lev]->nComp(),nghost,
                                                          MFInfo(), *ebfactory[lev]));
    yslopes_u[lev] = std::move(yslopes_u_new);
    yslopes_u[lev] -> setVal(0.);

    ng = yslopes_s[lev] -> nGrow();
    std::unique_ptr<MultiFab> yslopes_s_new(new  MultiFab(grids[lev], dmap[lev],yslopes_s[lev]->nComp(),nghost,
                                                          MFInfo(), *ebfactory[lev]));
    yslopes_s[lev] = std::move(yslopes_s_new);
    yslopes_s[lev] -> setVal(0.);

    // Slopes in z-direction
    ng = zslopes_u[lev] -> nGrow();
    std::unique_ptr<MultiFab> zslopes_u_new(new  MultiFab(grids[lev], dmap[lev],zslopes_u[lev]->nComp(),nghost,
                                                        MFInfo(), *ebfactory[lev]));
    zslopes_u[lev] = std::move(zslopes_u_new);
    zslopes_u[lev] -> setVal(0.);

    ng = zslopes_s[lev] -> nGrow();
    std::unique_ptr<MultiFab> zslopes_s_new(new  MultiFab(grids[lev], dmap[lev],zslopes_s[lev]->nComp(),nghost,
                                                        MFInfo(), *ebfactory[lev]));
    zslopes_s[lev] = std::move(zslopes_s_new);
    zslopes_s[lev] -> setVal(0.);

   /****************************************************************************
    * Face-based Arrays                                                        *
    ****************************************************************************/

    BoxArray x_ba = grids[lev];
    x_ba = x_ba.surroundingNodes(0);

    // MAC/diffusion coefficient on x-faces
    std::unique_ptr<MultiFab> bcx_mac_new(new MultiFab(x_ba,dmap[lev],1,nghost,MFInfo(), *ebfactory[lev]));
    bcoeff[lev][0] = std::move(bcx_mac_new);
    bcoeff[lev][0] -> setVal(0.0);

   //****************************************************************************

    BoxArray y_ba = grids[lev];
    y_ba = y_ba.surroundingNodes(1);

    // MAC/diffusion coefficient on y-faces
    std::unique_ptr<MultiFab> bcy_mac_new(new MultiFab(y_ba,dmap[lev],1,nghost,MFInfo(), *ebfactory[lev]));
    bcoeff[lev][1] = std::move(bcy_mac_new);
    bcoeff[lev][1] -> setVal(0.0);

   //****************************************************************************

    BoxArray z_ba = grids[lev];
    z_ba = z_ba.surroundingNodes(2);

    // MAC/diffusion coefficient on z-faces
    std::unique_ptr<MultiFab> bcz_mac_new(new MultiFab(z_ba,dmap[lev],1,nghost,MFInfo(), *ebfactory[lev]));
    bcoeff[lev][2] = std::move(bcz_mac_new);
    bcoeff[lev][2] -> setVal(0.0);

    // ********************************************************************************
    // Make sure we fill the ghost cells as appropriate -- this is copied from init_fluid
    // ********************************************************************************

      ep_g[lev]->FillBoundary(geom[lev].periodicity());
     ep_go[lev]->FillBoundary(geom[lev].periodicity());

      ro_g[lev]->FillBoundary(geom[lev].periodicity());
     ro_go[lev]->FillBoundary(geom[lev].periodicity());

        mu_g[lev]->FillBoundary(geom[lev].periodicity());
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

       std::unique_ptr<MultiFab> new_level_set(new MultiFab);

       if (level_sets[a_lev]->boxArray() == nd_ba)
       {       
           MFUtil::regrid(* new_level_set, nd_ba, dm, * level_sets[a_lev], true);
       }
       else
       {
           int nc = level_sets[a_lev]->nComp();
           int ng = level_sets[a_lev]->nGrow();
           const Periodicity& period = geom[a_lev].periodicity();
           new_level_set->define(nd_ba, dm, nc, ng);
           new_level_set->copy(*level_sets[a_lev], 0, 0, nc, 0, ng, period);
       }
       
       level_sets[a_lev] = std::move(new_level_set);

       //________________________________________________________________________
       // If we're operating in single-level mode, the level-set has a second
       // (refined) MultiFab that also needs to be regridded.
       if ((nlev == 1) && (a_lev == 0))
       {
           Print() << "Also regridding refined level-set" << std::endl;

           BoxArray ref_nd_ba = amrex::convert(ba, IntVect::TheNodeVector());
           ref_nd_ba.refine(levelset__refinement);

           std::unique_ptr<MultiFab> new_level_set(new MultiFab);

           if (level_sets[a_lev+1]->boxArray() == ref_nd_ba)
           {       
               MFUtil::regrid(* new_level_set, ref_nd_ba, dm, * level_sets[a_lev+1], true);
           }
           else
           {
               int nc = level_sets[a_lev+1]->nComp();
               int ng = level_sets[a_lev+1]->nGrow();
               const Periodicity& period = geom[a_lev].periodicity();
               new_level_set->define(ref_nd_ba, dm, nc, ng);
               new_level_set->copy(*level_sets[a_lev+1], 0, 0, nc, 0, ng, period);
           }
           
           level_sets[a_lev+1] = std::move(new_level_set);
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
