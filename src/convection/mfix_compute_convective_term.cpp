#include <mfix.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_EB_utils.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>

//
// Compute the three components of the convection term
//
void
mfix::mfix_compute_convective_term (const bool update_laplacians,
                                    Vector< MultiFab* >& conv_u_in,
                                    Vector< MultiFab* >& conv_s_in,
                                    Vector< MultiFab* >& conv_X_in,
                                    Vector< MultiFab* >& lap_T_out,
                                    Vector< MultiFab* >& lap_X_out,
                                    Vector< MultiFab* > const& vel_in,
                                    Vector< MultiFab* > const& ep_g_in,
                                    Vector< MultiFab* > const& ep_u_mac,
                                    Vector< MultiFab* > const& ep_v_mac,
                                    Vector< MultiFab* > const& ep_w_mac,
                                    Vector< MultiFab* > const& ro_g_in,
                                    Vector< MultiFab* > const& MW_g_in,
                                    Vector< MultiFab* > const& T_g_in,
                                    Vector< MultiFab* > const& cp_g_in,
                                    Vector< MultiFab* > const& k_g_in,
                                    Vector< MultiFab* > const& h_g_in,
                                    Vector< MultiFab* > const& T_g_on_eb_in,
                                    Vector< MultiFab* > const& k_g_on_eb_in,
                                    Vector< MultiFab* > const& trac_in,
                                    Vector< MultiFab* > const& X_gk_in,
                                    Vector< MultiFab* > const& D_gk_in,
                                    Vector< MultiFab* > const& h_gk_in,
                                    Vector< MultiFab* > const& txfr_in,
                                    Vector< MultiFab* > const& ro_gk_txfr_in,
                                    Real time)
{
    BL_PROFILE("mfix::mfix_compute_convective_term");


    int slopes_comp; int conv_comp; int state_comp; int num_comp;

    // First do FillPatch of {velocity, density, tracer, enthalpy} so we know
    // the ghost cells of these arrays are all filled
    for (int lev = 0; lev < nlev; lev++)
    {
        // State with ghost cells
        MultiFab Sborder_u(grids[lev], dmap[lev], vel_in[lev]->nComp(), nghost,
                           MFInfo(), *ebfactory[lev]);
        FillPatchVel(lev, time, Sborder_u, 0, Sborder_u.nComp(), bcs_u);

        // Copy each FAB back from Sborder_u into the vel array, complete with filled ghost cells
        MultiFab::Copy(*vel_in[lev], Sborder_u, 0, 0, vel_in[lev]->nComp(), vel_in[lev]->nGrow());

        MultiFab Sborder_s(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]);

        // We FillPatch density even if not advecting it because we need it in the projections
        state_comp =  0; // comp = 0 --> density
        num_comp = 1;
        FillPatchScalar(lev, time, Sborder_s, state_comp, num_comp, bcs_s);
        MultiFab::Copy(*ro_g_in[lev], Sborder_s, 0, 0, num_comp, ro_g_in[lev]->nGrow());

        if (advect_tracer)
        {
           state_comp =  1; // comp = 1 --> tracer
           num_comp = 1;
           FillPatchScalar(lev, time, Sborder_s, state_comp, num_comp, bcs_s);
           MultiFab::Copy(*trac_in[lev], Sborder_s, 0, 0, num_comp, trac_in[lev]->nGrow());
        }

        if (advect_enthalpy)
        {
           state_comp =  5; // comp = 1 --> enthalpy
           num_comp = 1;
           FillPatchScalar(lev, time, Sborder_s, state_comp, num_comp, bcs_s);
           MultiFab::Copy(*h_g_in[lev], Sborder_s, 0, 0, num_comp, h_g_in[lev]->nGrow());
        }

        if (advect_fluid_species)
        {
          MultiFab Sborder_X(grids[lev], dmap[lev], FLUID::nspecies, nghost,
              MFInfo(), *ebfactory[lev]);

          Sborder_X.setVal(0);

          state_comp = 0;
          num_comp = FLUID::nspecies;

          FillPatchSpecies(lev, time, Sborder_X, state_comp, num_comp, bcs_X);

          MultiFab::Copy(*X_gk_in[lev], Sborder_X, 0, 0, num_comp,
              X_gk_in[lev]->nGrow());
        }

    }

    // Do projection on all AMR levels in one shot -- note that the {u_mac, v_mac, w_mac}
    //    arrays returned from this call are in fact {ep * u_mac, ep * v_mac, ep * w_mac}
    //    on face CENTROIDS
    compute_MAC_projected_velocities(time, vel_in, ep_u_mac, ep_v_mac, ep_w_mac,
       ep_g_in, ro_g_in, MW_g_in, T_g_in, cp_g_in, k_g_in, T_g_on_eb_in,
       k_g_on_eb_in, X_gk_in, D_gk_in, h_gk_in, txfr_in, ro_gk_txfr_in,
       update_laplacians, lap_T_out, lap_X_out);


    bool already_on_centroids = true;

    // Temporaries to store fluxes
    Vector< MultiFab* > fx(nlev);
    Vector< MultiFab* > fy(nlev);
    Vector< MultiFab* > fz(nlev);

    // Temporaries to store species fluxes
    Vector< MultiFab* > fx_X(nlev);
    Vector< MultiFab* > fy_X(nlev);
    Vector< MultiFab* > fz_X(nlev);

    Array<MultiFab*,3> fluxes;
    Array<MultiFab*,3> fluxes_X = {nullptr, nullptr, nullptr};

    for (int lev=0; lev < nlev; ++lev)
    {


        // We make these with ncomp = 3 so they can hold all three velocity components at once;
        //    note we can also use them to just hold the single density or tracer comp or enthalpy
        fx[lev] = new MultiFab(m_leveldata[lev]->u_mac->boxArray(), dmap[lev], 3, 2,
                               MFInfo(), *ebfactory[lev]);
        fy[lev] = new MultiFab(m_leveldata[lev]->v_mac->boxArray(), dmap[lev], 3, 2,
                               MFInfo(), *ebfactory[lev]);
        fz[lev] = new MultiFab(m_leveldata[lev]->w_mac->boxArray(), dmap[lev], 3, 2,
                               MFInfo(), *ebfactory[lev]);

        fx[lev]->setVal(covered_val);
        fy[lev]->setVal(covered_val);
        fz[lev]->setVal(covered_val);

        fluxes[0] = fx[lev];
        fluxes[1] = fy[lev];
        fluxes[2] = fz[lev];

        if (advect_fluid_species) {
          // We make these with ncomp = FLUID::nspecies so they can hold all
          // fluid species at once;
          fx_X[lev] = new MultiFab(m_leveldata[lev]->u_mac->boxArray(), dmap[lev],
              FLUID::nspecies, 2, MFInfo(), *ebfactory[lev]);
          fy_X[lev] = new MultiFab(m_leveldata[lev]->v_mac->boxArray(), dmap[lev],
              FLUID::nspecies, 2, MFInfo(), *ebfactory[lev]);
          fz_X[lev] = new MultiFab(m_leveldata[lev]->w_mac->boxArray(), dmap[lev],
              FLUID::nspecies, 2, MFInfo(), *ebfactory[lev]);

          fx_X[lev]->setVal(covered_val);
          fy_X[lev]->setVal(covered_val);
          fz_X[lev]->setVal(covered_val);

          fluxes_X[0] = fx_X[lev];
          fluxes_X[1] = fy_X[lev];
          fluxes_X[2] = fz_X[lev];
        }

        // We make this with ncomp = 3 so it can hold all three velocity
        // components at once; note we can also use it to just hold the single
        // density, enthalpy or tracer comp.
        // We note that it needs two ghost cells for the redistribution step.
        MultiFab conv_tmp(grids[lev], dmap[lev], 3, 2, MFInfo(), *ebfactory[lev]);
        conv_tmp.setVal(0.);

        // This MultiFab will be used for fluid species
        MultiFab* conv_X_tmp(nullptr);

        if (advect_fluid_species) {
          // We make this with ncomp = FLUID::nspecies so it can hold all the
          // fluid species at once.  We note that it needs two ghost cells for
          // the redistribution step.
          conv_X_tmp = new MultiFab(grids[lev], dmap[lev], FLUID::nspecies, 2,
              MFInfo(), *ebfactory[lev]);

          conv_X_tmp->setVal(0.);
        }

        // Compute slopes of velocity, density, enthalpy and tracer

        slopes_comp = 0;
        mfix_compute_slopes(lev, *vel_in[lev],
                            get_xslopes_u(), get_yslopes_u(), get_zslopes_u(),
                            slopes_comp, m_vel_g_bc_types);

        if (advect_density)
        {
           slopes_comp = 0;
           mfix_compute_slopes(lev, *ro_g_in[lev],
                               get_xslopes_s(), get_yslopes_s(), get_zslopes_s(),
                               slopes_comp, m_ro_g_bc_types);
        }

        if (advect_enthalpy)
        {
            // Convert enthalpy h_g to (rho * h_g) so we can use conservative update
            MultiFab::Multiply(*h_g_in[lev], *ro_g_in[lev], 0, 0, 1, h_g_in[lev]->nGrow());

            slopes_comp = 1;
            mfix_compute_slopes(lev, *h_g_in[lev],
                                get_xslopes_s(), get_yslopes_s(), get_zslopes_s(),
                                slopes_comp, m_T_g_bc_types);
        }

        if (advect_tracer)
        {
            // Convert tracer to (rho * tracer) so we can use conservative update
            MultiFab::Multiply(*trac_in[lev], *ro_g_in[lev], 0, 0, 1, trac_in[lev]->nGrow());

            slopes_comp = 2;
            mfix_compute_slopes(lev, *trac_in[lev],
                                get_xslopes_s(), get_yslopes_s(), get_zslopes_s(),
                                slopes_comp, m_trac_g_bc_types);
        }

        // Compute slopes of fluid species mass fractions
        if (advect_fluid_species)
        {
          // Convert mass fraction X_gk to (rho * X_gk) so we can use conservative update
          for (int n(0); n < FLUID::nspecies; n++)
            MultiFab::Multiply(*X_gk_in[lev], *ro_g_in[lev], 0, n, 1,
                X_gk_in[lev]->nGrow());

          slopes_comp = 0;

          mfix_compute_slopes(lev, *X_gk_in[lev], get_xslopes_X_gk(),
              get_yslopes_X_gk(), get_zslopes_X_gk(), slopes_comp,
              m_X_gk_bc_types);
        }

        // Initialize conv_s to 0 for both density, enthlapy and tracer
        conv_s_in[lev]->setVal(0, 0, conv_s_in[lev]->nComp(), conv_s_in[lev]->nGrow());

        if (advect_fluid_species)
          conv_X_in[lev]->setVal(0, 0, conv_X_in[lev]->nComp(),
              conv_X_in[lev]->nGrow());




        const EBFArrayBoxFactory* ebfact = &EBFactory(lev);

        const GpuArray<int, 3> bc_types =
          {bc_list.get_minf(), bc_list.get_pinf(), bc_list.get_pout()};

        // Array4<int> const& bct_ilo = bc_ilo[lev]->array();
        // Array4<int> const& bct_ihi = bc_ihi[lev]->array();
        // Array4<int> const& bct_jlo = bc_jlo[lev]->array();
        // Array4<int> const& bct_jhi = bc_jhi[lev]->array();
        // Array4<int> const& bct_klo = bc_klo[lev]->array();
        // Array4<int> const& bct_khi = bc_khi[lev]->array();



        // **************************************************
        // Compute div (ep_g u u) -- the update for velocity
        // **************************************************
        conv_comp = 0; state_comp = 0; num_comp = 3; slopes_comp = 0;

        mfix_compute_fluxes(lev, fx, fy, fz, vel_in, state_comp, num_comp,
                            get_xslopes_u(), get_yslopes_u(), get_zslopes_u(),
                            slopes_comp, ep_u_mac, ep_v_mac, ep_w_mac,
                            nghost, covered_val, bc_types,
                            bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                            bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                            bc_klo[lev]->array(), bc_khi[lev]->array(),
                            ebfact, geom);

        EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);

        single_level_weighted_redistribute(conv_tmp, *conv_u_in[lev],
               *ep_g_in[lev], conv_comp, num_comp, geom[lev]);

        // **************************************************
        // Compute div (ep_g rho u) -- the update for density
        // **************************************************
        if (advect_density)
        {
            conv_comp = 0; state_comp = 0; num_comp = 1; slopes_comp = 0;
            mfix_compute_fluxes(lev, fx, fy, fz, ro_g_in, state_comp, num_comp,
                                get_xslopes_s(), get_yslopes_s(), get_zslopes_s(),
                                slopes_comp, ep_u_mac, ep_v_mac, ep_w_mac,
                                nghost, covered_val, bc_types,
                                bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                bc_klo[lev]->array(), bc_khi[lev]->array(),
                                ebfact, geom);

            EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);

            single_level_weighted_redistribute(conv_tmp, *conv_s_in[lev],
                   *ep_g_in[lev], conv_comp, num_comp, geom[lev]);
        }

        // **********************************************************
        // Compute div (ep_g rho h_g u) -- the update for (rho*enthlapy)
        // **********************************************************
        if (advect_enthalpy)
        {
            conv_comp = 1; state_comp = 0; num_comp = 1; slopes_comp = 1;
            mfix_compute_fluxes(lev, fx, fy, fz, h_g_in, state_comp, num_comp,
                                get_xslopes_s(), get_yslopes_s(), get_zslopes_s(),
                                slopes_comp, ep_u_mac, ep_v_mac, ep_w_mac,
                                nghost, covered_val, bc_types,
                                bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                bc_klo[lev]->array(), bc_khi[lev]->array(),
                                ebfact, geom);

            EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);

            single_level_weighted_redistribute(conv_tmp, *conv_s_in[lev],
                   *ep_g_in[lev], conv_comp, num_comp, geom[lev]);

            // Convert (rho * enthalpy) back to enthalpy
            MultiFab::Divide(*h_g_in[lev],*ro_g_in[lev],0,0,1,h_g_in[lev]->nGrow());
        }

        // **********************************************************
        // Compute div (ep_g rho trac u) -- the update for (rho*trac)
        // **********************************************************
        if (advect_tracer)
        {
            conv_comp = 2; state_comp = 0; num_comp = 1; slopes_comp = 2;
            mfix_compute_fluxes(lev, fx, fy, fz, trac_in, state_comp, num_comp,
                                get_xslopes_s(), get_yslopes_s(), get_zslopes_s(),
                                slopes_comp, ep_u_mac, ep_v_mac, ep_w_mac,
                                nghost, covered_val, bc_types,
                                bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                                bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                                bc_klo[lev]->array(), bc_khi[lev]->array(),
                                ebfact, geom);


            EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);

            single_level_weighted_redistribute(conv_tmp, *conv_s_in[lev],
                   *ep_g_in[lev], conv_comp, num_comp, geom[lev]);

            // Convert (rho * tracer) back to tracer
            MultiFab::Divide(*trac_in[lev],*ro_g_in[lev],0,0,1,trac_in[lev]->nGrow());
        }

        // **************************************************
        // Compute div (ep_g rho X u) -- the update for species mass fraction
        // **************************************************
        if (advect_fluid_species)
        {
          conv_comp = 0;
          state_comp = 0;
          num_comp = FLUID::nspecies;
          slopes_comp = 0;

          mfix_compute_fluxes(lev, fx_X, fy_X, fz_X, X_gk_in, state_comp, num_comp,
              get_xslopes_X_gk(), get_yslopes_X_gk(), get_zslopes_X_gk(),
                              slopes_comp, ep_u_mac, ep_v_mac, ep_w_mac,
                              nghost, covered_val, bc_types,
                              bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                              bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                              bc_klo[lev]->array(), bc_khi[lev]->array(),
                              ebfact, geom);

          EB_computeDivergence(*conv_X_tmp, GetArrOfConstPtrs(fluxes_X), geom[lev], already_on_centroids);

          single_level_weighted_redistribute(*conv_X_tmp, *conv_X_in[lev],
              *ep_g_in[lev], conv_comp, num_comp, geom[lev]);

          // Convert (rho * mass_fractions) back to mass_fractions
          for (int n(0); n < FLUID::nspecies; n++)
            MultiFab::Divide(*X_gk_in[lev], *ro_g_in[lev], 0, n, 1,
                X_gk_in[lev]->nGrow());
        }

        // **********************************************************
        // Return the negative
        // **********************************************************
        conv_u_in[lev]->mult(-1.0);
        conv_s_in[lev]->mult(-1.0);

        if (advect_fluid_species) {
          conv_X_in[lev]->mult(-1.0);
          delete conv_X_tmp;
        }
    } // lev

    for (int lev(0); lev < nlev; ++lev) {
      delete fx[lev];
      delete fy[lev];
      delete fz[lev];

      if (advect_fluid_species) {
        delete fx_X[lev];
        delete fy_X[lev];
        delete fz_X[lev];
      }
    }
}
