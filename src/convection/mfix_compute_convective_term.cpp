#include <mfix.H>
#include <MOL.H>

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


    int conv_comp; int state_comp; int num_comp;

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

    Array<MultiFab*,3> fluxes;

    for (int lev=0; lev < nlev; ++lev)
    {

      const int ncomp = amrex::max(FLUID::nspecies, 3);

        // We make these with ncomp = 3 so they can hold all three velocity components at once;
        //    note we can also use them to just hold the single density or tracer comp or enthalpy
        fx[lev] = new MultiFab(m_leveldata[lev]->u_mac->boxArray(), dmap[lev], ncomp, 2,
                               MFInfo(), *ebfactory[lev]);
        fy[lev] = new MultiFab(m_leveldata[lev]->v_mac->boxArray(), dmap[lev], ncomp, 2,
                               MFInfo(), *ebfactory[lev]);
        fz[lev] = new MultiFab(m_leveldata[lev]->w_mac->boxArray(), dmap[lev], ncomp, 2,
                               MFInfo(), *ebfactory[lev]);

        fx[lev]->setVal(covered_val);
        fy[lev]->setVal(covered_val);
        fz[lev]->setVal(covered_val);

        fluxes[0] = fx[lev];
        fluxes[1] = fy[lev];
        fluxes[2] = fz[lev];

        // We make this with ncomp = 3 so it can hold all three velocity
        // components at once; note we can also use it to just hold the single
        // density, enthalpy or tracer comp.
        // We note that it needs two ghost cells for the redistribution step.
        MultiFab conv_tmp(grids[lev], dmap[lev], ncomp, 2, MFInfo(), *ebfactory[lev]);
        conv_tmp.setVal(0.);




        // Initialize conv_s to 0 for both density, enthlapy and tracer
        conv_s_in[lev]->setVal(0, 0, conv_s_in[lev]->nComp(), conv_s_in[lev]->nGrow());




        const EBFArrayBoxFactory* ebfact = &EBFactory(lev);

        const GpuArray<int, 2> bc_types =
          {bc_list.get_minf(), bc_list.get_pinf()};




        // **************************************************
        // Compute div (ep_u_mac * (u)) -- the update for velocity
        // **************************************************
        conv_comp = 0; state_comp = 0; num_comp = 3;

        mol::mfix_compute_fluxes(lev, fx, fy, fz, vel_in, state_comp, num_comp,
                            ep_u_mac, ep_v_mac, ep_w_mac,
                            nghost, covered_val, bc_types,
                            bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                            bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                            bc_klo[lev]->array(), bc_khi[lev]->array(),
                            ebfact, geom);

        EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);

        single_level_weighted_redistribute(conv_tmp, *conv_u_in[lev],
               *ep_g_in[lev], conv_comp, num_comp, geom[lev]);

        // **************************************************
        // Compute div (ep_u_mac * (rho)) -- the update for density
        // **************************************************
        if (advect_density)
        {
            conv_comp = 0; state_comp = 0; num_comp = 1;
            mol::mfix_compute_fluxes(lev, fx, fy, fz, ro_g_in, state_comp, num_comp,
                                ep_u_mac, ep_v_mac, ep_w_mac,
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
        // Compute div (ep_u_mac * (rho*h_g)) -- the update for (rho*enthlapy)
        // **********************************************************
        if (advect_enthalpy)
        {
          // Convert enthalpy h_g to (rho * h_g) so we can use conservative update
          MultiFab::Multiply(*h_g_in[lev], *ro_g_in[lev], 0, 0, 1, h_g_in[lev]->nGrow());

          conv_comp = 1; state_comp = 0; num_comp = 1;
            mol::mfix_compute_fluxes(lev, fx, fy, fz, h_g_in, state_comp, num_comp,
                                ep_u_mac, ep_v_mac, ep_w_mac,
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
        // Compute div (ep_u_mac * (rho*trac)) -- the update for (rho*trac)
        // **********************************************************
        if (advect_tracer)
        {
          // Convert tracer to (rho * tracer) so we can use conservative update
          MultiFab::Multiply(*trac_in[lev], *ro_g_in[lev], 0, 0, 1, trac_in[lev]->nGrow());

          conv_comp = 2; state_comp = 0; num_comp = 1;
          mol::mfix_compute_fluxes(lev, fx, fy, fz, trac_in, state_comp, num_comp,
                                ep_u_mac, ep_v_mac, ep_w_mac,
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
        // Compute div (ep_u_mac * (rho*X)) -- the update for species mass fraction
        // **************************************************
        if (advect_fluid_species)
        {

          // Convert mass fraction X_gk to (rho * X_gk) so we can use conservative update
          for (int n(0); n < FLUID::nspecies; n++)
            MultiFab::Multiply(*X_gk_in[lev], *ro_g_in[lev], 0, n, 1,
                X_gk_in[lev]->nGrow());

          conv_comp = 0;
          state_comp = 0;
          num_comp = FLUID::nspecies;

          mol::mfix_compute_fluxes(lev, fx, fy, fz, X_gk_in, state_comp, num_comp,
                              ep_u_mac, ep_v_mac, ep_w_mac,
                              nghost, covered_val, bc_types,
                              bc_ilo[lev]->array(), bc_ihi[lev]->array(),
                              bc_jlo[lev]->array(), bc_jhi[lev]->array(),
                              bc_klo[lev]->array(), bc_khi[lev]->array(),
                              ebfact, geom);

          EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);

          single_level_weighted_redistribute(conv_tmp, *conv_X_in[lev],
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
        }
    } // lev

    for (int lev(0); lev < nlev; ++lev) {
      delete fx[lev];
      delete fy[lev];
      delete fz[lev];
    }
}
