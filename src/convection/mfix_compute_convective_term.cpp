#include <mfix.H>
#include <param_mod_F.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>

//
// Compute the three components of the convection term
//
void
mfix::mfix_compute_convective_term (Vector< std::unique_ptr<MultiFab> >& conv_u_in,
                                    Vector< std::unique_ptr<MultiFab> >& conv_s_in,
                                    Vector< std::unique_ptr<MultiFab> >& vel_in,
                                    Vector< std::unique_ptr<MultiFab> >& ep_g_in,
                                    Vector< std::unique_ptr<MultiFab> >& ro_g_in,
                                    Vector< std::unique_ptr<MultiFab> >& trac_in,
                                    Real time)
{
    BL_PROFILE("mfix::mfix_compute_convective_term");

    // Temporaries to store fluxes
    Vector< std::unique_ptr<MultiFab> > fx;
    Vector< std::unique_ptr<MultiFab> > fy;
    Vector< std::unique_ptr<MultiFab> > fz;

    fx.resize(nlev);
    fy.resize(nlev);
    fz.resize(nlev);

    int slopes_comp; int conv_comp; int state_comp; int num_comp;

    // First do FillPatch of {velocity, density, tracer} so we know the ghost cells of
    // these arrays are all filled
    for (int lev = 0; lev <= finest_level; lev++)
    {
        // State with ghost cells
        MultiFab Sborder_u(grids[lev], dmap[lev], vel_in[lev]->nComp(), nghost,
                           MFInfo(), *ebfactory[lev]);
        FillPatchVel(lev, time, Sborder_u, 0, Sborder_u.nComp(), bcs_u);

        // Copy each FAB back from Sborder_u into the vel array, complete with filled ghost cells
        MultiFab::Copy(*vel_in[lev], Sborder_u, 0, 0, vel_in[lev]->nComp(), vel_in[lev]->nGrow());

        MultiFab Sborder_s(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]);

        // We FillPatch density even if not advecting it because we need it in the projections
        state_comp =  0; num_comp = 1;
        FillPatchScalar(lev, time, Sborder_s, state_comp, num_comp, bcs_s);
        MultiFab::Copy(*ro_g_in[lev], Sborder_s, 0, 0, num_comp, ro_g_in[lev]->nGrow());

        if (advect_tracer)
        {
           state_comp =  1; num_comp = 1;
           FillPatchScalar(lev, time, Sborder_s, state_comp, num_comp, bcs_s);
           MultiFab::Copy(*trac_in[lev], Sborder_s, 0, 0, num_comp, trac_in[lev]->nGrow());
        }

        // We make these with ncomp = 3 so they can hold all three velocity components at once;
        //    note we can also use them to just hold the single density or tracer comp
        fx[lev].reset(new MultiFab(u_mac[lev]->boxArray(),dmap[lev],3,2,MFInfo(),*ebfactory[lev]));
        fy[lev].reset(new MultiFab(v_mac[lev]->boxArray(),dmap[lev],3,2,MFInfo(),*ebfactory[lev]));
        fz[lev].reset(new MultiFab(w_mac[lev]->boxArray(),dmap[lev],3,2,MFInfo(),*ebfactory[lev]));

        // We need this to avoid FPE
        u_mac[lev]->setVal(covered_val);
        v_mac[lev]->setVal(covered_val);
        w_mac[lev]->setVal(covered_val);

        fx[lev]->setVal(covered_val);
        fy[lev]->setVal(covered_val);
        fz[lev]->setVal(covered_val);

        // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
        //    arrays returned from this call are in fact {ep * u_mac, ep * v_mac, ep * w_mac}
        //    on face CENTROIDS
        mfix_predict_vels_on_faces(lev, time, vel_in, u_mac, v_mac, w_mac, ep_g);
    }

    // Do projection on all AMR levels in one shot -- note that the {u_mac, v_mac, w_mac}
    //    arrays returned from this call are in fact {ep * u_mac, ep * v_mac, ep * w_mac}
    //    on face CENTROIDS
    apply_MAC_projection(u_mac, v_mac, w_mac, ep_g_in, ro_g_in, time);

    bool already_on_centroids = true;

    Array<MultiFab*,AMREX_SPACEDIM> fluxes;

    for (int lev=0; lev < nlev; ++lev)
    {

        fluxes[0] = fx[lev].get();
        fluxes[1] = fy[lev].get();
        fluxes[2] = fz[lev].get();

        // We make this with ncomp = 3 so it can hold all three velocity components at once;
        //    note we can also use it to just hold the single density or tracer comp
        // We note that it needs two ghost cells for the redistribution step.
        MultiFab conv_tmp(grids[lev], dmap[lev], 3, 2, MFInfo(), *ebfactory[lev]);
        conv_tmp.setVal(0.);

        if (advect_tracer)
        {
            // Convert tracer to (rho * tracer) so we can use conservative update
            MultiFab::Multiply(*trac_in[lev],*ro_g_in[lev],0,0,1,trac_in[lev]->nGrow());
        }

        // Compute slopes of density and tracer
        if (advect_density)
        {
           slopes_comp = 0;
           mfix_compute_slopes(lev, time, *ro_g_in[lev], xslopes_s, yslopes_s, zslopes_s, slopes_comp);
        }

        if (advect_tracer)
        {
           slopes_comp = 1;
           mfix_compute_slopes(lev, time, *trac_in[lev], xslopes_s, yslopes_s, zslopes_s, slopes_comp);
        }

        // Initialize conv_s to 0 for both density and tracer
        conv_s_in[lev]->setVal(0.,0,conv_s_in[lev]->nComp(),conv_s_in[lev]->nGrow());

        // **************************************************
        // Compute div (ep_g u u) -- the update for velocity
        // **************************************************
        conv_comp = 0; state_comp = 0; num_comp = 3; slopes_comp = 0;
        mfix_compute_fluxes(lev, fx, fy, fz, vel_in, state_comp, num_comp,
                            xslopes_u, yslopes_u, zslopes_u, slopes_comp,
                            u_mac, v_mac, w_mac);

        EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);
        single_level_weighted_redistribute(lev, conv_tmp, *conv_u_in[lev], *ep_g[lev], conv_comp, num_comp, geom);

        // **************************************************
        // Compute div (ep_g rho u) -- the update for density
        // **************************************************
        if (advect_density)
        {
            conv_comp = 0; state_comp = 0; num_comp = 1; slopes_comp = 0;
            mfix_compute_fluxes(lev, fx, fy, fz, ro_g_in, state_comp, num_comp,
                                xslopes_s, yslopes_s, zslopes_s, slopes_comp,
                                u_mac, v_mac, w_mac);
            EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);
            single_level_weighted_redistribute(lev, conv_tmp, *conv_s_in[lev], *ep_g[lev], conv_comp, num_comp, geom);
        }

        // **********************************************************
        // Compute div (ep_g rho trac u) -- the update for (rho*trac)
        // **********************************************************
        if (advect_tracer)
        {
            conv_comp = 1; state_comp = 0; num_comp = 1; slopes_comp = 1;
            mfix_compute_fluxes(lev, fx, fy, fz, trac_in, state_comp, num_comp,
                                xslopes_s, yslopes_s, zslopes_s, slopes_comp,
                                u_mac, v_mac, w_mac);
            EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);
            single_level_weighted_redistribute(lev, conv_tmp, *conv_s_in[lev], *ep_g[lev], conv_comp, num_comp, geom);
        }

        if (advect_tracer)
        {
           // Convert (rho * tracer) back to tracer
           MultiFab::Divide(*trac_in[lev],*ro_g_in[lev],0,0,1,trac_in[lev]->nGrow());
        }

        // Return the negative
        conv_u_in[lev]->mult(-1.0);
        conv_s_in[lev]->mult(-1.0);
    } // lev
}
