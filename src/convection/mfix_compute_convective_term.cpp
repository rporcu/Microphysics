#include <mfix.H>
#include <mfix_divop_conv.hpp>
#include <param_mod_F.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>

//
// Compute the three components of the convection term
//
void
mfix::mfix_compute_convective_term( Vector< std::unique_ptr<MultiFab> >& conv_u_in, 
                                    Vector< std::unique_ptr<MultiFab> >& conv_s_in,
                                    Vector< std::unique_ptr<MultiFab> >& vel_in,
                                    Vector< std::unique_ptr<MultiFab> >& ep_g_in,
                                    Vector< std::unique_ptr<MultiFab> >& ro_g_in,
                                    Vector< std::unique_ptr<MultiFab> >& trac_in,
                                    Real time)
{
    BL_PROFILE("mfix::mfix_compute_convective_term");

    // First do FillPatch of {velocity, density, tracer} so we know the ghost cells of
    // these arrays are all filled
    for (int lev = 0; lev <= finest_level; lev++)
    {
        // State with ghost cells
        MultiFab Sborder_u(grids[lev], dmap[lev], vel_in[lev]->nComp(), nghost,
                           MFInfo(), *ebfactory[lev]);
        FillPatchVel(lev, time, Sborder_u, 0, Sborder_u.nComp(), bcs_u);

        // Copy each FAB back from Sborder_u into the vel array, complete with filled ghost cells
        MultiFab::Copy (*vel_in[lev], Sborder_u, 0, 0, vel_in[lev]->nComp(), vel_in[lev]->nGrow());

        if (advect_density || advect_tracer)
        {
            MultiFab Sborder_s(grids[lev], dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]);

            if (advect_density)
            {
               int icomp = 0;
               FillPatchScalar(lev, time, Sborder_s, icomp, bcs_s);
               MultiFab::Copy (*ro_g_in[lev], Sborder_s, 0, 0, 1, ro_g_in[lev]->nGrow());
            }

            if (advect_tracer)
            {
               int icomp = 1;
               FillPatchScalar(lev, time, Sborder_s, icomp, bcs_s);
               MultiFab::Copy (*trac_in[lev], Sborder_s, 0, 0, 1, trac_in[lev]->nGrow());
            }
        }
    }

    // MAC velocity
    Vector< std::unique_ptr<MultiFab> > u_mac;
    Vector< std::unique_ptr<MultiFab> > v_mac;
    Vector< std::unique_ptr<MultiFab> > w_mac;

    // Temporaries to store fluxes and convective term
    Array< std::unique_ptr<MultiFab>,AMREX_SPACEDIM> fx;
    Array< std::unique_ptr<MultiFab>,AMREX_SPACEDIM> fy;
    Array< std::unique_ptr<MultiFab>,AMREX_SPACEDIM> fz;

    u_mac.resize(nlev);
    v_mac.resize(nlev);
    w_mac.resize(nlev);

    for (int lev=0; lev < nlev; ++lev)
    {
       BoxArray x_edge_ba = grids[lev];
       x_edge_ba.surroundingNodes(0);
       u_mac[lev].reset(new MultiFab(x_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));
       u_mac[lev]->setVal(covered_val);

       BoxArray y_edge_ba = grids[lev];
       y_edge_ba.surroundingNodes(1);
       v_mac[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));
       v_mac[lev]->setVal(covered_val);

       BoxArray z_edge_ba = grids[lev];
       z_edge_ba.surroundingNodes(2);
       w_mac[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));
       w_mac[lev]->setVal(covered_val);

       for (int n = 0; n < 3; ++n )
       {
           fx[n].reset(new MultiFab(x_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));
           fy[n].reset(new MultiFab(y_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));
           fz[n].reset(new MultiFab(z_edge_ba,dmap[lev],1,nghost,MFInfo(),*ebfactory[lev]));

           // We need this to avoid FPE
           fx[n]->setVal(covered_val);
           fy[n]->setVal(covered_val);
           fz[n]->setVal(covered_val);
       }
    }

    // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
    //    arrays returned from this call are in fact {ep * u_mac, ep * v_mac, ep * w_mac}
    //    on face CENTROIDS
    mfix_predict_vels_on_faces(time, vel_in, u_mac, v_mac, w_mac, ep_g);

    // Do projection on all AMR levels in one shot -- note that the {u_mac, v_mac, w_mac}
    //    arrays returned from this call are in fact {ep * u_mac, ep * v_mac, ep * w_mac}
    //    on face CENTROIDS
    apply_MAC_projection (u_mac, v_mac, w_mac, ep_g_in, ro_g_in, time, steady_state );

    int slopes_comp; int conv_comp; int state_comp; int num_comp;

    for (int lev=0; lev < nlev; ++lev)
    {
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

        conv_comp = 0; state_comp = 0; num_comp = 3; slopes_comp = 0;
        mfix_compute_fluxes(lev, conv_u_in, conv_comp, vel_in, state_comp, num_comp,
                            xslopes_u, yslopes_u, zslopes_u, slopes_comp,
                            u_mac, v_mac, w_mac, false);

        conv_comp = 0; state_comp = 0; num_comp = 1; slopes_comp = 0;
        if (advect_density)
            mfix_compute_fluxes(lev, conv_s_in, conv_comp, ro_g_in, state_comp, num_comp,
                                xslopes_s, yslopes_s, zslopes_s, slopes_comp,
                                u_mac, v_mac, w_mac, false);

        conv_comp = 1; state_comp = 0; num_comp = 1; slopes_comp = 1;
        if (advect_tracer)
            mfix_compute_fluxes(lev, conv_s_in, conv_comp, trac_in, state_comp, num_comp,
                                xslopes_s, yslopes_s, zslopes_s, slopes_comp,
                                u_mac, v_mac, w_mac, false);

        if (advect_tracer)
        {
           // Convert (rho * tracer) back to tracer
           MultiFab::Divide(*trac_in[lev],*ro_g_in[lev],0,0,1,trac_in[lev]->nGrow());
        }
    } // lev
}

