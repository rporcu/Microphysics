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

    // MAC velocity
    Vector< std::unique_ptr<MultiFab> > u_mac;
    Vector< std::unique_ptr<MultiFab> > v_mac;
    Vector< std::unique_ptr<MultiFab> > w_mac;

    // Temporaries to store fluxes 
    Vector< std::unique_ptr<MultiFab> > fx;
    Vector< std::unique_ptr<MultiFab> > fy;
    Vector< std::unique_ptr<MultiFab> > fz;

    u_mac.resize(nlev);
    v_mac.resize(nlev);
    w_mac.resize(nlev);

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

        BoxArray x_edge_ba = grids[lev];
        x_edge_ba.surroundingNodes(0);
 
        BoxArray y_edge_ba = grids[lev];
        y_edge_ba.surroundingNodes(1);
 
        BoxArray z_edge_ba = grids[lev];
        z_edge_ba.surroundingNodes(2);

        u_mac[lev].reset(new MultiFab(x_edge_ba,dmap[lev],1,2,MFInfo(),*ebfactory[lev]));
        v_mac[lev].reset(new MultiFab(y_edge_ba,dmap[lev],1,2,MFInfo(),*ebfactory[lev]));
        w_mac[lev].reset(new MultiFab(z_edge_ba,dmap[lev],1,2,MFInfo(),*ebfactory[lev]));
 
        // We make these with ncomp = 3 so they can hold all three velocity components at once;
        //    note we can also use them to just hold the single density or tracer comp
        fx[lev].reset(new MultiFab(x_edge_ba,dmap[lev],3,2,MFInfo(),*ebfactory[lev]));
        fy[lev].reset(new MultiFab(y_edge_ba,dmap[lev],3,2,MFInfo(),*ebfactory[lev]));
        fz[lev].reset(new MultiFab(z_edge_ba,dmap[lev],3,2,MFInfo(),*ebfactory[lev]));
 
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
    apply_MAC_projection (u_mac, v_mac, w_mac, ep_g_in, ro_g_in, time);

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

        //  USE THIS ONCE MULTI-COMPONENT WORKS AGAIN
        // conv_comp = 0; state_comp = 0; num_comp = 3; slopes_comp = 0;
        // mfix_compute_fluxes(lev, fx, fy, fz, conv_u_in, conv_comp+i, vel_in, state_comp+ num_comp,
        //                       xslopes_u, yslopes_u, zslopes_u, slopes_comp,
        //                       u_mac, v_mac, w_mac, false);
        // EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);
        // mfix_redistribute(lev, conv_tmp, conv_u_in, conv_comp, num_comp);

     
        conv_comp = 0; state_comp = 0; num_comp = 3; slopes_comp = 0;
        mfix_compute_fluxes(lev, fx, fy, fz, vel_in, state_comp, num_comp,
                            xslopes_u, yslopes_u, zslopes_u, slopes_comp,
                            u_mac, v_mac, w_mac);
        EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);
        mfix_redistribute(lev, conv_tmp, conv_u_in, conv_comp, num_comp);

        if (advect_density)
        {
            conv_comp = 0; state_comp = 0; num_comp = 1; slopes_comp = 0;
            mfix_compute_fluxes(lev, fx, fy, fz, ro_g_in, state_comp, num_comp,
                                xslopes_s, yslopes_s, zslopes_s, slopes_comp,
                                u_mac, v_mac, w_mac);
            EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);
            mfix_redistribute(lev, conv_tmp, conv_s_in, conv_comp, num_comp);
        }

        if (advect_tracer)
        {
            conv_comp = 1; state_comp = 0; num_comp = 1; slopes_comp = 1;
            mfix_compute_fluxes(lev, fx, fy, fz, trac_in, state_comp, num_comp,
                                xslopes_s, yslopes_s, zslopes_s, slopes_comp,
                                u_mac, v_mac, w_mac);
            EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);
            mfix_redistribute(lev, conv_tmp, conv_s_in, conv_comp, num_comp);
        }

        if (advect_tracer)
        {
           // Convert (rho * tracer) back to tracer
           MultiFab::Divide(*trac_in[lev],*ro_g_in[lev],0,0,1,trac_in[lev]->nGrow());
        }

        // Return the negative
        conv_u_in[lev] -> mult(-1.0);
        conv_s_in[lev] -> mult(-1.0);

    } // lev
#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif
}

void
mfix::mfix_redistribute( int lev, 
                         MultiFab& conv_tmp_in, 
                         Vector< std::unique_ptr<MultiFab> >& conv_out,
                         int conv_comp, int ncomp) 
{
    Box domain(geom[lev].Domain());

    EB_set_covered(conv_tmp_in, covered_val);
    conv_tmp_in.FillBoundary(geom[lev].periodicity());

    // Get EB geometric info
    Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
    const amrex::MultiFab*                    volfrac;

    areafrac  =   ebfactory[lev] -> getAreaFrac();
    volfrac   = &(ebfactory[lev] -> getVolFrac());

    for (MFIter mfi(*conv_out[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
          // Tilebox
          Box bx = mfi.tilebox ();

          // this is to check efficiently if this tile contains any eb stuff
          const EBFArrayBox&  conv_fab = static_cast<EBFArrayBox const&>((*conv_out[lev])[mfi]);
          const EBCellFlagFab&  flags = conv_fab.getEBCellFlagFab();

          if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
          {
                 // If tile is completely covered by EB geometry, set slopes
             // value to some very large number so we know if
             // we accidentally use these covered slopes later in calculations
             (*conv_out[lev])[mfi].setVal( get_my_huge(), bx, conv_comp, ncomp);
          }
          else
          {
             // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
             if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
             {
                (*conv_out[lev])[mfi].copy(conv_tmp_in[mfi],conv_comp,0,ncomp);
             }
             else
             {
                const int cyclic_x = geom[0].isPeriodic(0) ? 1 : 0;
                const int cyclic_y = geom[0].isPeriodic(1) ? 1 : 0;
                const int cyclic_z = geom[0].isPeriodic(2) ? 1 : 0;

                // Compute div(tau) with EB algorithm
                mfix_apply_eb_redistribution(bx, *conv_out[lev], conv_tmp_in, *ep_g[lev], &mfi,
                                             conv_comp, ncomp, flags, volfrac, domain,
                                             cyclic_x, cyclic_y, cyclic_z,
                                             geom[lev].CellSize());

             }
          }
    }
}

