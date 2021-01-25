#include <mfix.H>
#include <MOL.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_EB_utils.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>

using namespace amrex;
//
// Compute the three components of the convection term
//
void
mol::compute_convective_term (const int lev,
                              amrex::Vector<amrex::BoxArray> grids,
                              DistributionMapping const& dmap,
                              Vector< MultiFab* >& conv_u_in,
                              Vector< MultiFab* >& conv_s_in,
                              Vector< MultiFab* >& conv_X_in,
                              const bool advect_density,
                              const bool advect_enthalpy,
                              const bool advect_tracer,
                              const bool advect_species,
                              Vector< MultiFab* >& fx,
                              Vector< MultiFab* >& fy,
                              Vector< MultiFab* >& fz,
                              Vector< MultiFab* > const& vel_in,
                              Vector< MultiFab* > const& ep_g_in,
                              Vector< MultiFab* > const& ep_u_mac,
                              Vector< MultiFab* > const& ep_v_mac,
                              Vector< MultiFab* > const& ep_w_mac,
                              Vector< MultiFab* > const& ro_g_in,
                              Vector< MultiFab* > const& h_g_in,
                              Vector< MultiFab* > const& trac_in,
                              Vector< MultiFab* > const& X_gk_in,
                              const GpuArray<int, 3> bc_types,
                              std::map<std::string, Gpu::DeviceVector<int>>& velocity_bcs,
                              std::map<std::string, Gpu::DeviceVector<int>>& density_bcs,
                              std::map<std::string, Gpu::DeviceVector<int>>& enthalpy_bcs,
                              std::map<std::string, Gpu::DeviceVector<int>>& tracer_bcs,
                              std::map<std::string, Gpu::DeviceVector<int>>& species_bcs,
                              Array4<int const> const& bct_ilo,
                              Array4<int const> const& bct_ihi,
                              Array4<int const> const& bct_jlo,
                              Array4<int const> const& bct_jhi,
                              Array4<int const> const& bct_klo,
                              Array4<int const> const& bct_khi,
                              EBFArrayBoxFactory const* ebfact,
                              Vector<Geometry> geom)
{
    BL_PROFILE("mol::compute_advection_term");
    const Real covered_val = 1.e40;

    bool already_on_centroids = true;

    const int maxncomp = amrex::max(FLUID::nspecies, 3);

    BoxArray ba = grids[lev];

     Array<MultiFab*,3> fluxes;
     fluxes[0] = fx[lev];
     fluxes[1] = fy[lev];
     fluxes[2] = fz[lev];

    fx[lev]->setVal(covered_val);
    fy[lev]->setVal(covered_val);
    fz[lev]->setVal(covered_val);

    // We make this with ncomp = 3 so it can hold all three velocity
    // components at once; note we can also use it to just hold the single
    // density, enthalpy or tracer comp.
    // We note that it needs two ghost cells for the redistribution step.
    MultiFab conv_tmp(grids[lev], dmap, maxncomp, 2, MFInfo(), *ebfact);
    conv_tmp.setVal(0.);

    // Initialize conv_s to 0 for both density, enthlapy and tracer
    conv_s_in[lev]->setVal(0, 0, conv_s_in[lev]->nComp(), conv_s_in[lev]->nGrow());

    int conv_comp; int state_comp; int num_comp;

    // Four ghost cells are required when using EB
    const int nghost = 4;

    // **************************************************
    // Compute div (ep_u_mac * (u)) -- the update for velocity
    // **************************************************
    conv_comp = 0; state_comp = 0; num_comp = 3;

    mol::compute_convective_fluxes(lev, fx, fy, fz, vel_in, state_comp, num_comp,
            ep_u_mac, ep_v_mac, ep_w_mac, nghost, covered_val, bc_types, velocity_bcs,
            bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi, ebfact, geom);

    EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes),
                         geom[lev], already_on_centroids);

    single_level_weighted_redistribute(conv_tmp, *conv_u_in[lev],
           *ep_g_in[lev], conv_comp, num_comp, geom[lev]);


    // **************************************************
    // Compute div (ep_u_mac * (rho)) -- the update for density
    // **************************************************
    if (advect_density)
    {
        conv_comp = 0; state_comp = 0; num_comp = 1;
        mol::compute_convective_fluxes(lev, fx, fy, fz, ro_g_in, state_comp, num_comp,
                ep_u_mac, ep_v_mac, ep_w_mac, nghost, covered_val, bc_types, density_bcs,
                bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi, ebfact, geom);

        EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes),
                             geom[lev], already_on_centroids);

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
        mol::compute_convective_fluxes(lev, fx, fy, fz, h_g_in, state_comp, num_comp,
                ep_u_mac, ep_v_mac, ep_w_mac, nghost, covered_val, bc_types, enthalpy_bcs,
                bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi, ebfact, geom);

        EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes),
                             geom[lev], already_on_centroids);

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
      mol::compute_convective_fluxes(lev, fx, fy, fz, trac_in, state_comp, num_comp,
              ep_u_mac, ep_v_mac, ep_w_mac, nghost, covered_val, bc_types, tracer_bcs,
              bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi, ebfact, geom);


        EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes), geom[lev], already_on_centroids);

        single_level_weighted_redistribute(conv_tmp, *conv_s_in[lev],
               *ep_g_in[lev], conv_comp, num_comp, geom[lev]);

        // Convert (rho * tracer) back to tracer
        MultiFab::Divide(*trac_in[lev],*ro_g_in[lev],0,0,1,trac_in[lev]->nGrow());
    }

    // **************************************************
    // Compute div (ep_u_mac * (rho*X)) -- the update for species mass fraction
    // **************************************************
    if (advect_species)
    {

      // Convert mass fraction X_gk to (rho * X_gk) so we can use conservative update
      for (int n(0); n < FLUID::nspecies; n++)
        MultiFab::Multiply(*X_gk_in[lev], *ro_g_in[lev], 0, n, 1,
            X_gk_in[lev]->nGrow());

      conv_comp = 0;
      state_comp = 0;
      num_comp = FLUID::nspecies;

      mol::compute_convective_fluxes(lev, fx, fy, fz, X_gk_in, state_comp, num_comp,
              ep_u_mac, ep_v_mac, ep_w_mac, nghost, covered_val, bc_types, species_bcs,
              bct_ilo, bct_ihi, bct_jlo, bct_jhi, bct_klo, bct_khi, ebfact, geom);

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

    if (advect_species) {
      conv_X_in[lev]->mult(-1.0);
    }

}
