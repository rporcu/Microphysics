#include <mfix.H>

#include <AMReX_VisMF.H>
#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>

using namespace amrex;

void
mfix::compute_laps (const bool update_lapT,
                    const bool update_lapTrac,
                    const bool update_flux,
                    Vector< MultiFab*      > const& lapT,
                    Vector< MultiFab*      > const& lapTrac,
                    Vector< Array< MultiFab*, AMREX_SPACEDIM> > const& J_gk,
                    Vector< MultiFab*      > const& T_g,
                    Vector< MultiFab*      > const& trac,
                    Vector< MultiFab*      > const& X_gk,
                    Vector< MultiFab const*> const& ep_g,
                    Vector< MultiFab const*> const& ro_g)
{
  if (update_lapT) {
    diffusion_op->ComputeLapT(lapT, T_g, ep_g, get_T_g_on_eb_const());
  }

  if (update_lapTrac) {
    diffusion_op->ComputeLapS(lapTrac, trac, ro_g, ep_g, mfix::mu_s);
  }

  if (update_flux) {
    diffusion_op->ComputeFlux(J_gk, X_gk, ro_g, ep_g);
  }

}
