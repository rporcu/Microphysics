#include <mfix.H>

#include <AMReX_VisMF.H>
#include <mfix_mf_helpers.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_calc_fluid_coeffs.H>


void
mfix::compute_laps (const bool update_lapT,
                    const bool update_lapTrac,
                    const bool update_lapX,
                    amrex::Vector< amrex::MultiFab*      >& lapT,
                    amrex::Vector< amrex::MultiFab*      >& lapTrac,
                    amrex::Vector< amrex::MultiFab*      >& lapX,
                    amrex::Vector< amrex::MultiFab*      > const& T_g,
                    amrex::Vector< amrex::MultiFab*      > const& trac,
                    amrex::Vector< amrex::MultiFab*      > const& X_gk,
                    amrex::Vector< amrex::MultiFab const*> const& ep_g,
                    amrex::Vector< amrex::MultiFab const*> const& ro_g)

{
  if (update_lapT) {
    diffusion_op->ComputeLapT(lapT, T_g, ep_g, get_k_g_const(),
                              get_T_g_on_eb_const(), get_k_g_on_eb_const());
  }

  if (update_lapTrac) {
    diffusion_op->ComputeLapS(lapTrac, trac, ro_g, ep_g, mfix::mu_s);
  }

  if (update_lapX) {
    diffusion_op->ComputeLapX(lapX, X_gk, ro_g, ep_g, get_D_gk_const());
  }

}
