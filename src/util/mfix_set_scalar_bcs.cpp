#include <mfix.H>

#include <eos_mod.H>

#include <param_mod_F.H>
#include <MFIX_FLUID_Parms.H>

using namespace amrex;

//
// Set the BCs for all the variables EXCEPT pressure or velocity.
//
void
mfix::mfix_set_scalar_bcs (Real time,
                           Vector< MultiFab* > const& cp_g_in,
                           Vector< MultiFab* > const& k_g_in,
                           Vector< MultiFab* > const& mu_g_in)
{
  BL_PROFILE("mfix::mfix_set_scalar_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*(m_leveldata[lev]->ep_g), false); mfi.isValid(); ++mfi)
     {
        set_thermal_conductivity_bcs(time, lev, (*k_g_in[lev])[mfi], domain);
        set_specific_heat_bcs(time, lev, (*cp_g_in[lev])[mfi], domain);
        set_viscosity_bcs(time, lev, (*mu_g_in[lev])[mfi], domain);
     }

     k_g_in[lev] -> FillBoundary (geom[lev].periodicity());
     cp_g_in[lev] -> FillBoundary (geom[lev].periodicity());
     mu_g_in[lev] -> FillBoundary (geom[lev].periodicity());

     EB_set_covered(*k_g_in[lev], 0, k_g_in[lev]->nComp(), k_g_in[lev]->nGrow(), covered_val);
     EB_set_covered(*cp_g_in[lev], 0, cp_g_in[lev]->nComp(), cp_g_in[lev]->nGrow(), covered_val);
     EB_set_covered(*mu_g_in[lev], 0, mu_g_in[lev]->nComp(), mu_g_in[lev]->nGrow(), covered_val);

  }
}
