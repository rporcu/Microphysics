#include <mfix.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>

using namespace amrex;

//
// Set the BCs for all the variables EXCEPT pressure or velocity.
//
void
mfix::mfix_set_scalar_bcs (Real time,
                           Vector< MultiFab* > const& mu_g_in,
                           Vector< MultiFab* > const& cp_g_in,
                           Vector< MultiFab* > const& k_g_in,
                           Vector< MultiFab* > const& MW_g_in)
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
        set_viscosity_bcs(time, lev, (*mu_g_in[lev])[mfi], domain);
        set_molecular_weight_bcs(time, lev, (*MW_g_in[lev])[mfi], domain);
        
        if (advect_enthalpy) {
          set_specific_heat_bcs(time, lev, (*cp_g_in[lev])[mfi], domain);
          set_thermal_conductivity_bcs(time, lev, (*k_g_in[lev])[mfi], domain);
        }
     }

     mu_g_in[lev] -> FillBoundary (geom[lev].periodicity());
     MW_g_in[lev] -> FillBoundary (geom[lev].periodicity());

     if (advect_enthalpy) {
       cp_g_in[lev] -> FillBoundary (geom[lev].periodicity());
       k_g_in[lev] -> FillBoundary (geom[lev].periodicity());
     }

     EB_set_covered(*mu_g_in[lev], 0, mu_g_in[lev]->nComp(),
         mu_g_in[lev]->nGrow(), covered_val);
     EB_set_covered(*MW_g_in[lev], 0, MW_g_in[lev]->nComp(),
         MW_g_in[lev]->nGrow(), covered_val);

     if (advect_enthalpy) {
       EB_set_covered(*cp_g_in[lev], 0, cp_g_in[lev]->nComp(),
           cp_g_in[lev]->nGrow(), covered_val);
       EB_set_covered(*k_g_in[lev], 0, k_g_in[lev]->nComp(),
           k_g_in[lev]->nGrow(), covered_val);
     }

  }
}
