#include <eos_mod.H>

AMREX_GPU_HOST_DEVICE
amrex::Real idealgas (const amrex::Real& mw, const amrex::Real& pg, const amrex::Real& tg)
{
  return (pg*mw)/(8314.56*tg);
}



AMREX_GPU_HOST_DEVICE
amrex::Real sutherland (const amrex::Real& tg)
{
   // Dummy arguments
   //---------------------------------------------------------------------//

   // Gas viscosity   (in Poise or Pa.s)
   // Calculating gas viscosity using Sutherland's formula with
   // Sutherland's constant (C) given by Vogel's equation C = 1.47*Tb.
   // For air  C = 110 (Tb=74.82)
   //         mu = 1.71*10-4 poise at T = 273K

   const amrex::Real tg_scaled = tg/273.0;

   return 1.7e-5 * tg_scaled * std::sqrt(tg_scaled) * (383./(tg+110.));
}
