#include <eos_mod.hpp>

AMREX_GPU_HOST_DEVICE
amrex::Real sutherland(const amrex::Real& tg)
{
   // Dummy arguments
   //---------------------------------------------------------------------//

   // Gas viscosity   (in Poise or Pa.s)
   // Calculating gas viscosity using Sutherland's formula with
   // Sutherland's constant (C) given by Vogel's equation C = 1.47*Tb.
   // For air  C = 110 (Tb=74.82)
   //         mu = 1.71*10-4 poise at T = 273K

   return 1.7e-5 * std::pow(tg/273.0, 1.5) * (383./(tg+110.));
}
