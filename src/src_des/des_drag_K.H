#include <AMReX_REAL.H>
#include <MFIX_DEM_Parms.H>

/**
 *  \brief Function for computing drag according to the BVK drag model.
 *  
 *  Returns: the calculated drag coefficient.
 *
 *  Inputs:
 *    EPg  - gas volume fraction
 *    Mug  - gas laminar viscosity
 *    ROpg - gas density * EP_g
 *    vrel - magnitude of gas-solids relative velocity
 *    DPM  - particle diamater of solids phase M
 *    DPA  - average particle diameter
 *    PHIS - solids volume fraction of solids phases
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real compute_drag_bvk (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                              amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS)
{
    amrex::Real lDga;  // result - the drag coefficient.
    
    amrex::Real RE;  // Reynolds number
    if (Mug > 0.0)
    {
        RE = DPA*vrel*ROPg/Mug;
    }
    else
    {
        RE = DEMParams::large_number;
    }

    if (RE > DEMParams::eps)
    {
        amrex::Real F_Stokes = 18.0*Mug*EPg/DPM/DPM;
        amrex::Real F = 10.0*PHIS/EPg/EPg + EPg*EPg*(1.0 + 1.5*sqrt(PHIS));
        F += 0.413*RE/(24.0*EPg*EPg) *
            (1.0/EPg + 3.0*EPg*PHIS + 8.4/std::pow(RE, 0.343)) /
            (1.0 + std::pow(10.0, 3.0*PHIS)/std::pow(RE, 0.5 + 2.0*PHIS));

        lDga = F*F_Stokes;
    }
    else
    {
        lDga = 0.0;
    }

    return lDga;
}

/**
 *  \brief Function for computing drag according to the Koch-Hill drag model.
 *  
 *  Returns: the calculated drag coefficient.
 *
 *  Inputs:
 *    EPg  - gas volume fraction
 *    Mug  - gas laminar viscosity
 *    ROpg - gas density * EP_g
 *    vrel - magnitude of gas-solids relative velocity
 *    DPM  - particle diamater of solids phase M
 *    DPA  - average particle diameter
 *    PHIS - solids volume fraction of solids phases
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real compute_drag_koch_hill (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg,
                                    amrex::Real vrel, amrex::Real DPM, amrex::Real DPA,
                                    amrex::Real PHIS)
{
    amrex::Abort("Not implemented");

    return 0.0;
}
