#ifndef _MFIX_DES_CONV_COEFF_K_H_
#define _MFIX_DES_CONV_COEFF_K_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>


struct NullConvectionCoeff
{
    // Constructor
    NullConvectionCoeff () {};

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real operator() (amrex::Real /*EPg*/, amrex::Real /*Mug*/, amrex::Real /*Kg*/, amrex::Real /*Cpg*/,
                            amrex::Real /*ROPg*/, amrex::Real /*Vrel*/, amrex::Real /*dp*/,
                            int /*i*/, int /*j*/, int /*k*/, int /*pid*/) const
    { return 0; }
};


struct ComputeConvRanzMarshall
{
    amrex::Real small_number;
    amrex::Real large_number;
    amrex::Real eps;

    // Constructor
    ComputeConvRanzMarshall (amrex::Real small_, amrex::Real large_, amrex::Real eps_)
        : small_number(small_), large_number(large_), eps(eps_) {}

    /*
     *  \brief Returns: the calculated heat transfer coefficient.
     *
     *  Inputs:
     *    EPg     - gas volume fraction
     *    Mug     - gas laminar viscosity
     *    Kg      - gas thermal conductivity
     *    Cpg     - gas specific heat
     *    ROPg    - gas density * EP_g
     *    Vrel    - magnitude of gas-solids relative velocity
     *    dp      - particle diameter of solids phase M
     *    i, j, k - particle cell indices
     *    pid     - particle id number
     *
     * Literature/Document References:
     * Ranz, W.E. and Marshall, W.R., "Friction and transfer coefficients for
     *     single particles and packed beds," Chemical Engineering Science,
     *     Vol. 48, No. 5, pp 247-253, 1925.
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real operator() (amrex::Real /*EPg*/, amrex::Real Mug, amrex::Real Kg, amrex::Real Cpg,
                            amrex::Real ROPg, amrex::Real Vrel, amrex::Real dp,
                            int /*i*/, int /*j*/, int /*k*/, int /*pid*/) const
    {
        // Calculate the Reynolds Number
        amrex::Real N_Re = (Mug > 0.0) ? dp*Vrel*ROPg/Mug : large_number;

        // Calculate the Prandtl Number
        amrex::Real N_Pr = (Kg  > 0.0) ? Cpg*Mug/Kg : large_number;

        // Calculate the Nusselt Number
        amrex::Real N_Nu = 2.0 + 0.6 * std::sqrt(N_Re) * std::pow(N_Pr, 0.333);

        return N_Nu * Kg / dp;
    }
};


struct ComputeConvGunn
{
    amrex::Real small_number;
    amrex::Real large_number;
    amrex::Real eps;

    // Constructor
    ComputeConvGunn (amrex::Real small_, amrex::Real large_, amrex::Real eps_)
        : small_number(small_), large_number(large_), eps(eps_) {}

    /*
     *  \brief Returns: the calculated heat transfer coefficient.
     *
     *  Inputs:
     *    EPg     - gas volume fraction
     *    Mug     - gas laminar viscosity
     *    Kg      - gas thermal conductivity
     *    Cpg     - gas specific heat
     *    ROPg    - gas density * EP_g
     *    Vrel    - magnitude of gas-solids relative velocity
     *    dp      - particle diameter of solids phase M
     *    i, j, k - particle cell indices
     *    pid     - particle id number
     *
     * Literature/Document References:
     * Ranz, W.E. and Marshall, W.R., "Friction and transfer coefficients for
     *     single particles and packed beds," Chemical Engineering Science,
     *     Vol. 48, No. 5, pp 247-253, 1925.
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real Kg, amrex::Real Cpg,
                            amrex::Real ROPg, amrex::Real Vrel, amrex::Real dp,
                            int /*i*/, int /*j*/, int /*k*/, int /*pid*/) const
    {
        // Calculate the Reynolds Number
        amrex::Real N_Re = (Mug > 0.0) ? dp*Vrel*ROPg/Mug : large_number;

        // Calculate the Prandtl Number
        amrex::Real N_Pr = (Kg  > 0.0) ? Cpg*Mug/Kg : large_number;

        // Calculate the Nusselt Number
        amrex::Real N_Nu = 
          (7 - 10*EPg + 5*EPg*EPg)*(1 + .7*std::pow(N_Re, 0.2)*std::cbrt(N_Pr))
          + (1.33 - 2.4*EPg + 1.2*EPg*EPg)*std::pow(N_Re, 0.7)*std::cbrt(N_Pr);


        return N_Nu * Kg / dp;
    }
};

#endif
