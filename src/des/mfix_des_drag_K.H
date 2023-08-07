#ifndef MFIX_DES_DRAG_K_H_
#define MFIX_DES_DRAG_K_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

#include <cmath>

struct ComputeDragWenYu
{
    amrex::Real small_number, large_number, eps;
    ComputeDragWenYu (amrex::Real small_, amrex::Real large_, amrex::Real eps_)
        : small_number(small_), large_number(large_), eps(eps_) {}

     // Literature/Document References:                                     C
     // Wen CY, Yu YH (1966). Chemical Engineering Progress Symposium    C
     //    Series 62: 100-111.

  /**
     *  \brief Returns: the calculated drag coefficient.
     *
     *  Inputs:
     *    EPg     - gas volume fraction
     *    Mug     - gas laminar viscosity
     *    ROpg    - gas density * EP_g
     *    vrel    - magnitude of gas-solids relative velocity
     *    DPM     - particle diameter of solids phase M
     *    DPA     - average particle diameter
     *    PHIS    - solids volume fraction of solids phases
     *    fvelx   - x component of the fluid velocity at the particle position
     *    fvely   - y component of the fluid velocity at the particle position
     *    fvelz   - z component of the fluid velocity at the particle position
     *    i, j, k - particle cell indices
     *    pid     - particle id number
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                            amrex::Real DPM, amrex::Real /*DPA*/, amrex::Real /*PHIS*/,
                            amrex::Real /*fvelx*/, amrex::Real /*fvely*/, amrex::Real /*fvelz*/,
                            int /*i*/, int /*j*/, int /*k*/, int /*pid*/) const
    {
        amrex::Real RE = (Mug > 0.0) ? DPM*vrel*ROPg/Mug : large_number;
        if (RE < eps) return 0.0;

        amrex::Real C_d;
        if (RE <= 1000.0)
        {
            C_d = (24.0/(RE+small_number)) * (1.0 + 0.15*std::pow(RE, 0.687));
        }
        else
        {
            C_d = 0.44;
        }

        return 0.75 * C_d * vrel * ROPg * std::pow(EPg, -2.65) / DPM;
    }
};


struct ComputeDragGidaspow
{
    amrex::Real small_number, large_number, eps;
    ComputeDragGidaspow (amrex::Real small_, amrex::Real large_, amrex::Real eps_)
        : small_number(small_), large_number(large_), eps(eps_) {}

    // Literature/Document References:
    //  Ding J, Gidaspow D (1990). AIChE Journal 36: 523-538.
    //  original source unknown:
    //  Lathouwers D, Bellan J (2000). Proceedings of the 2000 U.S. DOE
    //     Hydrogen Program Review NREL/CP-570-28890. Available from
    //  http://www.eere.energy.gov/hydrogenandfuelcells/pdfs/28890k.pdf.

    /**
     *  \brief Returns: the calculated drag coefficient.
     *
     *  Inputs:
     *    EPg     - gas volume fraction
     *    Mug     - gas laminar viscosity
     *    ROpg    - gas density * EP_g
     *    vrel    - magnitude of gas-solids relative velocity
     *    DPM     - particle diameter of solids phase M
     *    DPA     - average particle diameter
     *    PHIS    - solids volume fraction of solids phases
     *    fvelx   - x component of the fluid velocity at the particle position
     *    fvely   - y component of the fluid velocity at the particle position
     *    fvelz   - z component of the fluid velocity at the particle position
     *    i, j, k - particle cell indices
     *    pid     - particle id number
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                            amrex::Real DPM, amrex::Real /*DPA*/, amrex::Real /*PHIS*/,
                            amrex::Real /*fvelx*/, amrex::Real /*fvely*/, amrex::Real /*fvelz*/,
                            int /*i*/, int /*j*/, int /*k*/, int /*pid*/) const
    {
        amrex::Real ROg = ROPg / EPg;

        amrex::Real RE = (Mug > 0.0) ? DPM*vrel*ROPg/Mug : large_number;
        if (RE < eps) return 0.0;

        // Dense phase - EPg <= 0.8
        amrex::Real Ergun = 150.0*(1.0 - EPg)*Mug / (EPg*DPM*DPM) + 1.75*ROg*vrel/DPM;

        // Dilute phase - EPg > 0.8
        amrex::Real C_d;
        if (RE <= 1000.0)
        {
          C_d = (24.0/(RE+small_number)) * (1.0 + 0.15*std::pow(RE, 0.687));
        }
        else
        {
          C_d = 0.44;
        }
        amrex::Real WenYu = 0.75*C_d*vrel*ROPg*std::pow(EPg, -2.65) / DPM;

        // switch function
        amrex::Real PHI_gs = std::atan(150.0*1.75*(EPg - 0.8))/M_PI + 0.5;

        // blend the models
        return (1.0 - PHI_gs)*Ergun + PHI_gs*WenYu;
    }
};


/**

   Subroutine: DRAG_BVK2
   Purpose: Calculate the gas-solids drag coefficient

   Literature/Document References:
   van der Hoef, Beetstra, Kuipers, J. Fluid. Mech., 528 (2005)

   Beetstra, van der Hoef, Kuipers, Chem. Eng. Science, 62
     (Jan 2007)

   Tenneti, Garg, Subramaniam, Int. J. Multiphase Flow, 37(9)
     (2011)  

   Tang, Peters, Kuipers, Kriebitzsch, van der Hoef, AIChEJ,
       61(2) (Feb 2015)

 */
struct ComputeDragBVK2
{
    amrex::Real small_number, large_number, eps;
    ComputeDragBVK2 (amrex::Real small_, amrex::Real large_, amrex::Real eps_)
        : small_number(small_), large_number(large_), eps(eps_) {}

    // Literature/Document References:
    //  Beetstra, van der Hoef, Kuipers, Chem. Eng. Science 62
    //  (Jan 2007)
    //
    //  Tang, Peters, Kuipers, Kriebitzsch, van der Hoef, AIChEJ,
    //  61(2) (Feb 2015)

    /**
     *  \brief Returns: the calculated drag coefficient.
     *
     *  Inputs:
     *    EPg     - gas volume fraction
     *    Mug     - gas laminar viscosity
     *    ROpg    - gas density * EP_g
     *    vrel    - magnitude of gas-solids relative velocity
     *    DPM     - particle diameter of solids phase M
     *    DPA     - average particle diameter
     *    PHIS    - solids volume fraction of solids phases
     *    fvelx   - x component of the fluid velocity at the particle position
     *    fvely   - y component of the fluid velocity at the particle position
     *    fvelz   - z component of the fluid velocity at the particle position
     *    i, j, k - particle cell indices
     *    pid     - particle id number
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                            amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS,
                            amrex::Real /*fvelx*/, amrex::Real /*fvely*/, amrex::Real /*fvelz*/,
                            int /*i*/, int /*j*/, int /*k*/, int /*pid*/) const
    {
        amrex::Real RE = (Mug > 0.0) ? DPA*vrel*ROPg/Mug : large_number;
        if (RE < eps) return 0.0;

        amrex::Real oEPgfour = 1.0/EPg/EPg/EPg/EPg;

        // (this F_Stokes is /= of Koch_Hill by a factor of ep_g)
        // see eq(5) J. Fluid. Mech. 528, 2005
        amrex::Real F_Stokes = 18.0*Mug*EPg/DPM/DPM;

        // BVK2:
        amrex::Real F = 10.0*PHIS/EPg/EPg + EPg*EPg*(1.0 + 1.5*sqrt(PHIS));

        F += RE*(0.11*PHIS*(1.0+PHIS) - 4.56e-3*oEPgfour +
             std::pow(RE, -0.343)*(0.169*EPg + 6.44e-2*oEPgfour));

        // BVK: 
        // amrex::Real F = 10.0*PHIS/EPg/EPg + EPg*EPg*(1.0 + 1.5*sqrt(PHIS));
        // F += 0.413*RE/(24.0*EPg*EPg) *
        //     (1.0/EPg + 3.0*EPg*PHIS + 8.4/std::pow(RE, 0.343)) /
        //     (1.0 + std::pow(10.0, 3.0*PHIS)/std::pow(RE, 0.5 + 2.0*PHIS));

        //TGS: 
        // amrex::Real F = (amrex::max(1. + 0.15*std::pow(RE, 0.687), 0.44*RE/24.)
        //                 + 5.81*PHIS + 0.48*std::pow(PHIS, 1./3.)/EPg)/EPg/EPg
        //                 + EPg*PHIS*PHIS*PHIS*RE*(0.95 + 0.61*PHIS*PHIS*PHIS/EPg/EPg);

        return F*F_Stokes;
    }
};


/**********************************************************************************
*
*********************************************************************************/
struct ComputeDragSyamOBrien1988
{
    amrex::Real small_number, large_number, eps;
    amrex::Real drag_c1, drag_d1;

    ComputeDragSyamOBrien1988 (amrex::Real small_, amrex::Real large_, amrex::Real eps_,
                               amrex::Real c1_, amrex::Real d1_)
        : small_number(small_)
        , large_number(large_)
        , eps(eps_)
        , drag_c1(c1_)
        , drag_d1(d1_)
        {}

    // Literature/Document References:
    // Syamlal M, O'Brien TJ (1988). International Journal of Multiphase Flow 14: 473-481.

    /**
     *  \brief Returns: the calculated drag coefficient.
     *
     *  Inputs:
     *    EPg     - gas volume fraction
     *    Mug     - gas laminar viscosity
     *    ROpg    - gas density * EP_g
     *    vrel    - magnitude of gas-solids relative velocity
     *    DPM     - particle diameter of solids phase M
     *    DPA     - average particle diameter
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                            amrex::Real DPM, amrex::Real DPA, amrex::Real /*PHIS*/,
                            amrex::Real /*fvelx*/, amrex::Real /*fvely*/, amrex::Real /*fvelz*/,
                            int /*i*/, int /*j*/, int /*k*/, int /*pid*/) const
    {

        amrex::Real RE = (Mug > 0.0) ? DPA*vrel*ROPg/Mug : large_number;

        if (RE < eps) return 0.0;

        // Calculate V_rm
        const amrex::Real A = std::pow(EPg, 4.14);

        const amrex::Real B = (EPg <= 0.85) ?
            drag_c1*std::pow(EPg, 1.28) :
                    std::pow(EPg, drag_d1);


        const amrex::Real V_RM = 0.5*(A - 0.06*RE +
           std::sqrt(0.0036*RE*RE + 0.12*RE*(2.0*B - A) + A*A));


        // Calculate single sphere drag correlation of Dalla Valle (1948)
        // multiplied by the Reynolds number.
        amrex::Real C_DSXRE_DV = (0.63*std::sqrt(RE/V_RM) + 4.8);
        C_DSXRE_DV *= C_DSXRE_DV;

        return 0.75*Mug*EPg*C_DSXRE_DV/(V_RM*DPM*DPM);


    }
};


struct ComputeDragUser
{
    amrex::Real small_number, large_number, eps;
    ComputeDragUser (amrex::Real small_, amrex::Real large_, amrex::Real eps_)
        : small_number(small_), large_number(large_), eps(eps_) {}

    /**
     *  \brief Returns: the calculated drag coefficient.
     *
     *  Inputs:
     *    EPg     - gas volume fraction
     *    Mug     - gas laminar viscosity
     *    ROpg    - gas density * EP_g
     *    vrel    - magnitude of gas-solids relative velocity
     *    DPM     - particle diameter of solids phase M
     *    DPA     - average particle diameter
     *    PHIS    - solids volume fraction of solids phases
     *    fvelx   - x component of the fluid velocity at the particle position
     *    fvely   - y component of the fluid velocity at the particle position
     *    fvelz   - z component of the fluid velocity at the particle position
     *    i, j, k - particle cell indices
     *    pid     - particle id number
     */
    AMREX_GPU_HOST_DEVICE
#if defined(AMREX_USE_DPCPP)
    SYCL_EXTERNAL
#endif
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                            amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS,
                            amrex::Real fvelx, amrex::Real fvely, amrex::Real fvelz,
                            int i, int j, int k, int pid) const;
};

#endif
