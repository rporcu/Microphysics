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
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                            amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS,
                            amrex::Real fvelx, amrex::Real fvely, amrex::Real fvelz,
                            int i, int j, int k, int pid) const
    {
        amrex::Real RE = (Mug > 0.0) ? DPM*vrel*ROPg/Mug : large_number;

        amrex::Real C_d;
        if (RE <= 1000.0)
        {
            C_d = (24.0/(RE+small_number)) * (1.0 + 0.15*std::pow(RE, 0.687));
        }
        else
        {
            C_d = 0.44;
        }

        if (RE < eps) return 0.0;
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
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                            amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS,
              amrex::Real fvelx, amrex::Real fvely, amrex::Real fvelz,
              int i, int j, int k, int pid) const
    {
  amrex::Real ROg = ROPg / EPg;

        amrex::Real RE = (Mug > 0.0) ? DPM*vrel*ROPg/Mug : large_number;

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
  if (RE < eps) return 0.0;
  return (1.0 - PHI_gs)*Ergun + PHI_gs*WenYu;
    }
};

/**

   Subroutine: DRAG_BVK2
   Purpose: Calculate the gas-solids drag coefficient

   Literature/Document References:
       Beetstra, van der Hoef, Kuipers, Chem. Eng. Science 62
     (Jan 2007)

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
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                            amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS,
              amrex::Real fvelx, amrex::Real fvely, amrex::Real fvelz,
              int i, int j, int k, int pid) const
    {
        amrex::Real RE = (Mug > 0.0) ? DPA*vrel*ROPg/Mug : large_number;

        if (RE > eps)
        {
      amrex::Real oEPgfour = 1.0/EPg/EPg/EPg/EPg;

      // eq(9) BVK J. fluid. Mech. 528, 2005
      // (this F_Stokes is /= of Koch_Hill by a factor of ep_g)
            amrex::Real F_Stokes = 18.0*Mug*EPg/DPM/DPM;

            amrex::Real F = 10.0*PHIS/EPg/EPg + EPg*EPg*(1.0 + 1.5*sqrt(PHIS));

      F += RE*(0.11*PHIS*(1.0+PHIS) - 4.56e-3*oEPgfour +
           std::pow(RE, -0.343)*(0.169*EPg + 6.44e-2*oEPgfour));

            // F += 0.413*RE/(24.0*EPg*EPg) *
            //     (1.0/EPg + 3.0*EPg*PHIS + 8.4/std::pow(RE, 0.343)) /
            //     (1.0 + std::pow(10.0, 3.0*PHIS)/std::pow(RE, 0.5 + 2.0*PHIS));

            return F*F_Stokes;
        }
        else
        {
            return 0.0;
        }
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
    AMREX_GPU_EXTERNAL
    AMREX_GPU_HOST_DEVICE
    amrex::Real operator() (amrex::Real EPg, amrex::Real Mug, amrex::Real ROPg, amrex::Real vrel,
                            amrex::Real DPM, amrex::Real DPA, amrex::Real PHIS,
                            amrex::Real fvelx, amrex::Real fvely, amrex::Real fvelz,
                            int i, int j, int k, int pid) const;
};
