#ifndef MFIX_PIC_PARMS_H_
#define MFIX_PIC_PARMS_H_

#include <AMReX_REAL.H>

#define MFP_DISABLED 0
#define SCALE_MFP_DISABLED 0
#define SCALE_MFP_CHAPMAN  0
#define SCALE_MFP_CUTCHIS  0
#define SCALE_MFP_MUSSER   1


class MFIXPIC {

  public:
    struct InitialStepType { enum {Invalid=-1, zero_eps, nth_eps, taylor_approx}; };

    MFIXPIC ();

    void Initialize ();

    int solve () const { return m_solve; }
    void set_solve (const int val) { m_solve = val; }

    int NPHASE () const { return m_NPHASE; }

    int verbose () const { return m_verbose; }

    // Coefficients for computing the solid stress
    amrex::Real Ps () const { return m_Ps; }
    amrex::Real beta () const { return m_beta; }
    amrex::Real ep_cp () const { return m_ep_cp; }

    amrex::Real vel_ref_frame () const { return m_vel_ref_frame; }

    amrex::Real small_number () const { return m_small_number; }
    amrex::Real damping_factor () const { return m_damping_factor; }

    amrex::Real damping_factor_wall_normal () const { return m_damping_factor_wall_normal; }
    amrex::Real damping_factor_wall_tangent () const { return m_damping_factor_wall_tangent; }

    amrex::Real advance_vel_p () const { return m_advance_vel_p; }

    int max_iter () const { return m_max_iter; }

    int initial_step () const { return m_initial_step; }

    int restart_refinement () const { return m_restart_refinement; }

  private:
    int m_solve;

    int m_NPHASE;

    int m_verbose;

    // Coefficients for computing the solid stress
    amrex::Real m_Ps;
    amrex::Real m_beta;
    amrex::Real m_ep_cp;

    amrex::Real m_vel_ref_frame;

    amrex::Real m_small_number;
    amrex::Real m_damping_factor;

    amrex::Real m_damping_factor_wall_normal;
    amrex::Real m_damping_factor_wall_tangent;

    amrex::Real m_advance_vel_p;

    int m_max_iter;

    int m_initial_step;

    int m_restart_refinement;
};

#endif
