#ifndef _MFIX_TIMER_H_
#define _MFIX_TIMER_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParallelDescriptor.H>

#include <limits>
#include <string>

struct MFIXRunStatusType { enum {OK, TimeIsOver, IsFinalStep, RuntimeIsOver, UserStop}; };

class MFIXTimer
{
  public:
    struct TimestepType { enum {Invalid=-1, Fixed, Dynamic}; };

    // Constructor
    MFIXTimer ()
      : m_runtime_start(0.)
      , m_walltime_limit(-1.)
      , m_avg_step_runtime(0.)
      , m_max_write_chkpt_time(0.)
      , m_start_time(0.)
      , m_time(0.)
      , m_stop_time(std::numeric_limits<amrex::Real>::infinity())
      , m_timestep_type(TimestepType::Invalid)
      , m_dt(-1.)
      , m_dt_min(0.)
      , m_dt_max(1.e14)
      , m_first_step(0)
      , m_nstep(0)
      , m_max_step(-1)
      , m_usr_max_step(-1)
      , m_clean_exit("")
      , m_run_status_type(MFIXRunStatusType::OK)
    {
      m_start_time = system_time();
    }

    void Initialize ();

    void reset (const MFIXTimer& other);

    amrex::Real runtime_start () const { return m_runtime_start; }
    amrex::Real walltime_limit () const { return m_walltime_limit; }
    amrex::Real& walltime_limit () { return m_walltime_limit; }
    amrex::Real avg_step_runtime () const { return m_avg_step_runtime; }
    amrex::Real max_write_chkpt_time () const { return m_max_write_chkpt_time; }

    AMREX_FORCE_INLINE
    void update_avg_step_runtime (const amrex::Real step_time)
    {
      m_avg_step_runtime = (m_avg_step_runtime*(m_nstep-1) + step_time) / m_nstep;
      amrex::ParallelDescriptor::ReduceRealMax(&m_avg_step_runtime, 1);

      if ((m_max_step > 0) && (m_walltime_limit > 0.))
        this->update_max_step();
    }

    AMREX_FORCE_INLINE
    void update_max_write_chkpt_time (const amrex::Real write_time)
    {
      m_max_write_chkpt_time = amrex::max(m_max_write_chkpt_time, write_time);
      amrex::ParallelDescriptor::ReduceRealMax(&m_max_write_chkpt_time, 1);
    }

    AMREX_FORCE_INLINE
    amrex::Real system_time () const { return amrex::ParallelDescriptor::second(); }

    AMREX_FORCE_INLINE
    amrex::Real elapsed_runtime () const
    { return elapsed_runtime(m_start_time); }

    AMREX_FORCE_INLINE
    amrex::Real elapsed_runtime (const amrex::Real start) const
    { return (system_time() - start); }

    amrex::Real start_time () const { return m_start_time; }
    amrex::Real time () const { return m_time; }
    amrex::Real& time () { return m_time; }
    amrex::Real stop_time () const { return m_stop_time; }
    int overstep_end_time () const { return m_overstep_end_time; }

    AMREX_FORCE_INLINE
    void advance_time (const amrex::Real dt)
    { m_time += dt; }

    int timestep_type () const { return m_timestep_type; }
    amrex::Real dt () const { return m_dt; }
    amrex::Real& dt () { return m_dt; }
    amrex::Real dt_min () const { return m_dt_min; }
    amrex::Real dt_max () const { return m_dt_max; }

    int first_step () const { return m_first_step; }
    int nstep () const { return m_nstep; }
    int& nstep () { return m_nstep; }
    int max_step () const { return m_max_step; }
    int usr_max_step () const { return m_usr_max_step; }

    AMREX_FORCE_INLINE
    void advance_nstep (const int n)
    { m_nstep += n; }

    const std::string& clean_exit () const { return m_clean_exit; }

    int run_status_type () const { return m_run_status_type; }

    int ok ();
    int runtime_left_is_sufficient ();

  private:
    amrex::Real m_runtime_start;
    amrex::Real m_walltime_limit;
    amrex::Real m_avg_step_runtime;
    amrex::Real m_max_write_chkpt_time;

    amrex::Real m_start_time;
    amrex::Real m_time;
    amrex::Real m_stop_time;
    int m_overstep_end_time = 0;

    int m_timestep_type;
    amrex::Real m_dt;
    amrex::Real m_dt_min;
    amrex::Real m_dt_max;

    int m_first_step;
    int m_nstep;
    int m_max_step;
    int m_usr_max_step;

    std::string m_clean_exit;

    int m_run_status_type;

    AMREX_FORCE_INLINE
    void update_max_step ()
    { 
      if (m_walltime_limit > 0. && m_avg_step_runtime > 0.) {
        m_max_step = std::min(m_usr_max_step, int(amrex::Math::floor(m_walltime_limit / m_avg_step_runtime)));
      }
    }
};

#endif
