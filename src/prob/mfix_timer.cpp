#include <mfix_timer.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>
#include <fstream>

using namespace amrex;

void
MFIXTimer::Initialize ()
{
  ParmParse pp("mfix");

  pp.query("fixed_dt", m_dt);

  if(m_dt > 0.)
    m_timestep_type = TimestepType::Fixed;
  else
    m_timestep_type = TimestepType::Dynamic;

  pp.query("dt_min", m_dt_min);
  pp.query("dt_max", m_dt_max);

  std::string walltime_limit_in;
  if (pp.query("walltime_limit", walltime_limit_in)) {
    int HH(0), MM(0), SS(0);
    if (sscanf(walltime_limit_in.c_str(), "%d:%d:%d", &HH, &MM, &SS) >= 2) {
      m_walltime_limit = static_cast<Real>(HH*3600 + MM*60 + SS);
    } else {
      std::string message =
        " Error: Unable to correctly parse walltime limit "
        + walltime_limit_in + "\n" + " The correct format is HH:MM:SS\n";
      Print() << message;
      Abort(message);
    }
  }

  pp.query("stop_time", m_stop_time);
  pp.query("overstep_end_time", m_overstep_end_time);
  pp.query("max_step", m_max_step);
  m_usr_max_step = m_max_step; // save the original max step selected by the user

  pp.query("clean_exit", m_clean_exit);
}


void
MFIXTimer::reset (const MFIXTimer& other)
{
  m_runtime_start = other.runtime_start();
  m_walltime_limit = other.walltime_limit();
  m_avg_step_runtime = other.avg_step_runtime();
  m_start_time = other.start_time();
  m_time = other.time();
  m_stop_time = other.stop_time();
  m_timestep_type = other.timestep_type();
  m_dt = other.dt();
  m_dt_min = other.dt_min();
  m_dt_max = other.dt_max();
  m_first_step = other.first_step();
  m_nstep = other.nstep();
  m_max_step = other.max_step();
  m_usr_max_step = other.usr_max_step();
  m_clean_exit = other.clean_exit();
  m_run_status_type = other.run_status_type();
}


int
MFIXTimer::ok ()
{
  if ((m_stop_time >= 0.) && ((m_time+0.1*m_dt) >= m_stop_time))
    m_run_status_type = MFIXRunStatusType::TimeIsOver;

  if ((m_max_step >= 0) && (m_nstep >= m_max_step))
    m_run_status_type = MFIXRunStatusType::IsFinalStep;

  if (ParallelDescriptor::IOProcessor()) {

    if ((m_walltime_limit > 0.) && (!runtime_left_is_sufficient()))
      m_run_status_type = MFIXRunStatusType::RuntimeIsOver;

    if((m_clean_exit != "") && std::ifstream(m_clean_exit.c_str()).good())
      m_run_status_type = MFIXRunStatusType::UserStop;
  }

  if (((m_walltime_limit > 0.)) || (m_clean_exit != "")) {
    ParallelDescriptor::Bcast(&m_run_status_type, 1, ParallelDescriptor::IOProcessorNumber());
  }

  return (m_run_status_type == MFIXRunStatusType::OK);
}


int
MFIXTimer::runtime_left_is_sufficient ()
{
  const Real needed_time = 1.2*(m_max_write_chkpt_time + m_avg_step_runtime);

  const Real missing_time = m_walltime_limit - elapsed_runtime();

  if (needed_time < missing_time)
    return 1;

  return 0;
}
