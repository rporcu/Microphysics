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

  std::string walltime_buffer_in;
  if (pp.query("walltime_buffer", walltime_buffer_in)) {
    int HH(0), MM(0), SS(0);
    if (sscanf(walltime_buffer_in.c_str(), "%d:%d:%d", &HH, &MM, &SS) >= 2) {
      m_walltime_buffer = static_cast<Real>(HH*3600 + MM*60 + SS);
    } else {
      std::string message =
        " Error: Unable to correctly parse walltime buffer "
        + walltime_buffer_in + "\n" + " The correct format is HH:MM:SS\n";
      Print() << message;
      Abort(message);
    }
  }

  pp.query("stop_time", m_stop_time);
  pp.query("overstep_end_time", m_overstep_end_time);
  pp.query("max_step", m_max_step);

  pp.query("clean_exit", m_clean_exit);
}


int
MFIXTimer::ok ()
{
  if ((m_stop_time >= 0.) && ((m_time+0.1*m_dt) >= m_stop_time))
    m_run_status_type = MFIXRunStatusType::TimeIsOver;

  if ((m_max_step >= 0) && (m_nstep >= m_max_step))
    m_run_status_type = MFIXRunStatusType::IsFinalStep;

  if (ParallelDescriptor::IOProcessor()) {

    if ((m_walltime_limit > 0.) && ((elapsed_runtime()+m_avg_step_runtime) > m_walltime_limit))
      m_run_status_type = MFIXRunStatusType::RuntimeIsOver;

    if((m_clean_exit != "") && std::ifstream(m_clean_exit.c_str()).good())
      m_run_status_type = MFIXRunStatusType::UserStop;
  }

  if (((m_walltime_limit > 0.)) || (m_clean_exit != "")) {

    ParallelDescriptor::Bcast(&m_run_status_type, 1, ParallelDescriptor::IOProcessorNumber());
  }

  return (m_run_status_type == MFIXRunStatusType::OK);
}
