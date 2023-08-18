#include <AMReX.H>
#include <AMReX_Print.H>

#include <pm_fluid.H>

using namespace amrex;

pm_fluid::
pm_fluid ( std::string a_plot_file_name )
  : m_verbose(0)
  , m_error(0)
  , m_pf_name(a_plot_file_name)
  , m_pf(a_plot_file_name)
  , m_pf_size(0)
{

  if (m_pf.finestLevel() != 0) {
    amrex::Print() << "The tool is only setup for single level.\n";
    m_error = 1; return;
  }

  m_pf_size = m_pf.nComp();

  m_pf_variables = m_pf.varNames();

  for (int lc(0); lc < m_pf_size; ++lc) {
    m_variable_map[m_pf_variables[lc]] = lc;
  }
}

void
pm_fluid::
show_map ( )
{
  Print() << "\n Fluid variable count: " << m_pf_size << "\n";
  for (int lc(0); lc<m_pf_size; ++lc) {
    Print() << "  Variable: " << std::setw(12) << std::setfill(' ') << std::left << m_pf_variables[lc]
            << "  Mapped index: " <<  m_variable_map[m_pf_variables[lc]] << "\n";
  }
  Print() << "\n";
}


#if 0
  if ( a_input_size_fluid > 0 ) {

    if(!m_variable_map_fluid.count("ep_g")) {
      Print() << "Volume fraction (ep_g) not found in plot file!\n";
      m_error = 1; return;
    }

    for (int lc(0); lc < m_input_variables_fluid.size(); lc++) {
      if(!m_variable_map_fluid.count(m_input_variables_fluid[lc])) {
        Print() << m_input_variables_fluid[lc] << " not found in plot file!\n";
        m_error = 1; return;
      } else if (m_verbose) {
        Print() << "Found " << m_input_variables_fluid[lc] << " in plot file.\n";
      }
    }
  } else if (m_verbose) {
    Print() << "No fluid variable names provided!\n";
  }
#endif
