#include <AMReX.H>
#include <AMReX_Print.H>

#include <post_mfix.H>

using namespace amrex;

// Destructor
post_mfix::~post_mfix ()
{
}

// Constructor
post_mfix::post_mfix (std::string a_plot_file_name)
  : m_verbose(0)
  , m_error(0)
  , m_print_vars(0)
  , m_no_work(1)
  , m_plot_file_name(a_plot_file_name)
  , m_fluid(a_plot_file_name)
  , m_particles(a_plot_file_name, "particles")
  , m_input_size_fluid(0)
  , m_input_size_particles(0)
{

  const int narg = amrex::command_argument_count();

  int farg = 1;
  while (farg <= narg) {
    const std::string& name = amrex::get_command_argument(farg);
    if (name == "-v" || name == "--verbose") {
      m_verbose = 1;
    } else if (name == "-f" || name == "--fluid") {
      m_input_variables_fluid = input_list_to_vector(get_command_argument(++farg));
    } else if (name == "-p" || name == "--particle") {
      m_input_variables_particles = input_list_to_vector(get_command_argument(++farg));
    } else if (name == "--print-vars") {
      m_print_vars = 1;
    } else {
      break;
    }
    ++farg;
  }

  if (m_fluid.error()) { m_error = 1; return; }
  //if (m_particles.error()) { m_error = 1; return; }

  m_input_size_fluid = m_input_variables_fluid.size();
  m_fluid.set_verbose(m_verbose);


  // Loop over the variables that the user provided and make
  // sure we have them in the plot file.
  if ( m_input_size_fluid ) {

    m_no_work = 0;

    if(!m_fluid.contains("ep_g")) {
      Print() << "Volume fraction (ep_g) not found in plot file!\n";
      m_error = 1; return;
    }

    for (int lc(0); lc < m_input_variables_fluid.size(); lc++) {
      if(!m_fluid.contains(m_input_variables_fluid[lc])) {
        Print() << m_input_variables_fluid[lc] << " not found in plot file!\n";
        m_error = 1; return;
      } else if (m_verbose) {
        Print() << "Found " << m_input_variables_fluid[lc] << " in plot file.\n";
      }
    }
  } else if (m_verbose) {
    Print() << "No fluid variable names provided!\n";
  }

  m_input_size_particles = m_input_variables_particles.size();
  m_particles.set_verbose(m_verbose);

  if ( m_input_size_particles ) {

    m_no_work = 0;

    for (int lc(0); lc < m_input_size_particles; lc++) {
      if(!m_particles.contains(m_input_variables_particles[lc])) {
        Print() << m_input_variables_particles[lc] << " not found in plot file!\n";
        m_error = 1; return;
      } else if (m_verbose) {
        Print() << "Found " << m_input_variables_particles[lc] << " in plot file.\n";
      }
    }
  } else if (m_verbose) {
    Print() << "No particle variable names provided!\n";
  }

}


amrex::Vector<std::string>
post_mfix::
input_list_to_vector( std::string const a_list )
{
  // Count the number of variables by counting commas.
  int const size = 1 + std::count(a_list.begin(), a_list.end(), ',');

  amrex::Vector<std::string> var_vect(size);
  std::istringstream var_stream{a_list};

  int lc(0);
  for(std::string var; std::getline(var_stream, var,','); ) {
    var_vect[lc++] = var;
  }
  return var_vect;
}

void
post_mfix::
print_vars ()
{
  if (m_print_vars) {
    m_fluid.show_map();
    m_particles.show_map();
  }

}
