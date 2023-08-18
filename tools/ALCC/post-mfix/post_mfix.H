#ifndef POST_MFIX_
#define POST_MFIX_

#include <AMReX_PlotFileUtil.H>

#include <pm_fluid.H>
#include <pm_particles.H>

class post_mfix
{
  public:

    post_mfix (std::string a_plot_file_name);

    ~post_mfix ();

    int error () noexcept { return m_error; }

    //pm_fluid* fluid () noexcept { return &m_fluid; }
    //pm_particles* particles () noexcept { return &m_particles; }

    int size_fluid () noexcept { return m_input_size_fluid; }

    std::string variable_fluid (int const a_lc) noexcept
    { return m_input_variables_fluid[a_lc]; }

    int size_particles () noexcept { return m_input_size_particles; }

    std::string variable_particles (int const a_lc) noexcept
    { return m_input_variables_particles[a_lc]; }

    int size () noexcept { return m_input_size_fluid + m_input_size_particles; }

    void print_vars ();

    void report ();

  private:

    int const m_lev = 0;

    int m_verbose;

    int m_error;

    int m_print_vars;

    int m_no_work;

    std::string m_plot_file_name;

    pm_fluid m_fluid;
    pm_particles m_particles;

    // Input arguments
    int m_input_size_fluid;
    amrex::Vector< std::string > m_input_variables_fluid;

    int m_input_size_particles;
    amrex::Vector< std::string > m_input_variables_particles;

    // Function to convert comma separated list to vector
    amrex::Vector<std::string> input_list_to_vector( std::string const a_list );

};
#endif
