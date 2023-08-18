#ifndef POST_MFIX_FLUID_
#define POST_MFIX_FLUID_

#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

class pm_fluid
{
  public:

    pm_fluid ( std::string a_pf_name );

    int error () noexcept { return m_error; }

    void set_verbose ( int const a_verbose) noexcept { m_verbose = a_verbose; }

    void show_map ();

    int size () noexcept { return m_pf_size; }

    int contains ( std::string a_var ) noexcept { return m_variable_map.count(a_var); }

    int index ( std::string a_var ) noexcept {
      if (!m_variable_map.count(a_var)) {
        std::string message;
        message = "Variable " + a_var + " not found in particle plot file!";
        amrex::Print() << message << "\n";
        amrex::Abort(message);
      }
      return m_variable_map[a_var];
    }

    // Total number of grid cells in the domain.
    int ncells () noexcept {
      amrex::Box domain = m_pf.probDomain(m_lev);
      return domain.length(0)*domain.length(1)*domain.length(2);
    }

    amrex::MultiFab get () noexcept { return m_pf.get(m_lev); }

  private:

    int m_verbose;
    int m_error;

    int const m_lev = 0;

    std::string m_pf_name;

    amrex::PlotFileData m_pf;

    int m_pf_size;
    amrex::Vector<std::string> m_pf_variables;
    std::map<std::string, int> m_variable_map;

};

#endif

