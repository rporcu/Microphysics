#ifndef POST_MFIX_PARTICLES_
#define POST_MFIX_PARTICLES_

#include <AMReX_Gpu.H>
#include <AMReX_Print.H>

class pm_particles {

  public:

    pm_particles ( const std::string& a_plt_file,
                   const std::string& a_particle_type);

    void set_verbose ( int const a_verbose ) noexcept { m_verbose = a_verbose; }

    void show_map ();

    int contains ( std::string a_var ) noexcept {
      return (m_variable_map_r.count(a_var)
            + m_variable_map_i.count(a_var));
    }

    int index ( std::string a_var ) noexcept {
      if (!m_variable_map_r.count(a_var)) {
        std::string message;
        message = "Variable " + a_var + " not found in particle plot file!\n";
        amrex::Abort(message);
      }
      return m_variable_map_r[a_var];
    }

    int size () noexcept { return m_pf_size; }

    int rank_particle_count () noexcept { return m_np_all_grids; }
    int particle_count () noexcept { return m_nparticles; }

    amrex::Real sum ( std::string a_var,
                      bool a_local = false);

    amrex::Real fluct ( std::string const a_var,
                        amrex::Real const a_avg,
                        bool a_local = false );


    //amrex::Gpu::DeviceVector<amrex::Real>& get ( std::string a_var )
    amrex::Real* get ( std::string a_var )
    {
      int var_index = index(a_var);
      get_pf_data_r(var_index);

      return m_d_data_r.data();
    }

  private:

    void get_pf_data_r ( int const a_index );

    std::string const m_directions{"xyz"};
    int const m_lev = 0;

    int m_last_var_read;

    int m_verbose;

    std::string m_pf_name;

    int m_np_all_grids;

    std::string m_pf_version;

    // dimension
    int m_pf_dim;

    int m_pf_size_r;       // number of reals, including position
    int m_pf_size_r_extra; // number of reals, excluding position

    // Variable names (reals)
    std::vector<std::string> m_pf_variables_r;

    int m_pf_size_i;       // number of ints, excluding id and cpu

    // Variable names (ints)
    std::vector<std::string> m_pf_variables_i;

    bool m_is_checkpoint;

    int m_nparticles;
    int m_next_id;
    int m_finest_level;

    std::vector<int> m_size_grids;

    std::vector<std::vector<int> > m_file_nums;
    std::vector<std::vector<int> > m_particle_counts;
    std::vector<std::vector<int> > m_offsets;

    // This is additional metadata derived from the above
    int m_pf_size;

    std::string hdr_file_name;
    //std::vector<std::string> comp_names;

    std::map<std::string, int> m_variable_map_r;
    std::map<std::string, int> m_variable_map_i;

    int m_ibegin;
    int m_iend;

    amrex::Gpu::HostVector<amrex::Real>   m_h_data_r;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_data_r;

    std::string getDataFileName(const std::string& a_prefix, int a_level, int a_file_num);

    void getDataBuffer( std::vector<char>& a_buffer, const std::string& a_file,
                        size_t a_buffer_size, int a_offset);

    // These operators only use the data actually in the headers, not the derived stuff.
    //friend std::ostream& operator<< (std::ostream& a_stream, const pm_particles& a_header);
    friend std::istream& operator>> (std::istream& a_stream, pm_particles& header);
};

#endif
