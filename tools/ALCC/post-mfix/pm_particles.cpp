#include <AMReX_Reduce.H>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <cstring>
#include <cmath>
#include <cstdlib>

#ifdef BL_USE_MPI
#include <mpi.h>
#endif

#include <pm_particles.H>

using namespace amrex;

///
/// Initializes the struct to store the header information associated with
/// the given particle header file. For example, if plt_file = "plt00000"
/// and particle_type = "Tracer", this will read the Header file at
/// "plt00000/Tracer/Header"
///
pm_particles::
pm_particles(const std::string& a_plot_file,
             const std::string& a_particle_type)
  : m_verbose(0)
  , m_pf_name( a_plot_file + "/" + a_particle_type)
  , m_last_var_read(-1)
  , m_np_all_grids(0)
{


  hdr_file_name = m_pf_name + "/Header";

  std::ifstream file(hdr_file_name.c_str());
  if ( ! file.is_open() ) {
#ifdef BL_USE_MPI
    int myproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
#else
    int myproc = 0;
#endif
    if ( myproc == 0) {
      std::cout << "Error! Could not open file " << hdr_file_name << "." << std::endl;
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif
    exit(1);
  }

  file >> *this;
  if (! file.is_open() ) {
#ifdef BL_USE_MPI
    int myproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
#else
    int myproc = 0;
#endif
    if ( myproc == 0) {
        std::cout << "File header is corrupt." << std::endl;
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif
    exit(1);
  }

  file.close();


  // Prepend variable names for position to the reals list
  if (m_is_checkpoint) {
    m_pf_size_i += 2;

    m_pf_variables_i.insert(m_pf_variables_i.begin(), "cpu");
    m_pf_variables_i.insert(m_pf_variables_i.begin(), "id");
  }

  for (int i = 0; i < m_pf_size_i; ++i) {
    m_variable_map_i[m_pf_variables_i[i]] = i;
  }

  // Prepend variable names for position to the reals list
  m_pf_size_r = m_pf_dim + m_pf_size_r_extra;

  for (int i = 0; i < m_pf_dim; ++i) {
    std::stringstream ss;
    ss << "pos" << m_directions[i];
    m_pf_variables_r.insert(m_pf_variables_r.begin()+i, ss.str());
  }

  for (int i = 0; i < m_pf_size_r; ++i) {
    m_variable_map_r[m_pf_variables_r[i]] = i;
  }

  m_pf_size = m_pf_size_r + m_pf_size_i;

  int nprocs = ParallelDescriptor::NProcs();
  int myproc = ParallelDescriptor::MyProc();

  // cut the list of grids based on the number of MPI tasks
  int navg = m_file_nums[m_lev].size()/nprocs;
  int nleft = m_file_nums[m_lev].size() - navg * nprocs;

  m_ibegin = (myproc < nleft) ? myproc*(navg+1) : myproc*navg + nleft;
  m_iend   = (myproc < nleft) ? m_ibegin+(navg+1) : m_ibegin + navg;

  // Sum up the number of particles on the grids that this
  // processor mangages.
  for (int grid = m_ibegin; grid < m_iend; ++grid) {
    m_np_all_grids += m_particle_counts[m_lev][grid];
  }

  m_h_data_r.resize(m_np_all_grids);
  m_d_data_r.resize(m_np_all_grids);

#ifdef BL_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

void
pm_particles::
show_map ( )
{
  Print() << "\nVariable count: " << m_pf_size_r << "\n";
  for (int lc(0); lc<m_pf_size_r; ++lc) {
    Print() << "  Variable: " << std::setw(12) << std::setfill(' ') << std::left << m_pf_variables_r[lc]
            << "  Mapped index: " <<  m_variable_map_r[m_pf_variables_r[lc]] << "\n";
  }
  Print() << "\n";
}


void
pm_particles::
get_pf_data_r ( int const a_index )
{
  //amrex::Print() << " Getting particle data: " << m_pf_variables_r[a_index];

  // Do nothing. This variable is already in memory.
  if (a_index == m_last_var_read) return;

  std::size_t const found = m_pf_version.find("single");
  int const single_precision = (found!=std::string::npos) ? 1 : 0;

  std::size_t const rsize = single_precision ? sizeof(float) : sizeof(double);

  // Size of int and real data for a single particle.
  int const idata_size = m_pf_size_i*sizeof(int);
  int const rdata_size = m_pf_size_r*rsize;

  // Total size of data for a single particle.
  int const pdata_size = rdata_size + idata_size;

  int sum_np = 0;
  for (int grid = m_ibegin; grid < m_iend; ++grid) {

    int const np = m_particle_counts[m_lev][grid];

    if (np == 0) continue;

    //Print() << "Particle count: " << np << "\n";

    std::string read_file = getDataFileName(m_pf_name, m_lev, m_file_nums[m_lev][grid]);

    // Total size of particle data in hile.
    size_t buffer_size = pdata_size * np;

    std::vector<char> read_data(buffer_size);
    getDataBuffer(read_data, read_file, buffer_size, m_offsets[m_lev][grid]);

    // data is stored with all the ints for all the particles first, then all the reals
    // we convert this to array-of-structs for sorting.

    std::vector<char> data(rdata_size);
    {
      char* read_ptr = read_data.data();

      for (int p = 0; p < np; ++p) {

        double val_r;

        char* data_ptr = data.data();

        // Copy real particle data from buffer into data_ptr
        read_ptr = read_data.data() + idata_size*np + rdata_size*p;
        std::memcpy(data_ptr, read_ptr, rdata_size);

        // position the data pointer to the variable we want
        data_ptr = data.data() + a_index*rsize;
        std::memcpy(&val_r, data_ptr, rsize);

        m_h_data_r[p + sum_np] = val_r;

        //if (p <= -1) { Print() << " Particle: " << p << "  " << val_r << "\n"; }

      } // particles loop

      sum_np += np;

    } // scope block

  } // grid loop

  Gpu::copyAsync(Gpu::hostToDevice, m_h_data_r.begin(),
    m_h_data_r.end(), m_d_data_r.begin());

  Gpu::synchronize();
#ifdef BL_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  m_last_var_read = a_index;
}


amrex::Real
pm_particles::
sum ( std::string const a_var,
      bool const a_local )
{
  amrex::Print() << " computing sum(" << a_var << ")";
  int const var_index = index(a_var);

  get_pf_data_r(var_index);

  Real* p_data = m_d_data_r.data();

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  reduce_op.eval(m_np_all_grids, reduce_data, [p_data]
  AMREX_GPU_DEVICE (int p) -> ReduceTuple
  { return p_data[p]; });

  Real sm = amrex::get<0>(reduce_data.value(reduce_op));

  if (!a_local) { ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub()); }

  amrex::Print() << " done.\n";
  return sm;

}


amrex::Real
pm_particles::
fluct ( std::string const a_var,
        Real const a_avg,
        bool const a_local )
{
  amrex::Print() << " computing fluct(" << a_var << ")";
  int const var_index = index(a_var);

  get_pf_data_r(var_index);

  Real* p_data = m_d_data_r.data();

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  reduce_op.eval(m_np_all_grids, reduce_data, [p_data, a_avg]
  AMREX_GPU_DEVICE (int p) -> ReduceTuple
  { Real t = p_data[p] - a_avg;
    return t*t;
  });

  Real sm = amrex::get<0>(reduce_data.value(reduce_op));

  if (!a_local) { ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub()); }

  amrex::Print() << " done.\n";
  return sm;

}


std::istream& operator>> (std::istream& a_stream, pm_particles& a_header) {

  a_stream >> a_header.m_pf_version;
  a_stream >> a_header.m_pf_dim;

  a_stream >> a_header.m_pf_size_r_extra;
  for (int i = 0; i < a_header.m_pf_size_r_extra; ++i) {
    std::string comp;
    a_stream >> comp;
    a_header.m_pf_variables_r.push_back(comp);
  }

  a_stream >> a_header.m_pf_size_i;
  for (int i = 0; i < a_header.m_pf_size_i; ++i) {
    std::string comp;
    a_stream >> comp;
    a_header.m_pf_variables_i.push_back(comp);
  }

  a_stream >> a_header.m_is_checkpoint;
  a_stream >> a_header.m_nparticles;
  a_stream >> a_header.m_next_id;
  a_stream >> a_header.m_finest_level;

  for (int i = 0; i <= a_header.m_finest_level; ++i) {
    int grid_count;
    a_stream >> grid_count;
    a_header.m_size_grids.push_back(grid_count);
  }

  a_header.m_file_nums.resize(a_header.m_finest_level+1);
  a_header.m_particle_counts.resize(a_header.m_finest_level+1);
  a_header.m_offsets.resize(a_header.m_finest_level+1);

  for (int i = 0; i <= a_header.m_finest_level; ++i) {
    for (int j = 0; j < a_header.m_size_grids[i]; ++j) {
      int num;
      a_stream >> num;
      a_header.m_file_nums[i].push_back(num);
      a_stream >> num;
      a_header.m_particle_counts[i].push_back(num);
      a_stream >> num;
      a_header.m_offsets[i].push_back(num);
    }
  }
  return a_stream;
}


std::string
pm_particles::
getDataFileName( const std::string& prefix, int level, int file_num) {
    std::stringstream ss;
    ss << prefix << "/Level_" << level << "/DATA_";
    ss << std::setfill('0');
    ss << std::setw(5);
    ss << file_num;
    return ss.str();
}

void
pm_particles::
getDataBuffer( std::vector<char>& buffer, const std::string& file,
               size_t buffer_size, int offset) {
    std::ifstream is(file.c_str(), std::ifstream::binary);
    assert(is.is_open());
    is.seekg(offset);
    is.read(buffer.data(), buffer_size);
    is.close();
}
