#include <AMReX.H>
#include <AMReX_Print.H>

#include <post_mfix.H>

using namespace amrex;

void
post_mfix::
report( )
{

  if (m_no_work) { return; }

  // Total number of particles in plot file
  int const pcount = m_particles.particle_count();

  const MultiFab& fluid_data = m_fluid.get();

  amrex::Vector<Real> par_avgs(size_particles(), 0.);
  amrex::Vector<Real> par_flct(size_particles(), 0.);

  for (int lc = 0; lc < size_particles(); ++lc) {
    std::string var = variable_particles(lc);

    par_avgs[lc] = m_particles.sum(var);
    par_avgs[lc] /= static_cast<Real>(pcount);

    par_flct[lc] = m_particles.fluct(var, par_avgs[lc]);
    par_flct[lc] /= static_cast<Real>(pcount);

  }


  // Step 1: Compute <ep_g> = sum(ep_g) / Ncells
  // ---------------------------------------------------------------------------
  Real sum_epg(0.);

  sum_epg = fluid_data.sum(m_fluid.index("ep_g"));
  ParallelAllReduce::Sum(sum_epg, ParallelContext::CommunicatorSub());

  Real const r_ncells = static_cast<Real>(m_fluid.ncells());
  Real const avg_epg = sum_epg / r_ncells;



  // Step 2: Compute average ep_g' = sqrt( sum((ep_g(i,j,k) - <ep_g>)^2) / N )
  // ---------------------------------------------------------------------------

  int pos(0);
  amrex::Vector<Real> sums(size()+1, 0.);

  {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    auto const& xma = fluid_data.const_arrays();
    int  const epg_comp = m_fluid.index("ep_g");
    sums[pos++] = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, fluid_data, IntVect(0),
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
    {
      auto const& xfab = xma[box_no];
      Real t = xfab(i,j,k,epg_comp) - avg_epg;
      return t*t;
    });
  }

  // Step 3: Compute average ep_g' = sqrt( sum((ep_g(i,j,k) - <ep_g>)^2) / N )
  // ---------------------------------------------------------------------------

  // We set the local flag in the call to Dot to true so that no MPI calls are
  // made. We do this all at once after the sums vector is filled.

  for (int lc = 0; lc < size_fluid(); ++lc) {

    std::string const var = variable_fluid(lc);

    Print() << " computing Dot(ep_g," << var << ")";

    sums[pos++] = amrex::Dot(fluid_data, m_fluid.index("ep_g"),
                             fluid_data, m_fluid.index(var),
                             1, IntVect(0), true);
    Print() << " done.\n";
  }

  // Collect to IO proc
  ParallelDescriptor::ReduceRealSum(sums.data(), sums.size());


  // Step 3: Compute averages
  // ---------------------------------------------------------------------------

  amrex::Vector<Real> avgs(size()+1, 0.);

  // <ep_g'>
  avgs[0] = std::sqrt(sums[0] / r_ncells);

  int const offset = 1;

  // <ep_g * var> = sum(ep_g*var) / (sum ep_g)
  for (int lc = 0; lc < size(); ++lc) {
    avgs[lc+offset] = sums[lc+offset]/sum_epg;
  }


  Print() << "\n\nResults:\n";
  Print() << "==================================================\n";

  Print() << std::setw(12) << std::setfill(' ') << std::left << "<ep_g>";
  Print() << std::setw(18) << std::setprecision(12) << std::scientific << avg_epg << "\n";

  Print() << std::setw(12) << std::setfill(' ') << std::left << "<ep_g'>";
  Print() << std::setw(18) << std::setprecision(12) << std::scientific << avgs[0] << "\n";

  for (int lc = 0; lc < size_fluid(); ++lc) {
    Print() << std::setw(12) << std::setfill(' ') << std::left << "<"  + variable_fluid(lc) + ">";
    Print() << std::setw(18) << std::setprecision(12) << std::scientific << avgs[lc+offset] << "\n";
  }


  for (int lc = 0; lc < size_particles(); ++lc) {
    std::stringstream ss;
    ss << "<" << variable_particles(lc) << ">";
    Print() << std::setw(12) << std::setfill(' ') << std::left << ss.str();
    Print() << std::setw(18) << std::setprecision(12) << std::scientific << par_avgs[lc] << "\n";
  }
  for (int lc = 0; lc < size_particles(); ++lc) {
    std::stringstream ss;
    ss << "<" << variable_particles(lc) << "'>";
    Print() << std::setw(12) << std::setfill(' ') << std::left << ss.str();
    Print() << std::setw(18) << std::setprecision(12) << std::scientific << par_flct[lc] << "\n";
  }

  Print() << "==================================================\n\n\n";

}
