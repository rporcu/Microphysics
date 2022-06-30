#include <mfix_monitors.H>

#include <AMReX_FabArrayUtility.H>
#include <AMReX_Loop.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>

#include <cstdlib>
#include <sys/stat.h>

using namespace amrex;
using MFIXParIter = MFIXParticleContainer::MFIXParIter;


Monitor::Monitor (const std::array<std::string,2> specs,
                  const std::string pparse,
                  Vector<std::unique_ptr<EBFArrayBoxFactory>>& ebfactory,
                  Regions& regions)
  : m_nlev(0)
  , m_specs(specs)
  , m_monitoring_results(0)
  , m_ebfactory(ebfactory)
  , m_regions(regions)
  , m_filename(std::string())
  , m_plot_int(-1)
  , m_plot_per_approx(0.)
  , m_region_name(std::string())
  , m_region(nullptr)
  , m_variables(0)
  , m_boxes(0)
  , m_direction(-1)
  , m_plane_name(std::string())
  , m_plane(nullptr)
{
  ParmParse ppMonitor(pparse.c_str());

  ppMonitor.get("plot_file", m_filename);
  m_filename.append(".csv");
  // Remove any file having the same name to avoid appending values to an older
  // existing .csv file
  std::string cmd = "rm -f " + m_filename;
  std::system(cmd.c_str());

  ppMonitor.getarr("variables", m_variables);
  m_variables_nb = m_variables.size();

  const int plot_int = ppMonitor.query("plot_int", m_plot_int);
  const int per_approx = ppMonitor.query("plot_per_approx", m_plot_per_approx);

  // XOR operation
  AMREX_ALWAYS_ASSERT(plot_int ^ per_approx);

  ppMonitor.get("region", m_region_name);

  if (m_specs[0].compare("flowrate") == 0) {
    ppMonitor.get("plane", m_plane_name);
  }
}


void
Monitor::setup ()
{
  m_nlev = m_ebfactory.size();

  AMREX_ALWAYS_ASSERT(m_nlev > 0);

  m_monitoring_results.clear();
  m_monitoring_results.resize(m_nlev, Vector<Real>());

  AMREX_ALWAYS_ASSERT(m_regions.is_initialized());
  m_region = m_regions.get_region(m_region_name);

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_region != nullptr, "Region does not exist");

  AMREX_ALWAYS_ASSERT(m_region->ok());

  for (int lev(0); lev < m_nlev; ++lev) {

    Box box = calc_box(m_ebfactory[lev]->Geom(), *m_region);
    AMREX_ALWAYS_ASSERT(box.ok());

    m_boxes.push_back(box);
  }
}


void
Monitor::write_csv (const Real& time,
                    const Real& dt)
// NOTE: right now we only print out lev 0
{
  this->setup();

  this->monitor(dt);

  // MPI synchronization
  ParallelDescriptor::Barrier();

  // Write csv output file headers
  if (ParallelDescriptor::IOProcessor()) {

    // Check if .csv output file does not exists
    // If file does not exist, write headers
    struct stat buf;
    if (stat(m_filename.c_str(), &buf) == -1) {

      std::ofstream output_file;
      output_file.open(m_filename.c_str(), std::ios::out | std::ios::trunc);

      output_file << "time";
      for (const std::string& variable: m_variables)
        output_file << "," << variable;
      output_file << std::endl;

      const int lev = 0;
//    for (int lev(0); lev < m_nlev; ++lev) {
        output_file << time;
        for (int var(0); var < m_variables_nb; ++var) {
          output_file << "," << m_monitoring_results[lev][var];
        }
        output_file << std::endl;
//    }

      output_file.close();

    } else { // file existed, append monitored values

      std::ofstream output_file;
      output_file.open(m_filename.c_str(), std::ios::out | std::ios::app);

      const int lev = 0;
//    for (int lev(0); lev < m_nlev; ++lev) {
        output_file << time;
        for (int var(0); var < m_variables_nb; ++var) {
          output_file << "," << m_monitoring_results[lev][var];
        }
        output_file << std::endl;
//    }

      output_file.close();
    }
  }
}


Box
Monitor::calc_box (const Geometry& geometry,
                   const RealBox& realbox) const
{
  const GpuArray<Real,3> dxi = geometry.InvCellSizeArray();

  const GpuArray<Real,3> prob_lo = geometry.ProbLoArray();
  const GpuArray<Real,3> prob_hi = geometry.ProbHiArray();

  const Real* realbox_lo = realbox.lo();
  const Real* realbox_hi = realbox.hi();

  for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
    AMREX_ALWAYS_ASSERT(!(realbox_lo[dir] < prob_lo[dir]));
    AMREX_ALWAYS_ASSERT(!(realbox_hi[dir] > prob_hi[dir]));
  }

  const IntVect box_type(AMREX_D_DECL(
        static_cast<int>(std::abs(realbox_hi[0] - realbox_lo[0]) < 1.e-15),
        static_cast<int>(std::abs(realbox_hi[1] - realbox_lo[1]) < 1.e-15),
        static_cast<int>(std::abs(realbox_hi[2] - realbox_lo[2]) < 1.e-15)));

  const int sum = AMREX_D_TERM(box_type[0], + box_type[1], + box_type[2]);
  const int prod = AMREX_D_TERM(box_type[0], * box_type[1], * box_type[2]);

  AMREX_ALWAYS_ASSERT((sum <= 1) && (prod == 0));

  IntVect box_lo, box_hi;

  for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {

    if (box_type[dir] == 0) {

      box_lo[dir] = static_cast<int>(Math::ceil((realbox_lo[dir]-prob_lo[dir])*dxi[dir]));
      box_hi[dir] = static_cast<int>(Math::floor((realbox_hi[dir]-prob_lo[dir])*dxi[dir]));
    
    } else if (box_type[dir] == 1) {

      box_lo[dir] = static_cast<int>(std::round((realbox_lo[dir]-prob_lo[dir])*dxi[dir]));
      box_hi[dir] = box_lo[dir];

    } else {
      amrex::Abort("How did we arrive here?");
    }
  }

  for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
    AMREX_ALWAYS_ASSERT(box_lo[dir] <= box_hi[dir]);
  }

  return Box(box_lo, box_hi, box_type);
}


RealBox
Monitor::calc_realbox (const Geometry& geometry,
                       const Box& box) const
{
  const GpuArray<Real,3> dx = geometry.CellSizeArray();

  const GpuArray<Real,3> prob_lo = geometry.ProbLoArray();
  const GpuArray<Real,3> prob_hi = geometry.ProbHiArray();

  const IntVect box_lo = box.smallEnd();
  const IntVect box_hi = box.bigEnd();

  RealBox realbox;

  realbox.setLo({AMREX_D_DECL(prob_lo[0]+box_lo[0]*dx[0],
                              prob_lo[1]+box_lo[1]*dx[1],
                              prob_lo[2]+box_lo[2]*dx[2])});

  realbox.setHi({AMREX_D_DECL(prob_lo[0]+box_hi[0]*dx[0],
                              prob_lo[1]+box_hi[1]*dx[1],
                              prob_lo[2]+box_hi[2]*dx[2])});

  const Real* realbox_lo = realbox.lo();
  const Real* realbox_hi = realbox.hi();

  for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
    AMREX_ALWAYS_ASSERT(!(realbox_lo[dir] < prob_lo[dir]));
    AMREX_ALWAYS_ASSERT(!(realbox_hi[dir] > prob_hi[dir]));
  }

  return realbox;
}


RealBox
Monitor::realboxes_intersection (const RealBox& realbox_a,
                                 const RealBox& realbox_b) const
{
  RealBox realbox;

  const Real* a_lo = realbox_a.lo();
  const Real* b_lo = realbox_b.lo();

  realbox.setLo({AMREX_D_DECL(amrex::max(a_lo[0], b_lo[0]),
                              amrex::max(a_lo[1], b_lo[1]),
                              amrex::max(a_lo[2], b_lo[2]))});

  const Real* a_hi = realbox_a.hi();
  const Real* b_hi = realbox_b.hi();

  realbox.setHi({AMREX_D_DECL(amrex::min(a_hi[0], b_hi[0]),
                              amrex::min(a_hi[1], b_hi[1]),
                              amrex::min(a_hi[2], b_hi[2]))});

  return realbox;
}


namespace EulerianMonitor {

void
BaseMonitor::set_vars_and_comps (const int lev,
                                 Vector<const MultiFab*>& vars,
                                 Vector<int>& comps)
{
  vars.clear();
  comps.clear();

  Vector<std::string> variables_names;
  variables_names.clear();

  for (int i(0); i < m_variables_nb; ++i) {
  
    const std::string& var = m_variables[i];

    if (var.compare("ep_g") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->ep_g);
      comps.push_back(0);

    } else if (var.compare("p_g") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->p_g);
      comps.push_back(0);

    } else if (var.compare("ro_g") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->ro_g);
      comps.push_back(0);

    } else if (var.compare("trac") == 0) {

      const int ntrac = m_leveldata[lev]->trac->nComp();

      if (ntrac == 1) {
        variables_names.push_back(var);
        vars.push_back(m_leveldata[lev]->trac);
        comps.push_back(0);
      } else {
        
        for (int n(0); n < ntrac; ++n) {
          variables_names.push_back("trac_"+std::to_string(n+1));
          vars.push_back(m_leveldata[lev]->trac);
          comps.push_back(n);
        }
      }

    } else if (var.compare("vel_g") == 0) {

      variables_names.push_back("vel_g_x");
      vars.push_back(m_leveldata[lev]->vel_g);
      comps.push_back(0);

      variables_names.push_back("vel_g_y");
      vars.push_back(m_leveldata[lev]->vel_g);
      comps.push_back(1);

      variables_names.push_back("vel_g_z");
      vars.push_back(m_leveldata[lev]->vel_g);
      comps.push_back(2);

    } else if (var.compare("vel_g_x") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->vel_g);
      comps.push_back(0);

    } else if (var.compare("vel_g_y") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->vel_g);
      comps.push_back(1);

    } else if (var.compare("vel_g_z") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->vel_g);
      comps.push_back(2);

    } else if (var.compare("gp") == 0) {

      variables_names.push_back("gp_x");
      vars.push_back(m_leveldata[lev]->gp);
      comps.push_back(0);

      variables_names.push_back("gp_y");
      vars.push_back(m_leveldata[lev]->gp);
      comps.push_back(1);

      variables_names.push_back("gp_z");
      vars.push_back(m_leveldata[lev]->gp);
      comps.push_back(2);

    } else if (var.compare("gp_x") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->gp);
      comps.push_back(0);

    } else if (var.compare("gp_y") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->gp);
      comps.push_back(1);

    } else if (var.compare("gp_z") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->gp);
      comps.push_back(2);

    } else if (var.compare("T_g") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->T_g);
      comps.push_back(0);

    } else if (var.compare("h_g") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->h_g);
      comps.push_back(0);

    } else if (var.compare("X_gk") == 0) {

      for (int n_g(0); n_g < m_fluid.nspecies; ++n_g) {
        variables_names.push_back("X_gk_"+m_fluid.species_names[n_g]);
        vars.push_back(m_leveldata[lev]->X_gk);
        comps.push_back(n_g);
      }

    } else if (var.substr(0,5).compare("X_gk_") == 0) {

      const int var_name_size = var.size();
      const std::string var_species = var.substr(5,var_name_size-5);

      for (int n_g(0); n_g < m_fluid.nspecies; ++n_g) {
        if (var_species.compare(m_fluid.species_names[n_g]) == 0) {
          variables_names.push_back("X_gk_"+m_fluid.species_names[n_g]);
          vars.push_back(m_leveldata[lev]->X_gk);
          comps.push_back(n_g);
          break;
        }
      }

    } else if (var.compare("vort") == 0) {

      variables_names.push_back("vort_x");
      vars.push_back(m_leveldata[lev]->vort);
      comps.push_back(0);

      variables_names.push_back("vort_y");
      vars.push_back(m_leveldata[lev]->vort);
      comps.push_back(1);

      variables_names.push_back("vort_z");
      vars.push_back(m_leveldata[lev]->vort);
      comps.push_back(2);

    } else if (var.compare("vort_x") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->vort);
      comps.push_back(0);

    } else if (var.compare("vort_y") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->vort);
      comps.push_back(1);

    } else if (var.compare("vort_z") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->vort);
      comps.push_back(2);

    } else if (var.compare("txfr_velocity") == 0) {

      variables_names.push_back("txfr_vel_x");
      vars.push_back(m_leveldata[lev]->txfr);
      comps.push_back(Transfer::velx);

      variables_names.push_back("txfr_vel_y");
      vars.push_back(m_leveldata[lev]->txfr);
      comps.push_back(Transfer::vely);

      variables_names.push_back("txfr_vel_z");
      vars.push_back(m_leveldata[lev]->txfr);
      comps.push_back(Transfer::velz);

    } else if (var.compare("txfr_vel_x") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->txfr);
      comps.push_back(Transfer::velx);

    } else if (var.compare("txfr_vel_y") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->txfr);
      comps.push_back(Transfer::vely);

    } else if (var.compare("txfr_vel_z") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->txfr);
      comps.push_back(Transfer::velz);

    } else if (var.compare("txfr_beta") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->txfr);
      comps.push_back(Transfer::beta);

    } else if (var.compare("txfr_gammaTp") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->txfr);
      comps.push_back(Transfer::gammaTp);

    } else if (var.compare("txfr_gamma") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->txfr);
      comps.push_back(Transfer::gamma);

    } else if (var.compare("chem_txfr_X_gk") == 0) {

      const LevelData& ld = *(m_leveldata[lev]);
      ChemTransfer chem_transfer(ld.fluid->nspecies, ld.reactions->nreactions);

      for (int n_g(0); n_g < m_fluid.nspecies; ++n_g) {
        variables_names.push_back("chem_txfr_X_gk_"+m_fluid.species_names[n_g]);
        vars.push_back(m_leveldata[lev]->chem_txfr);
        comps.push_back(chem_transfer.ro_gk_txfr+n_g);
      }

    } else if (var.substr(0,15).compare("chem_txfr_X_gk_") == 0) {

      const LevelData& ld = *(m_leveldata[lev]);
      ChemTransfer chem_transfer(ld.fluid->nspecies, ld.reactions->nreactions);

      const int var_name_size = var.size();
      const std::string var_species = var.substr(15,var_name_size-15);

      for (int n_g(0); n_g < m_fluid.nspecies; ++n_g) {
        if (var_species.compare(m_fluid.species_names[n_g]) == 0) {
          variables_names.push_back("chem_txfr_X_gk_"+m_fluid.species_names[n_g]);
          vars.push_back(m_leveldata[lev]->chem_txfr);
          comps.push_back(chem_transfer.ro_gk_txfr+n_g);
          break;
        }
      }

    } else if (var.compare("chem_txfr_velocity") == 0) {

      const LevelData& ld = *(m_leveldata[lev]);
      ChemTransfer chem_transfer(ld.fluid->nspecies, ld.reactions->nreactions);

      variables_names.push_back("chem_txfr_vel_x");
      vars.push_back(m_leveldata[lev]->chem_txfr);
      comps.push_back(chem_transfer.vel_g_txfr+0);

      variables_names.push_back("chem_txfr_vel_y");
      vars.push_back(m_leveldata[lev]->chem_txfr);
      comps.push_back(chem_transfer.vel_g_txfr+1);

      variables_names.push_back("chem_txfr_vel_z");
      vars.push_back(m_leveldata[lev]->chem_txfr);
      comps.push_back(chem_transfer.vel_g_txfr+2);

    } else if (var.compare("chem_txfr_vel_x") == 0) {

      const LevelData& ld = *(m_leveldata[lev]);
      ChemTransfer chem_transfer(ld.fluid->nspecies, ld.reactions->nreactions);

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->chem_txfr);
      comps.push_back(chem_transfer.vel_g_txfr+0);

    } else if (var.compare("chem_txfr_vel_y") == 0) {

      const LevelData& ld = *(m_leveldata[lev]);
      ChemTransfer chem_transfer(ld.fluid->nspecies, ld.reactions->nreactions);

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->chem_txfr);
      comps.push_back(chem_transfer.vel_g_txfr+1);

    } else if (var.compare("chem_txfr_vel_z") == 0) {

      const LevelData& ld = *(m_leveldata[lev]);
      ChemTransfer chem_transfer(ld.fluid->nspecies, ld.reactions->nreactions);

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->chem_txfr);
      comps.push_back(chem_transfer.vel_g_txfr+2);

    } else if (var.compare("chem_txfr_h") == 0) {

      const LevelData& ld = *(m_leveldata[lev]);
      ChemTransfer chem_transfer(ld.fluid->nspecies, ld.reactions->nreactions);

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->chem_txfr);
      comps.push_back(chem_transfer.h_g_txfr);

    } else if (var.compare("diveu") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->diveu);
      comps.push_back(0);

    } else if (var.compare("mac_phi") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->mac_phi);
      comps.push_back(0);

    } else if (var.compare("divtau") == 0) {

      variables_names.push_back("divtau_x");
      vars.push_back(m_leveldata[lev]->divtau_o);
      comps.push_back(0);

      variables_names.push_back("divtau_y");
      vars.push_back(m_leveldata[lev]->divtau_o);
      comps.push_back(1);

      variables_names.push_back("divtau_z");
      vars.push_back(m_leveldata[lev]->divtau_o);
      comps.push_back(2);

    } else if (var.compare("divtau_x") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->divtau_o);
      comps.push_back(0);

    } else if (var.compare("divtau_y") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->divtau_o);
      comps.push_back(1);

    } else if (var.compare("divtau_z") == 0) {

      variables_names.push_back(var);
      vars.push_back(m_leveldata[lev]->divtau_o);
      comps.push_back(2);

    } else if (var.compare("none") == 0) {

      if (m_specs[1].compare("area") == 0) {

        variables_names.push_back("area");
        vars.push_back(nullptr);
        comps.push_back(0);

      } else if (m_specs[1].compare("massflowrate") == 0) {

        variables_names.push_back("mass_flow_rate");
        vars.push_back(nullptr);
        comps.push_back(0);

      } else if (m_specs[1].compare("volumeflowrate") == 0) {
        variables_names.push_back("volume_flow_rate");
        vars.push_back(nullptr);
        comps.push_back(0);

      } else if (m_specs[1].compare("volume") == 0) {

        variables_names.push_back("volume");
        vars.push_back(nullptr);
        comps.push_back(0);

      } else {

        amrex::Abort("inputs error");
      }

    } else {

      Print() << "var = " << var << "\n";
      amrex::Abort("Unrecognized variable to monitor. Fix inputs");
    }
  }

  m_variables.swap(variables_names);
  m_variables_nb = m_variables.size();
}


void
PointRegion::check_boxes_are_ok () const
{
  for (int lev(0); lev < m_nlev; ++lev) {

    const IndexType& idx_type = m_boxes[lev].ixType();

    if (!idx_type.cellCentered())
      amrex::Abort("Box is not cell-centered");

    for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
      if (m_boxes[lev].length(dir) != 1)
        amrex::Abort("PointRegion monitor Box is not made of one single cell");
    }
  }
}


bool
PointRegion::check_mf_is_ok (const MultiFab& mf) const
{
  const BoxArray& box_array = mf.boxArray();
  const IndexType& idx_type = box_array.ixType();

  if (!(idx_type.cellCentered() || idx_type.nodeCentered())) {
    return false;
  }

  return true;
}


void
PointRegion::monitor (const Real& /*dt*/)
{

  if (m_specs[1].compare("value") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {

      Vector<const MultiFab*> mf(0);
      Vector<int> components(0);
      set_vars_and_comps(lev, mf, components);

      m_monitoring_results[lev] = value(lev, mf, components);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
PointRegion::value (const int lev,
                    const Vector<const MultiFab*>& mf,
                    const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0.};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*frac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_value = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_value, 1);

    result[var] = l_value;
  }

  return result;
}


void
PointRegion::convert_mf_if_needed (Vector<int>& mf_ok,
                                   Vector<const MultiFab*>& mf,
                                   Vector<int>& components,
                                   Vector<int>& conversion_flags,
                                   const int lev)
{
  AMREX_ALWAYS_ASSERT(mf.size() == mf_ok.size());
  AMREX_ALWAYS_ASSERT(mf.size() == components.size());
  AMREX_ALWAYS_ASSERT(mf.size() == conversion_flags.size());

  const int var_nb = mf.size();

  BoxArray box_array = this->get_box_array(lev);
  const DistributionMapping& d_map = this->get_d_map(lev);

  for (int var(0); var < var_nb; ++var) {

    if (mf_ok[var]) {
      if (mf[var]->boxArray().ixType().nodeCentered()) {
        MultiFab* mf_cc = new MultiFab(box_array, d_map, 1, 0);
        amrex::average_node_to_cellcenter(*mf_cc, 0, *mf[var], components[var], 1, 0);
        mf[var] = mf_cc;
        components[var] = 0;
        conversion_flags[var] = 1;
      }
    }
  }
}


void
AreaMonitor::set_direction ()
{
  const IndexType& idx_type = m_boxes[0].ixType();

  for (int dir(0); dir < AMREX_SPACEDIM; ++dir)
    if (idx_type.nodeCentered(dir)) {
      m_direction = dir;
      break;
    }

  AMREX_ALWAYS_ASSERT(AMREX_D_TERM(m_boxes[0].length(m_direction) == 1, &&
                                   m_boxes[0].length((m_direction+1)%AMREX_SPACEDIM) > 1, &&
                                   m_boxes[0].length((m_direction+2)%AMREX_SPACEDIM) > 1));
}


void
AreaMonitor::check_boxes_are_ok () const
{
  for (int lev(0); lev < m_nlev; ++lev) {

    const int min_box_length = amrex::min(AMREX_D_DECL(m_boxes[lev].length(0),
                                                       m_boxes[lev].length(1),
                                                       m_boxes[lev].length(2)));

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(min_box_length == 1,
        "AreaMonitor Box is over-dimensionated");
  }
}


bool
AreaMonitor::check_mf_is_ok (const MultiFab& mf) const
{
  const BoxArray& box_array = mf.boxArray();
  const IndexType& idx_type = box_array.ixType();

  if (!idx_type.nodeCentered(m_direction))
    return false;

  for (int dir(1); dir < AMREX_SPACEDIM; ++dir)
    if (!idx_type.cellCentered((m_direction+dir)%AMREX_SPACEDIM))
      return false;

  return true;
}


void
AreaMonitor::convert_mf_if_needed (Vector<int>& mf_ok,
                                   Vector<const MultiFab*>& mf,
                                   Vector<int>& components,
                                   Vector<int>& conversion_flags,
                                   const int lev)
{
  AMREX_ALWAYS_ASSERT(mf.size() == mf_ok.size());
  AMREX_ALWAYS_ASSERT(mf.size() == components.size());
  AMREX_ALWAYS_ASSERT(mf.size() == conversion_flags.size());

  const int var_nb = mf.size();

  const DistributionMapping& d_map = this->get_d_map(lev);

  for (int var(0); var < var_nb; ++var) {

    if (mf_ok[var]) {
      if (mf[var]->boxArray().ixType().cellCentered()) {
        Array<MultiFab*, AMREX_SPACEDIM> mf_fc;

        for(int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          BoxArray edge_ba = m_ebfactory[lev]->boxArray();
          edge_ba.surroundingNodes(dir);
          mf_fc[dir] = new MultiFab(edge_ba, d_map, 1, 0);
          mf_fc[dir]->setVal(0.);
        }

        const auto& geom = m_ebfactory[lev]->Geom();
        Vector<BCRec> bc_rec(1, BCRec());

        EB_interp_CellCentroid_to_FaceCentroid(*mf[var], mf_fc, components[var], 0, 1, geom, bc_rec);

        mf[var] = mf_fc[m_direction];
        components[var] = 0;
        conversion_flags[var] = 1;

        for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
          if (dir != m_direction)
            delete mf_fc[dir];
        }
      }
    }
  }
}


void
AreaRegion::monitor (const Real& /*dt*/)
{
  // Setup the variables and components
  Vector<Vector<const MultiFab*>> mf(m_nlev, Vector<const MultiFab*>());
  Vector<Vector<int>> components(m_nlev, Vector<int>());
  Vector<Vector<int>> conversion_flags(m_nlev, Vector<int>());

  for (int lev(0); lev < m_nlev; ++lev) {
    set_vars_and_comps(lev, mf[lev], components[lev]);
  }

  // Select the monitor
  if (m_specs[1].compare("sum") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = sum(lev, mf[lev], components[lev]);
    }

  } else if (m_specs[1].compare("min") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = min(lev, mf[lev], components[lev]);
    }

  } else if (m_specs[1].compare("max") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = max(lev, mf[lev], components[lev]);
    }

  } else if (m_specs[1].compare("average") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = average(lev, mf[lev], components[lev]);
    }

  } else if (m_specs[1].compare("standarddeviation") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = stddev(lev, mf[lev], components[lev]);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
AreaRegion::sum (const int lev,
                 const Vector<const MultiFab*>& mf,
                 const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0.};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*areafrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_sum = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_sum, 1);

    result[var] = l_sum;
  }

  return result;
}


Vector<Real>
AreaRegion::min (const int lev,
                 const Vector<const MultiFab*>& mf,
                 const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {std::numeric_limits<Real>::max()};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*areafrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_min = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMin(&l_min, 1);

    result[var] = l_min;
  }

  return result;
}


Vector<Real>
AreaRegion::max (const int lev,
                 const Vector<const MultiFab*>& mf,
                 const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {std::numeric_limits<Real>::min()};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*areafrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_max = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMax(&l_max, 1);

    result[var] = l_max;
  }

  return result;
}


Vector<Real>
AreaRegion::average (const int lev,
                     const Vector<const MultiFab*>& mf,
                     const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,long> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_values = {0., 0};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*areafrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return {mf_arr(ijk), 1}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_values);

    Real numerator = amrex::get<0>(host_tuple);
    long denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceLongSum(&denominator, 1);
  
    result[var] = numerator / Real(denominator);
  }

  return result;
}


Vector<Real>
AreaRegion::stddev (const int lev,
                    const Vector<const MultiFab*>& mf,
                    const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  const Vector<Real> avg = this->average(lev, mf, components);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,long> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_values = {0., 0};

    const Real d_avg = avg[var];

    auto R = [d_avg] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                         Array4<Real const> const& /*areafrac_arr*/,
                                         Array4<Real const> const& /*epsilon_arr*/,
                                         Array4<Real const> const& /*density_arr*/,
                                         Array4<Real const> const& /*velocity_arr*/,
                                         const IntVect& ijk) -> ReduceTuple
    {
      Real diff = mf_arr(ijk) - d_avg;
      return {diff*diff, 1};
    };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_values);

    Real numerator = amrex::get<0>(host_tuple);
    long denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceLongSum(&denominator, 1);

    result[var] = std::sqrt(numerator / Real(denominator));
  }

  return result;
}


void
VolumeMonitor::check_boxes_are_ok () const
{
  for (int lev(0); lev < m_nlev; ++lev) {

    const IndexType& idx_type = m_boxes[lev].ixType();

    if (!idx_type.cellCentered())
      amrex::Abort("Box is not cell-centered");

    for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
      if (m_boxes[lev].length(dir) < 1)
        amrex::Abort("VolumeMonitor monitor Box is not made of one single cell");
    }
  }
}


bool
VolumeMonitor::check_mf_is_ok (const MultiFab& mf) const
{
  const BoxArray& box_array = mf.boxArray();
  const IndexType& idx_type = box_array.ixType();

  if (!(idx_type.cellCentered() || idx_type.nodeCentered()))
    return false;

  return true;
}


void
VolumeMonitor::convert_mf_if_needed (Vector<int>& mf_ok,
                                     Vector<const MultiFab*>& mf,
                                     Vector<int>& components,
                                     Vector<int>& conversion_flags,
                                     const int lev)
{
  AMREX_ALWAYS_ASSERT(mf.size() == mf_ok.size());
  AMREX_ALWAYS_ASSERT(mf.size() == components.size());
  AMREX_ALWAYS_ASSERT(mf.size() == conversion_flags.size());

  const int var_nb = mf.size();

  BoxArray box_array = this->get_box_array(lev);
  const DistributionMapping& d_map = this->get_d_map(lev);

  for (int var(0); var < var_nb; ++var) {

    if (mf_ok[var]) {
      if (mf[var]->boxArray().ixType().nodeCentered()) {
        MultiFab* mf_cc = new MultiFab(box_array, d_map, 1, 0);
        amrex::average_node_to_cellcenter(*mf_cc, 0, *mf[var], components[var], 1, 0);
        mf[var] = mf_cc;
        components[var] = 0;
        conversion_flags[var] = 1;
      }
    }
  }
}


void
VolumeRegion::monitor (const Real& /*dt*/)
{
  // Setup the variables and components
  Vector<Vector<const MultiFab*>> mf(m_nlev, Vector<const MultiFab*>());
  Vector<Vector<int>> components(m_nlev, Vector<int>());
  Vector<Vector<int>> conversion_flags(m_nlev, Vector<int>());

  for (int lev(0); lev < m_nlev; ++lev) {
    set_vars_and_comps(lev, mf[lev], components[lev]);
  }

  // Select the monitor
  if (m_specs[1].compare("sum") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = sum(lev, mf[lev], components[lev]);
    }

  } else if (m_specs[1].compare("min") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = min(lev, mf[lev], components[lev]);
    }

  } else if (m_specs[1].compare("max") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = max(lev, mf[lev], components[lev]);
    }

  } else if (m_specs[1].compare("average") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = average(lev, mf[lev], components[lev]);
    }

  } else if (m_specs[1].compare("standarddeviation") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = stddev(lev, mf[lev], components[lev]);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
VolumeRegion::sum (const int lev,
                   const Vector<const MultiFab*>& mf,
                   const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0.};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*volfrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_sum = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_sum, 1);

    result[var] = l_sum;
  }

  return result;
}


Vector<Real>
VolumeRegion::min (const int lev,
                   const Vector<const MultiFab*>& mf,
                   const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {std::numeric_limits<Real>::max()};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*volfrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_min = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMin(&l_min, 1);

    result[var] = l_min;
  }

  return result;
}


Vector<Real>
VolumeRegion::max (const int lev,
                   const Vector<const MultiFab*>& mf,
                   const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {std::numeric_limits<Real>::min()};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*volfrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_max = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMax(&l_max, 1);

    result[var] = l_max;
  }

  return result;
}


Vector<Real>
VolumeRegion::average (const int lev,
                       const Vector<const MultiFab*>& mf,
                       const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,long> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_values = {0., 0};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*volfrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return {mf_arr(ijk), 1}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_values);

    Real numerator = amrex::get<0>(host_tuple);
    long denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceLongSum(&denominator, 1);

    result[var] = numerator / Real(denominator);
  }

  return result;
}


Vector<Real>
VolumeRegion::stddev (const int lev,
                      const Vector<const MultiFab*>& mf,
                      const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  const Vector<Real> avg = this->average(lev, mf, components);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,long> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_values = {0., 0};

    const Real d_avg = avg[var];

    auto R = [d_avg] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                       Array4<Real const> const& /*volfrac_arr*/,
                                       Array4<Real const> const& /*epsilon_arr*/,
                                       Array4<Real const> const& /*density_arr*/,
                                       Array4<Real const> const& /*velocity_arr*/,
                                       const IntVect& ijk) -> ReduceTuple
    {
      Real diff = mf_arr(ijk) - d_avg;
      return {diff*diff, 1};
    };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_values);

    Real numerator = amrex::get<0>(host_tuple);
    long denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceLongSum(&denominator, 1);

    result[var] = std::sqrt(numerator / Real(denominator));
  }

  return result;
}


void
SurfaceIntegral::set_dA ()
{
  for (int lev(0); lev < m_nlev; ++lev) {

    const Real* dx = m_ebfactory[lev]->Geom().CellSize();

    m_dA[lev] = 1.;

    for (int dir(1); dir < AMREX_SPACEDIM; ++dir) {
      const int i = (m_direction + dir) % AMREX_SPACEDIM;
      m_dA[lev] *= dx[i];
    }
  }
}


void
SurfaceIntegral::monitor (const Real& /*dt*/)
{
  // Setup the variables and components
  Vector<Vector<const MultiFab*>> mf(m_nlev, Vector<const MultiFab*>());
  Vector<Vector<int>> components(m_nlev, Vector<int>());

  for (int lev(0); lev < m_nlev; ++lev) {
    set_vars_and_comps(lev, mf[lev], components[lev]);
  }

  // Select the monitor
  if (m_specs[1].compare("area") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = area(lev);
    }

  } else if (m_specs[1].compare("areaweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = area_weighted_average(lev, mf[lev], components[lev]);
    }

  } else if (m_specs[1].compare("flowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      const MultiFab* epsilon = m_leveldata[lev]->ep_g;
      const MultiFab* density = m_leveldata[lev]->ro_g;
      const MultiFab* velocity = m_leveldata[lev]->vel_g;

      m_monitoring_results[lev] = flow_rate(lev, mf[lev], components[lev], epsilon, density, velocity);
    }

  } else if (m_specs[1].compare("massflowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      const MultiFab* epsilon = m_leveldata[lev]->ep_g;
      const MultiFab* density = m_leveldata[lev]->ro_g;
      const MultiFab* velocity = m_leveldata[lev]->vel_g;

      m_monitoring_results[lev] = mass_flow_rate(lev, epsilon, density, velocity);
    }

  } else if (m_specs[1].compare("massweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      const MultiFab* epsilon = m_leveldata[lev]->ep_g;
      const MultiFab* density = m_leveldata[lev]->ro_g;
      const MultiFab* velocity = m_leveldata[lev]->vel_g;

      m_monitoring_results[lev] = mass_weighted_average(lev, mf[lev], components[lev], epsilon, density, velocity);
    }

  } else if (m_specs[1].compare("volumeflowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      const MultiFab* epsilon = m_leveldata[lev]->ep_g;
      const MultiFab* velocity = m_leveldata[lev]->vel_g;

      m_monitoring_results[lev] = volume_flow_rate(lev, epsilon, velocity);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
SurfaceIntegral::area (const int lev)
{
  Vector<Real> result(1, 0.);

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);

  using ReduceTuple = typename decltype(reduce_data)::Type;

  ReduceTuple default_value = {0.};

  const Real dA = m_dA[lev];

  auto R = [dA] AMREX_GPU_DEVICE (Array4<Real const> const& /*mf_arr*/,
                                  Array4<Real const> const& areafrac_arr,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
  { return areafrac_arr(ijk)*dA; };

  const MultiFab* dummy_mf = nullptr;
  int dummy_component = 0;

  ReduceTuple host_tuple = apply(lev, dummy_mf, dummy_component, reduce_data, reduce_op, R, default_value);

  Real l_area = amrex::get<0>(host_tuple);

  ParallelDescriptor::ReduceRealSum(&l_area, 1);

  result[0] = l_area;

  return result;
}


Vector<Real>
SurfaceIntegral::area_weighted_average (const int lev,
                                        const Vector<const MultiFab*>& mf,
                                        const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_values = {0., 0.};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& areafrac_arr,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return {mf_arr(ijk)*areafrac_arr(ijk), areafrac_arr(ijk)}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_values);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    result[var] = numerator / denominator;
  }

  return result;
}


Vector<Real>
SurfaceIntegral::flow_rate (const int lev,
                            const Vector<const MultiFab*>& mf,
                            const Vector<int>& components,
                            const MultiFab* epsilon,
                            const MultiFab* density,
                            const MultiFab* velocity)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_values = {0.};

    const Real dA = m_dA[lev];

    auto R = [dA] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                    Array4<Real const> const& areafrac_arr,
                                    Array4<Real const> const& epsilon_arr,
                                    Array4<Real const> const& density_arr,
                                    Array4<Real const> const& velocity_arr,
                                    const IntVect& ijk) -> ReduceTuple
    { return {epsilon_arr(ijk)*density_arr(ijk)*mf_arr(ijk)*velocity_arr(ijk)*areafrac_arr(ijk)*dA}; };

    std::map<std::string, const MultiFab*> aux_mfs;
    aux_mfs["ep_g"] = epsilon;
    aux_mfs["ro_g"] = density;
    aux_mfs["vel_g"] = velocity;

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R,
                                   default_values, aux_mfs);

    Real l_flow_rate = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_flow_rate, 1);

    result[var] = l_flow_rate;
  }

  return result;
}


Vector<Real>
SurfaceIntegral::mass_flow_rate (const int lev,
                                 const MultiFab* epsilon,
                                 const MultiFab* density,
                                 const MultiFab* velocity)
{
  Vector<Real> result(1, 0.);

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);

  using ReduceTuple = typename decltype(reduce_data)::Type;

  ReduceTuple default_values = {0.};

  const Real dA = m_dA[lev];

  auto R = [dA] AMREX_GPU_DEVICE (Array4<Real const> const& /*mf_arr*/,
                                  Array4<Real const> const& areafrac_arr,
                                  Array4<Real const> const& epsilon_arr,
                                  Array4<Real const> const& density_arr,
                                  Array4<Real const> const& velocity_arr,
                                  const IntVect& ijk) -> ReduceTuple
  { return {epsilon_arr(ijk)*density_arr(ijk)*velocity_arr(ijk)*areafrac_arr(ijk)*dA}; };

  std::map<std::string, const MultiFab*> aux_mfs;
  aux_mfs["ep_g"] = epsilon;
  aux_mfs["ro_g"] = density;
  aux_mfs["vel_g"] = velocity;

  const MultiFab* dummy_mf = nullptr;
  int dummy_component = 0;

  ReduceTuple host_tuple = apply(lev, dummy_mf, dummy_component, reduce_data,
                                 reduce_op, R, default_values, aux_mfs);

  Real l_mass_flow_rate = amrex::get<0>(host_tuple);

  ParallelDescriptor::ReduceRealSum(&l_mass_flow_rate, 1);

  result[0] = l_mass_flow_rate;

  return result;
}


Vector<Real>
SurfaceIntegral::mass_weighted_average (const int lev,
                                        const Vector<const MultiFab*>& mf,
                                        const Vector<int>& components,
                                        const MultiFab* epsilon,
                                        const MultiFab* density,
                                        const MultiFab* velocity)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_values = {0., 0.};

    const Real dA = m_dA[lev];

    auto R = [dA] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                    Array4<Real const> const& areafrac_arr,
                                    Array4<Real const> const& epsilon_arr,
                                    Array4<Real const> const& density_arr,
                                    Array4<Real const> const& velocity_arr,
                                    const IntVect& ijk) -> ReduceTuple
    { return {epsilon_arr(ijk)*density_arr(ijk)*mf_arr(ijk)*velocity_arr(ijk)*areafrac_arr(ijk)*dA,
              epsilon_arr(ijk)*density_arr(ijk)*velocity_arr(ijk)*areafrac_arr(ijk)*dA}; };

    std::map<std::string, const MultiFab*> aux_mfs;
    aux_mfs["ep_g"] = epsilon;
    aux_mfs["ro_g"] = density;
    aux_mfs["vel_g"] = velocity;

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R,
                                   default_values, aux_mfs);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    result[var] = numerator / denominator;
  }

  return result;
}


Vector<Real>
SurfaceIntegral::volume_flow_rate (const int lev,
                                   const MultiFab* epsilon,
                                   const MultiFab* velocity)
{
  Vector<Real> result(1, 0.);

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);

  using ReduceTuple = typename decltype(reduce_data)::Type;

  ReduceTuple default_values = {0.};

  const Real dA = m_dA[lev];

  auto R = [dA] AMREX_GPU_DEVICE (Array4<Real const> const& /*mf_arr*/,
                                  Array4<Real const> const& areafrac_arr,
                                  Array4<Real const> const& epsilon_arr,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& velocity_arr,
                                  const IntVect& ijk) -> ReduceTuple
  { return {epsilon_arr(ijk)*velocity_arr(ijk)*areafrac_arr(ijk)*dA}; };

  std::map<std::string, const MultiFab*> aux_mfs;
  aux_mfs["ep_g"] = epsilon;
  aux_mfs["vel_g"] = velocity;

  const MultiFab* dummy_mf = nullptr;
  int dummy_component = 0;

  ReduceTuple host_tuple = apply(lev, dummy_mf, dummy_component, reduce_data,
                                 reduce_op, R, default_values, aux_mfs);

  Real l_volume_flow_rate = amrex::get<0>(host_tuple);

  ParallelDescriptor::ReduceRealSum(&l_volume_flow_rate, 1);

  result[0] = l_volume_flow_rate;

  return result;
}


void
VolumeIntegral::set_dV ()
{
  for (int lev(0); lev < m_nlev; ++lev) {

    m_dV[lev] = 1.;

    for (int dir(0); dir < AMREX_SPACEDIM; ++dir)
      m_dV[lev] *= m_ebfactory[lev]->Geom().CellSize(dir);
  }
}


void
VolumeIntegral::monitor (const Real& /*dt*/)
{
  // Setup the variables and components
  Vector<Vector<const MultiFab*>> mf(m_nlev, Vector<const MultiFab*>());
  Vector<Vector<int>> components(m_nlev, Vector<int>());

  for (int lev(0); lev < m_nlev; ++lev) {
    set_vars_and_comps(lev, mf[lev], components[lev]);
  }

  // Select the monitor
  if (m_specs[1].compare("volume") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {

      m_monitoring_results[lev] = volume(lev);
    }

  } else if (m_specs[1].compare("volumeintegral") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = volume_integral(lev, mf[lev], components[lev]);
    }

  } else if (m_specs[1].compare("volumeweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = volume_weighted_average(lev, mf[lev], components[lev]);
    }

  } else if (m_specs[1].compare("massweightedintegral") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      const MultiFab* epsilon = m_leveldata[lev]->ep_g;
      const MultiFab* density = m_leveldata[lev]->ro_g;

      m_monitoring_results[lev] = mass_weighted_integral(lev, mf[lev], components[lev], epsilon, density);
    }

  } else if (m_specs[1].compare("massweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      const MultiFab* epsilon = m_leveldata[lev]->ep_g;
      const MultiFab* density = m_leveldata[lev]->ro_g;

      m_monitoring_results[lev] = mass_weighted_average(lev, mf[lev], components[lev], epsilon, density);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
VolumeIntegral::volume (const int lev)
{
  Vector<Real> result(1, 0.);

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);

  using ReduceTuple = typename decltype(reduce_data)::Type;

  ReduceTuple default_value = {0.};

  const Real dV = m_dV[lev];

  auto R = [dV] AMREX_GPU_DEVICE (Array4<Real const> const& /*mf_arr*/,
                                  Array4<Real const> const& volfrac_arr,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
  { return volfrac_arr(ijk)*dV; };

  const MultiFab* dummy_mf = nullptr;
  int dummy_component = 0;

  ReduceTuple host_tuple = apply(lev, dummy_mf, dummy_component, reduce_data,
                                 reduce_op, R, default_value);

  Real l_volume = amrex::get<0>(host_tuple);

  ParallelDescriptor::ReduceRealSum(&l_volume, 1);

  result[0] = l_volume;

  return result;
}


Vector<Real>
VolumeIntegral::volume_integral (const int lev,
                                 const Vector<const MultiFab*>& mf,
                                 const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_values = {0.};

    const Real dV = m_dV[lev];

    auto R = [dV] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                    Array4<Real const> const& volfrac_arr,
                                    Array4<Real const> const& /*epsilon_arr*/,
                                    Array4<Real const> const& /*density_arr*/,
                                    Array4<Real const> const& /*velocity_arr*/,
                                    const IntVect& ijk) -> ReduceTuple
    { return {mf_arr(ijk)*volfrac_arr(ijk)*dV}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data,
                                   reduce_op, R, default_values);

    Real l_volume_integral = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_volume_integral, 1);

    result[var] = l_volume_integral;
  }

  return result;
}


Vector<Real>
VolumeIntegral::volume_weighted_average (const int lev,
                                         const Vector<const MultiFab*>& mf,
                                         const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_values = {0., 0.};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& volfrac_arr,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return {mf_arr(ijk)*volfrac_arr(ijk),
              volfrac_arr(ijk)}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data,
                                   reduce_op, R, default_values);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    result[var] = numerator / denominator;
  }

  return result;
}


Vector<Real>
VolumeIntegral::mass_weighted_integral (const int lev,
                                        const Vector<const MultiFab*>& mf,
                                        const Vector<int>& components,
                                        const MultiFab* epsilon,
                                        const MultiFab* density)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_values = {0.};

    const Real dV = m_dV[lev];

    auto R = [dV] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                    Array4<Real const> const& volfrac_arr,
                                    Array4<Real const> const& epsilon_arr,
                                    Array4<Real const> const& density_arr,
                                    Array4<Real const> const& /*velocity_arr*/,
                                    const IntVect& ijk) -> ReduceTuple
    { return epsilon_arr(ijk)*density_arr(ijk)*mf_arr(ijk)*volfrac_arr(ijk)*dV; };

    std::map<std::string, const MultiFab*> aux_mfs;
    aux_mfs["ep_g"] = epsilon;
    aux_mfs["ro_g"] = density;

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data,
                                   reduce_op, R, default_values, aux_mfs);

    Real l_integral = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_integral, 1);

    result[var] = l_integral;
  }

  return result;
}


Vector<Real>
VolumeIntegral::mass_weighted_average (const int lev,
                                       const Vector<const MultiFab*>& mf,
                                       const Vector<int>& components,
                                       const MultiFab* epsilon,
                                       const MultiFab* density)
{
  const int var_nb = mf.size();

  AMREX_ALWAYS_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_values = {0., 0.};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& volfrac_arr,
                                  Array4<Real const> const& epsilon_arr,
                                  Array4<Real const> const& density_arr,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return {epsilon_arr(ijk)*density_arr(ijk)*mf_arr(ijk)*volfrac_arr(ijk),
              epsilon_arr(ijk)*density_arr(ijk)*volfrac_arr(ijk)}; };

    std::map<std::string, const MultiFab*> aux_mfs;
    aux_mfs["ep_g"] = epsilon;
    aux_mfs["ro_g"] = density;

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data,
                                   reduce_op, R, default_values, aux_mfs);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    result[var] = numerator / denominator;
  }

  return result;
}

} // end namespace EulerianMonitor


namespace LagrangianMonitor {

void
BaseMonitor::set_indexes (Vector<int>& indexes)
{
  indexes.clear();

  Vector<std::string> variables_names;
  variables_names.clear();

  for (const std::string& var: m_variables) {

    if (var.compare("position") == 0) {

      variables_names.push_back("pos_x");
      indexes.push_back(0);

      variables_names.push_back("pos_y");
      indexes.push_back(1);

      variables_names.push_back("pos_z");
      indexes.push_back(2);

    } else if (var.compare("pos_x") == 0) {

      variables_names.push_back(var);
      indexes.push_back(0);

    } else if (var.compare("pos_y") == 0) {

      variables_names.push_back(var);
      indexes.push_back(1);

    } else if (var.compare("pos_z") == 0) {

      variables_names.push_back(var);
      indexes.push_back(2);

    } else if (var.compare("id") == 0) {

      variables_names.push_back(var);
      indexes.push_back(3);

    } else if (var.compare("cpu") == 0) {

      variables_names.push_back(var);
      indexes.push_back(4);

    } else if (var.compare("radius") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::radius);

    } else if (var.compare("volume") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::volume);

    } else if (var.compare("mass") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::mass);

    } else if (var.compare("density") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::density);

    } else if (var.compare("oneOverI") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::oneOverI);

    } else if (var.compare("velocity") == 0) {

      variables_names.push_back("vel_x");
      indexes.push_back(5+SoArealData::velx);

      variables_names.push_back("vel_y");
      indexes.push_back(5+SoArealData::vely);

      variables_names.push_back("vel_z");
      indexes.push_back(5+SoArealData::velz);

    } else if (var.compare("vel_x") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::velx);

    } else if (var.compare("vel_y") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::vely);

    } else if (var.compare("vel_z") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::velz);

    } else if (var.compare("omega") == 0) {

      variables_names.push_back("omega_x");
      indexes.push_back(5+SoArealData::omegax);

      variables_names.push_back("omega_y");
      indexes.push_back(5+SoArealData::omegay);

      variables_names.push_back("omega_z");
      indexes.push_back(5+SoArealData::omegaz);

    } else if (var.compare("omega_x") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::omegax);

    } else if (var.compare("omega_y") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::omegay);

    } else if (var.compare("omega_z") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::omegaz);

    } else if (var.compare("statwt") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::statwt);

    } else if (var.compare("dragcoeff") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::dragcoeff);

    } else if (var.compare("drag") == 0) {

      variables_names.push_back("drag_x");
      indexes.push_back(5+SoArealData::dragx);

      variables_names.push_back("drag_y");
      indexes.push_back(5+SoArealData::dragy);

      variables_names.push_back("drag_z");
      indexes.push_back(5+SoArealData::dragz);

    } else if (var.compare("drag_x") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::dragx);

    } else if (var.compare("drag_y") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::dragy);

    } else if (var.compare("drag_z") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::dragz);

    } else if (var.compare("cp_s") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::cp_s);

    } else if (var.compare("T_s") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::temperature);

    } else if (var.compare("convection") == 0) {

      variables_names.push_back(var);
      indexes.push_back(5+SoArealData::convection);

    } else if (var.compare("phase") == 0) {

      variables_names.push_back(var);
      const int idx = 5+SoArealData::count;
      indexes.push_back(idx+SoAintData::phase);

    } else if (var.compare("state") == 0) {

      variables_names.push_back(var);
      const int idx = 5+SoArealData::count;
      indexes.push_back(idx+SoAintData::state);

    } else if (var.compare("X_sn") == 0) {

      for (int n_s(0); n_s < m_solids.nspecies; ++n_s) {
        variables_names.push_back("X_sn_"+m_solids.species_names[n_s]);
        const int idx = 5+SoArealData::count+SoAintData::count;
        indexes.push_back(idx+n_s);
      }

    } else if (var.substr(0,5).compare("X_sn_") == 0) {

      const int var_name_size = var.size();
      const std::string var_species = var.substr(5,var_name_size-5);

      for (int n_s(0); n_s < m_solids.nspecies; ++n_s) {
        if (var_species.compare(m_solids.species_names[n_s]) == 0) {
          variables_names.push_back("X_sn_"+m_solids.species_names[n_s]);
          const int idx = 5+SoArealData::count+SoAintData::count;
          indexes.push_back(idx+n_s);
          break;
        }
      }

    } else if (var.compare("txfr_velocity") == 0) {

      const runtimeRealData& rrData = m_pc->m_runtimeRealData;

      {
        variables_names.push_back("txfr_vel_x");
        const int idx = 5+SoArealData::count+SoAintData::count;
        indexes.push_back(idx+rrData.vel_txfr+0);
      }

      {
        variables_names.push_back("txfr_vel_y");
        const int idx = 5+SoArealData::count+SoAintData::count;
        indexes.push_back(idx+rrData.vel_txfr+1);
      }

      {
        variables_names.push_back("txfr_vel_z");
        const int idx = 5+SoArealData::count+SoAintData::count;
        indexes.push_back(idx+rrData.vel_txfr+2);
      }

    } else if (var.compare("txfr_vel_x") == 0) {

      variables_names.push_back(var);
      const runtimeRealData& rrData = m_pc->m_runtimeRealData;
      const int idx = 5+SoArealData::count+SoAintData::count;
      indexes.push_back(idx+rrData.vel_txfr+0);

    } else if (var.compare("txfr_vel_y") == 0) {

      variables_names.push_back(var);
      const runtimeRealData& rrData = m_pc->m_runtimeRealData;
      const int idx = 5+SoArealData::count+SoAintData::count;
      indexes.push_back(idx+rrData.vel_txfr+1);

    } else if (var.compare("txfr_vel_z") == 0) {

      variables_names.push_back(var);
      const runtimeRealData& rrData = m_pc->m_runtimeRealData;
      const int idx = 5+SoArealData::count+SoAintData::count;
      indexes.push_back(idx+rrData.vel_txfr+2);

    } else if (var.compare("txfr_h") == 0) {

      variables_names.push_back(var);
      const runtimeRealData& rrData = m_pc->m_runtimeRealData;
      const int idx = 5+SoArealData::count+SoAintData::count;
      indexes.push_back(idx+rrData.h_txfr);

    } else if (var.compare("txfr_X_sn") == 0) {

      const runtimeRealData& rrData = m_pc->m_runtimeRealData;

      for (int n_s(0); n_s < m_solids.nspecies; ++n_s) {
        variables_names.push_back("txfr_X_sn_"+m_solids.species_names[n_s]);
        const int idx = 5+SoArealData::count+SoAintData::count;
        indexes.push_back(idx+rrData.mass_txfr+n_s);
      }

    } else if (var.substr(0,10).compare("txfr_X_sn_") == 0) {

      const runtimeRealData& rrData = m_pc->m_runtimeRealData;

      const int var_name_size = var.size();
      const std::string var_species = var.substr(10,var_name_size-10);

      for (int n_s(0); n_s < m_solids.nspecies; ++n_s) {
        if (var_species.compare(m_solids.species_names[n_s]) == 0) {
          variables_names.push_back("txfr_X_sn_"+m_solids.species_names[n_s]);
          const int idx = 5+SoArealData::count+SoAintData::count;
          indexes.push_back(idx+rrData.mass_txfr+n_s);
          break;
        }
      }

    } else {

      amrex::Abort("Unrecognized variable to monitor. Fix inputs");
    }
  }

  m_variables.swap(variables_names);
  m_variables_nb = m_variables.size();
}


void
BaseMonitor::check_realbox_is_ok () const
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_region->volume() > 0., "RealBox has non-positive volume");
}


void
GeneralProperty::monitor (const Real& /*dt*/)
{
  Vector<Vector<int>> components(m_nlev, Vector<int>());

  for (int lev(0); lev < m_nlev; ++lev) {
    set_indexes(components[lev]);
  }

  if (m_specs[1].compare("sum") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = sum(lev, components[lev]);
    }

  } else if (m_specs[1].compare("min") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = min(lev, components[lev]);
    }

  } else if (m_specs[1].compare("max") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = max(lev, components[lev]);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
GeneralProperty::sum (const int lev,
                      const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const Real p_value = get_value(ptile_data, index, i);
      return p_value;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, 0.);

    Real l_sum = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMax(&l_sum, 1);

    result[var] = l_sum;
  }

  return result;
}


Vector<Real>
GeneralProperty::min (const int lev,
                      const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const Real p_value = get_value(ptile_data, index, i);
      return p_value;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op,
                                   R, std::numeric_limits<Real>::max());

    Real l_min = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMax(&l_min, 1);

    result[var] = l_min;
  }

  return result;
}


Vector<Real>
GeneralProperty::max (const int lev,
                      const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const Real p_value = get_value(ptile_data, index, i);
      return p_value;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op,
                                   R, std::numeric_limits<Real>::min());

    Real l_max = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMax(&l_max, 1);

    result[var] = l_max;
  }

  return result;
}


void
AveragedProperty::monitor (const Real& /*dt*/)
{
  Vector<Vector<int>> components(m_nlev, Vector<int>());

  for (int lev(0); lev < m_nlev; ++lev) {
    set_indexes(components[lev]);
  }

  if (m_specs[1].compare("average") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = average(lev, components[lev]);
    }

  } else if (m_specs[1].compare("standarddeviation") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = stddev(lev, components[lev]);
    }

  } else if (m_specs[1].compare("massweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = mass_weighted_average(lev, components[lev]);
    }

  } else if (m_specs[1].compare("volumeweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = volume_weighted_average(lev, components[lev]);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
AveragedProperty::average (const int lev,
                           const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = ReduceData<Real,Real>::Type;

    ReduceTuple default_values = {0., 0.};

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);

      return {statwt*p_value, statwt};
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, default_values);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    result[var] = numerator / denominator;
  }

  return result;
}


Vector<Real>
AveragedProperty::stddev (const int lev,
                          const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  const Vector<Real> avg = average(lev, indexes);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = ReduceData<Real,Real>::Type;

    ReduceTuple default_values = {0., 0.};

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    const Real d_avg = avg[var];

    auto R = [get_value,d_avg] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                                 const int index,
                                                 const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);

      const Real diff = p_value - d_avg;

      return {statwt*diff*diff, statwt};
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, default_values);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    result[var] = std::sqrt(numerator / denominator);
  }

  return result;
}


Vector<Real>
AveragedProperty::mass_weighted_average (const int lev,
                                         const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = ReduceData<Real,Real>::Type;

    ReduceTuple default_values = {0., 0.};

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);
      const Real mass = particle.rdata(SoArealData::mass);

      return {statwt*mass*p_value, statwt*mass};
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, default_values);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    result[var] = numerator / denominator;
  }

  return result;
}


Vector<Real>
AveragedProperty::volume_weighted_average (const int lev,
                                           const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = ReduceData<Real,Real>::Type;

    ReduceTuple default_values = {0., 0.};

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);
      const Real volume = particle.rdata(SoArealData::volume);

      return {statwt*volume*p_value, statwt*volume};
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, default_values);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    result[var] = numerator / denominator;
  }

  return result;
}


void
FlowRate::set_direction ()
{
  for (int dir(0); dir < AMREX_SPACEDIM; ++dir)
    if(std::abs(m_plane->lo(dir) - m_plane->hi(dir)) < 1.e-15) {
      m_direction = dir;
      break;
    }

  AMREX_ALWAYS_ASSERT(AMREX_D_TERM(std::abs(m_plane->length(m_direction)) < 1.e-15, &&
                                   m_plane->length((m_direction+1)%AMREX_SPACEDIM) > 1.e-15, &&
                                   m_plane->length((m_direction+2)%AMREX_SPACEDIM) > 1.e-15));
}


void
FlowRate::check_plane_is_ok () const
{
  AMREX_ALWAYS_ASSERT(std::abs(m_plane->volume()) < 1.e-15);
}


void
FlowRate::set_coordinate ()
{
  AMREX_ALWAYS_ASSERT(std::abs(m_plane->lo(m_direction) - m_plane->hi(m_direction)) < 1.e-15);

  m_coordinate = m_plane->lo(m_direction);
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
int
FlowRate::CrossesFlowPlane::operator() (const Real& position,
                                        const Real& velocity,
                                        const Real& plane_coordinate,
                                        const Real& dt) const
{
  // If dt is not strictly positive, return false
  if (dt < 1.e-15) {
    return 0;
  }

  const Real min_velocity = (plane_coordinate - position) / dt;

  const Real abs_min_velocity = std::abs(min_velocity);
  const Real abs_velocity = std::abs(velocity);

  const int it_crosses = (min_velocity*velocity > 0.) &&
    (std::abs(abs_min_velocity - abs_velocity) < 1.e-15 || abs_velocity > abs_min_velocity);

  return it_crosses;
}


void
FlowRate::monitor (const Real& dt)
{
  Vector<Vector<int>> components(m_nlev, Vector<int>());

  for (int lev(0); lev < m_nlev; ++lev) {
    set_indexes(components[lev]);
  }

  if (m_specs[1].compare("flowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = flow_rate(lev, components[lev], dt);
    }

  } else if (m_specs[1].compare("massweightedflowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = mass_weighted_flow_rate(lev, components[lev], dt);
    }

  } else if (m_specs[1].compare("volumeweightedflowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = volume_weighted_flow_rate(lev, components[lev], dt);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
FlowRate::flow_rate (const int lev,
                     const Vector<int>& indexes,
                     const Real dt)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    const int direction = m_direction;
    const Real plane_coordinate = m_coordinate;

    CrossesFlowPlane crossing_check;

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value,crossing_check,direction,plane_coordinate,dt]
      AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                        const int index,
                        const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);

      const Real pos = particle.pos(direction);
      const Real vel = particle.rdata(SoArealData::velx+plane_coordinate);

      const Real sign = std::abs(vel) < 1.e-15 ? 0. : (vel > 0. ? 1. : -1.);

      if (crossing_check(pos, vel, plane_coordinate, dt))
        return statwt*p_value*sign;
      else
        return 0.;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, 0.);

    Real l_flow_rate = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_flow_rate, 1);

    result[var] = l_flow_rate;
  }
  
  return result;
}


Vector<Real>
FlowRate::mass_weighted_flow_rate (const int lev,
                                   const Vector<int>& indexes,
                                   const Real dt)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    const int direction = m_direction;
    const Real plane_coordinate = m_coordinate;

    CrossesFlowPlane crossing_check;

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value,crossing_check,direction,plane_coordinate,dt]
      AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                        const int index,
                        const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);
      const Real mass = particle.rdata(SoArealData::mass);

      const Real pos = particle.pos(direction);
      const Real vel = particle.rdata(SoArealData::velx+plane_coordinate);

      const Real sign = std::abs(vel) < 1.e-15 ? 0. : (vel > 0. ? 1. : -1.);

      if (crossing_check(pos, vel, plane_coordinate, dt))
        return statwt*mass*p_value*sign;
      else
        return 0.;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, 0.);

    Real l_flow_rate = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_flow_rate, 1);

    result[var] = l_flow_rate;
  }
  
  return result;
}


Vector<Real>
FlowRate::volume_weighted_flow_rate (const int lev,
                                     const Vector<int>& indexes,
                                     const Real dt)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    const int direction = m_direction;
    const Real plane_coordinate = m_coordinate;

    CrossesFlowPlane crossing_check;

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value,crossing_check,direction,plane_coordinate,dt]
      AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                        const int index,
                        const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);
      const Real volume = particle.rdata(SoArealData::volume);

      const Real pos = particle.pos(direction);
      const Real vel = particle.rdata(SoArealData::velx+plane_coordinate);

      const Real sign = std::abs(vel) < 1.e-15 ? 0. : (vel > 0. ? 1. : -1.);

      if (crossing_check(pos, vel, plane_coordinate, dt))
        return statwt*volume*p_value*sign;
      else
        return 0.;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, 0.);

    Real l_flow_rate = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_flow_rate, 1);

    result[var] =  l_flow_rate;
  }
  
  return result;
}

} // end namespace LagrangianMonitor
