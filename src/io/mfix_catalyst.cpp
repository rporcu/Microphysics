#ifdef MFIX_CATALYST

#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <mfix.H>
#include <mfix_rw.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>

#ifdef AMREX_USE_CONDUIT
#include <AMReX_Conduit_Blueprint.H>
#include "catalyst.hpp"
#include "conduit_cpp_to_c.hpp"
#endif

namespace internal
{

void EmptyParticleData(const std::vector<std::string>& real_fields,
                       const std::vector<std::string>& int_fields,
                       conduit::Node& node)
{
  const std::string topology_name = "particles";
  // make a dummy node for catalyst
  const std::string& patch_name = amrex::Concatenate("domain_x", 0, 6);
  conduit::Node &patch = node[patch_name];

  patch["state/domain_id"] = 0;

  std::string coordset_name = topology_name + "+coords";
  conduit::Node &n_coords = patch["coordsets"][coordset_name];
  n_coords["type"] = "explicit";

  // create an explicit points topology
  conduit::Node &n_topo = patch["topologies"][topology_name];
  n_topo["coordset"] = coordset_name;
  n_topo["type"] = "unstructured";
  n_topo["elements/shape"] = "point";
  n_topo["elements/connectivity"].set(conduit::DataType::c_int(0));

  n_coords["values/x"].set(conduit::DataType::c_int(0));
  n_coords["values/y"].set(conduit::DataType::c_int(0));
  n_coords["values/z"].set(conduit::DataType::c_int(0));

  conduit::Node &n_fields = patch["fields"];

  // id field
  conduit::Node &n_f_id = n_fields[topology_name + "_id"];
  n_f_id["topology"] = topology_name;
  n_f_id["association"] = "element";
  n_f_id["values"].set(conduit::DataType::c_int(0));

  // cpu field
  conduit::Node &n_f_cpu = n_fields[topology_name + "_cpu"];
  n_f_cpu["topology"] = topology_name;
  n_f_cpu["association"] = "element";
  n_f_cpu["values"].set(conduit::DataType::c_int(0));

  for (auto& rfield : real_fields)
  {
    conduit::Node & n_f = n_fields[rfield];
    n_f["topology"] = topology_name;
    n_f["association"] = "element";
    n_f["values"].set(conduit::DataType::c_int(0));
  }
  for (auto& ifield : int_fields)
  {
    conduit::Node & n_f = n_fields[ifield];
    n_f["topology"] = topology_name;
    n_f["association"] = "element";
    n_f["values"].set(conduit::DataType::c_int(0));
  }
}

} // namespace internal


void
MFIXReadWrite::RunCatalystAdaptor(int nstep, const Real time) const
{
#ifdef AMREX_USE_CONDUIT
  BL_PROFILE("mfix::RunCatalystAdaptor()");

  conduit::Node node;
  auto& state = node["catalyst/state"];
  state["timestep"].set(nstep);
  state["time"].set(time);

  amrex::Vector<std::string> pltFldNames;
  amrex::Vector< MultiFab* > mf(nlev);

  if (fluid.solve())
  {
    auto& channel = node["catalyst/channels/mesh"];
    channel["type"].set_string("amrmesh");
    auto& meshData = channel["data"];

    const int ngrow = 0;
    int ncomp = 5;

    pltFldNames.push_back("u_g");
    pltFldNames.push_back("v_g");
    pltFldNames.push_back("w_g");
    pltFldNames.push_back("ep_g");
    pltFldNames.push_back("volfrac");

    // Temperature in fluid
    if (fluid.solve_enthalpy()) {
      pltFldNames.push_back("T_g");
      ncomp += 1;
    }

    if ( fluid.solve_species()) {
      for (std::string specie: fluid.species_names()) {
        pltFldNames.push_back("Xg_"+specie);
        ncomp += 1;
      }
    }

    AMREX_ALWAYS_ASSERT(pltFldNames.size() == ncomp);

    amrex::Vector<int> level_steps(nlev);

    for (int lev = 0; lev < nlev; ++lev) {

      level_steps[lev] = nstep;

      mf[lev] = new MultiFab(grids[lev], dmap[lev], ncomp, ngrow,  MFInfo(), *ebfactory[lev]);

      MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->vel_g), 0, 0, 1, 0);
      MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->vel_g), 1, 1, 1, 0);
      MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->vel_g), 2, 2, 1, 0);
      MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->ep_g ), 0, 3, 1, 0);
      MultiFab::Copy(*mf[lev], ebfactory[lev]->getVolFrac(), 0, 4, 1, 0);

      int lc=5;
      if (fluid.solve_enthalpy()) {
        MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->T_g), 0, lc, 1, 0);
        lc += 1;
      }

      // Fluid species mass fractions
      if (fluid.solve_species()) {
        MultiFab::Copy(*mf[lev], *m_leveldata[lev]->X_gk, 0, lc, fluid.nspecies(), 0);
        lc += fluid.nspecies();
      }

      amrex::EB_set_covered(*mf[lev], 0.0);
    }

    //std::cout << "Getting mesh data" << std::endl;
    amrex::MultiLevelToBlueprint(nlev, amrex::GetVecOfConstPtrs(mf),
        pltFldNames, geom, time, level_steps, ref_ratio, meshData);
  }


  if ( m_dem.solve() || m_pic.solve() )
  {
    auto& channel = node["catalyst/channels/particles"];
    channel["type"].set_string("multimesh");
    auto& particleData = channel["data"];

    Vector<std::string> real_comp_names;
    Vector<std::string>  int_comp_names;

    real_comp_names.push_back("radius");
    real_comp_names.push_back("volume");
    real_comp_names.push_back("mass");
    real_comp_names.push_back("density");

    if(m_dem.solve()){
      real_comp_names.push_back("omoi");
    } else {
      real_comp_names.push_back("ep_s");
    }

    real_comp_names.push_back("velx");
    real_comp_names.push_back("vely");
    real_comp_names.push_back("velz");

    if(m_dem.solve()){
      real_comp_names.push_back("omegax");
      real_comp_names.push_back("omegay");
      real_comp_names.push_back("omegaz");
    } else {
      real_comp_names.push_back("grad_tau_x");
      real_comp_names.push_back("grad_tau_y");
      real_comp_names.push_back("grad_tau_z");
    }

    real_comp_names.push_back("statwt");
    real_comp_names.push_back("dragcoeff");
    real_comp_names.push_back("dragx");
    real_comp_names.push_back("dragy");
    real_comp_names.push_back("dragz");

    real_comp_names.push_back("c_ps");
    real_comp_names.push_back("temperature");
    real_comp_names.push_back("convection");

    if (solids.solve_species())
      for(auto species: solids.species_names())
        real_comp_names.push_back("X_"+species+"_s");

    if (solids.solve_species() && reactions.solve())
      for(auto species: solids.species_names())
        real_comp_names.push_back("chem_ro_txfr_"+species);

    if (reactions.solve()) {
      real_comp_names.push_back("chem_velx_txfr");
      real_comp_names.push_back("chem_vely_txfr");
      real_comp_names.push_back("chem_velz_txfr");
    }

    if (reactions.solve())
      real_comp_names.push_back("chem_h_txfr");

    int_comp_names.push_back("phase");
    int_comp_names.push_back("state");
#if MFIX_POLYDISPERSE
    int_comp_names.push_back("ptype");
#endif
    conduit::Node amrexParticles;
    amrex::ParticleContainerToBlueprint(*pc,
              real_comp_names, int_comp_names, amrexParticles);
    particleData.update(amrexParticles);
    if(!particleData.dtype().is_object())
    {
      internal::EmptyParticleData(real_comp_names, int_comp_names, particleData);
    }
  }

  // run catalyst
  catalyst_status err = catalyst_execute(conduit::c_node(&node));
  if (err != catalyst_status_ok)
  {
    std::cerr << "Failed to execute Catalyst: " << err << std::endl;
  }

  // cleanup intermediate mf
  if (fluid.solve())
  {
    for (int lev = 0; lev < nlev; ++lev) {
      delete mf[lev];
    }
  }
#endif
}

#endif // MFIX_CATALYST
