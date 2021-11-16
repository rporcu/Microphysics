#ifdef MFIX_CATALYST

#include <string>

#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <mfix.H>
#include <mfix_fluid_parms.H>
#include <mfix_solids_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>

#ifdef AMREX_USE_CONDUIT
#include <AMReX_Conduit_Blueprint.H>
#include "catalyst.hpp"
 #include "conduit_cpp_to_c.hpp"
#endif

void
mfix::RunCatalystAdaptor ( int nstep, Real time )
{
#ifdef AMREX_USE_CONDUIT
  BL_PROFILE("mfix::RunCatalystAdaptor()");

  if ( DEM::solve || PIC::solve ) {

    Vector<std::string> real_comp_names;
    Vector<std::string>  int_comp_names;

    real_comp_names.push_back("radius");
    real_comp_names.push_back("volume");
    real_comp_names.push_back("mass");
    real_comp_names.push_back("density");

    if(DEM::solve){
      real_comp_names.push_back("omoi");
    } else {
      real_comp_names.push_back("ep_s");
    }

    real_comp_names.push_back("velx");
    real_comp_names.push_back("vely");
    real_comp_names.push_back("velz");

    if(DEM::solve){
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

    if (solids.solve_species)
      for(auto species: solids.species)
        real_comp_names.push_back("X_"+species+"_s");

    if (solids.solve_species && reactions.solve)
      for(auto species: solids.species)
        real_comp_names.push_back("chem_ro_txfr_"+species);

    if (reactions.solve) {
      real_comp_names.push_back("chem_velx_txfr");
      real_comp_names.push_back("chem_vely_txfr");
      real_comp_names.push_back("chem_velz_txfr");
    }

    if (reactions.solve)
      real_comp_names.push_back("chem_h_txfr");

    int_comp_names.push_back("phase");
    int_comp_names.push_back("state");

    MFIXParticleContainer* pc = getParticleContainer();

    conduit::Node bp_particles;

    // add time/cycle information
    auto& state = bp_particles["catalyst/state"];
    state["timestep"].set(nstep);
    state["time"].set(time);

    auto& channel = bp_particles["catalyst/channels/particles"];
    channel["type"].set("multimesh");
    amrex::ParticleContainerToBlueprint(*pc, real_comp_names, int_comp_names, channel["data"]);

    bp_particles.print();
    catalyst_status err = catalyst_execute(conduit::c_node(&bp_particles));
    if (err != catalyst_status_ok)
    {
      std::cerr << "Failed to execute Catalyst: " << err << std::endl;
    }
  }
#endif
}

#endif // MFIX_CATALYST
