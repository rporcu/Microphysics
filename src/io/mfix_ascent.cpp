#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <mfix_rw.H>
#include <mfix_fluid_parms.H>
#include <mfix_solids_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>

#ifdef AMREX_USE_ASCENT
#include <AMReX_Conduit_Blueprint.H>
#include <ascent.hpp>
#endif

using namespace amrex;


namespace MfixIO {

void
MfixRW::WriteAscentFile (int nstep, const Real time) const
{
#ifdef AMREX_USE_ASCENT
  BL_PROFILE("mfix::WriteAscentFile()");

  //amrex::Print() << "Writing Ascent output\n";

  ParmParse pp("ascent");

  std::string ascent_fluid_actions_yaml {""};
  pp.query("fluid_actions", ascent_fluid_actions_yaml);

  std::string ascent_particle_actions_yaml {""};
  pp.query("particle_actions", ascent_particle_actions_yaml);


  if (fluid.solve && !ascent_fluid_actions_yaml.empty()) {

    Vector<std::string> pltFldNames;

    pltFldNames.push_back("u_g");
    pltFldNames.push_back("v_g");
    pltFldNames.push_back("w_g");

    pltFldNames.push_back("ep_g");
    pltFldNames.push_back("volfrac");

    const int ngrow = 0;
    int ncomp = 5;

    // Temperature in fluid
    if (fluid.solve_enthalpy) {
      pltFldNames.push_back("T_g");
      ncomp += 1;
    }

    if ( fluid.solve_species ) {
      for (std::string specie: fluid.species) {
        pltFldNames.push_back("Xg_"+specie);
        ncomp += 1;
      }
    }

    AMREX_ALWAYS_ASSERT(pltFldNames.size() == ncomp);

    amrex::Vector<int> level_steps;
    Vector< MultiFab* > mf(nlev);

    for (int lev = 0; lev < nlev; ++lev) {

      level_steps.push_back(nstep);

      mf[lev] = new MultiFab(grids[lev], dmap[lev], ncomp, ngrow,  MFInfo(), *ebfactory[lev]);

      MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->vel_g), 0, 0, 1, 0);
      MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->vel_g), 1, 1, 1, 0);
      MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->vel_g), 2, 2, 1, 0);
      MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->ep_g ), 0, 3, 1, 0);
      MultiFab::Copy(*mf[lev], ebfactory[lev]->getVolFrac(), 0, 4, 1, 0);

      int lc=5;
      if (fluid.solve_enthalpy) {
        MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->T_g), 0, lc, 1, 0);
        lc += 1;
      }

      // Fluid species mass fractions
      if (fluid.solve_species) {
        MultiFab::Copy(*mf[lev], *m_leveldata[lev]->X_gk, 0, lc, fluid.nspecies, 0);
        lc += fluid.nspecies;
      }

      amrex::EB_set_covered(*mf[lev], 0.0);
    }

    conduit::Node bp_mesh;

    amrex::MultiLevelToBlueprint(nlev, amrex::GetVecOfConstPtrs(mf),
        pltFldNames, geom, time, level_steps, ref_ratio, bp_mesh);



    // for the MPI case, provide the mpi comm
    ascent::Ascent ascent;
    conduit::Node opts;
    opts["exceptions"] = "catch";
    opts["actions_file"] = ascent_fluid_actions_yaml;
#ifdef BL_USE_MPI
    opts["mpi_comm"] = MPI_Comm_c2f(ParallelDescriptor::Communicator());
#endif
    ascent.open(opts);
    ascent.publish(bp_mesh);
    conduit::Node actions;
    ascent.execute(actions);
    ParallelDescriptor::Barrier();
    ascent.close();

    for (int lev = 0; lev < nlev; ++lev) {
      delete mf[lev];
    }
  }


  if (( DEM::solve || PIC::solve ) && !ascent_particle_actions_yaml.empty()) {

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

    // for the MPI case, provide the mpi comm
    ascent::Ascent ascent;
    conduit::Node opts;

    opts["actions_file"] = ascent_particle_actions_yaml;
    opts["exceptions"] = "catch";

#ifdef BL_USE_MPI
    opts["mpi_comm"] = MPI_Comm_c2f(ParallelDescriptor::Communicator());
#endif

    ascent.open(opts);

    conduit::Node bp_particles;

    amrex::ParticleContainerToBlueprint(*pc,
              real_comp_names, int_comp_names, bp_particles);

    ascent.publish(bp_particles);

    conduit::Node actions;
    ascent.execute(actions);
    ParallelDescriptor::Barrier();
    ascent.close();
  }
#else
  amrex::ignore_unused(nstep);
  amrex::ignore_unused(time);
#endif
}

} // end namespace MfixIO
