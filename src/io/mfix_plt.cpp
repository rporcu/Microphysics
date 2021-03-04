#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_ParmParse.H>

#include <mfix.H>
#include <mfix_fluid_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_pic_parms.H>

namespace
{
    const std::string level_prefix {"Level_"};
}

void
mfix::InitIOPltData ()
{
  if (ooo_debug) amrex::Print() << "InitIOPltData" << std::endl;

  // Variables to simplify checkpoint IO

  pltVarCount = 0;

  ParmParse pp("amr");
  if (FLUID::solve)
    {

      pp.query("plt_vel_g",   plt_vel_g  );
      pp.query("plt_ep_g",    plt_ep_g   );
      pp.query("plt_p_g",     plt_p_g    );
      pp.query("plt_ro_g",    plt_ro_g   );
      pp.query("plt_MW_g",    plt_MW_g   );
      pp.query("plt_h_g",     plt_h_g    );
      pp.query("plt_T_g",     plt_T_g    );
      pp.query("plt_trac",    plt_trac   );
      pp.query("plt_cp_g",    plt_cp_g   );
      pp.query("plt_k_g",     plt_k_g    );
      pp.query("plt_mu_g",    plt_mu_g   );
      pp.query("plt_diveu",   plt_diveu  );
      pp.query("plt_vort",    plt_vort   );
      pp.query("plt_volfrac", plt_volfrac);
      pp.query("plt_gradp_g", plt_gradp_g);
      pp.query("plt_X_g",     plt_X_gk   );
      pp.query("plt_D_g",     plt_D_gk   );
      pp.query("plt_cp_gk",   plt_cp_gk  );
      pp.query("plt_h_gk",    plt_h_gk   );
      pp.query("plt_ro_txfr", plt_ro_txfr);

      // Special test for CCSE regression test. Override all individual
      // flags and save all data to plot file.

      int plt_ccse_regtest = 0;
      pp.query("plt_regtest", plt_ccse_regtest);

      if (plt_ccse_regtest != 0) {
        plt_vel_g   = 1;
        plt_ep_g    = 1;
        plt_p_g     = 0;
        plt_ro_g    = 1;
        plt_MW_g    = 0;
        plt_h_g     = 1;
        plt_T_g     = 1;
        plt_trac    = 1;
        plt_cp_g    = 1;
        plt_k_g     = 1;
        plt_mu_g    = 1;
        plt_vort    = 1;
        plt_diveu   = 1;
        plt_volfrac = 1;
        plt_gradp_g = 1;
        plt_X_gk    = 1;
        plt_D_gk    = 1;
        plt_cp_gk   = 0;
        plt_h_gk    = 0;
        plt_ro_txfr = 0;
      }

      // Count the number of variables to save.
      if( plt_vel_g   == 1) pltVarCount += 3;
      if( plt_gradp_g == 1) pltVarCount += 3;
      if( plt_ep_g    == 1) pltVarCount += 1;
      if( plt_p_g     == 1) pltVarCount += 1;
      if( plt_ro_g    == 1) pltVarCount += 1;
      if( plt_MW_g    == 1) pltVarCount += 1;
      if( plt_trac    == 1) pltVarCount += 1;
      if( plt_mu_g    == 1) pltVarCount += 1;
      if( plt_vort    == 1) pltVarCount += 1;
      if( plt_diveu   == 1) pltVarCount += 1;
      if( plt_volfrac == 1) pltVarCount += 1;

      if (advect_enthalpy) {
        if( plt_T_g  == 1) pltVarCount += 1;
        if( plt_cp_g == 1) pltVarCount += 1;
        if( plt_k_g  == 1) pltVarCount += 1;
        if( plt_h_g  == 1) pltVarCount += 1;
      }

      if (FLUID::solve_species) {
        if( plt_X_gk == 1)  pltVarCount += FLUID::nspecies;
        if( plt_D_gk == 1)  pltVarCount += FLUID::nspecies;

        if (advect_enthalpy) {
        if( plt_cp_gk == 1) pltVarCount += FLUID::nspecies;
        if( plt_h_gk == 1)  pltVarCount += FLUID::nspecies;
        }
      }

      if (FLUID::solve_species and REACTIONS::solve) {
        if ( plt_ro_txfr == 1) pltVarCount += FLUID::nspecies;
      }
    }

    if(DEM::solve or PIC::solve)
    {
      int plt_ccse_regtest = 0;
      pp.query("plt_regtest", plt_ccse_regtest);

      if (SOLIDS::solve_species) {
        const int size = AoSrealData::count + SoArealData::count +
                         SoAspeciesData::count*SOLIDS::nspecies;
        write_real_comp.resize(size, 0);
      }

      if (SOLIDS::solve_species and REACTIONS::solve) {
        const int size = AoSrealData::count + SoArealData::count +
                         SoAspeciesData::count*SOLIDS::nspecies +
                         SoAreactionsData::count*SOLIDS::nspecies*REACTIONS::nreactions;
        write_real_comp.resize(size, 0);
      }

      // All flags are true by default so we only need to turn off the
      // variables we don't want if not doing CCSE regression tests.
      if (plt_ccse_regtest == 0)
      {
        int gap = AoSrealData::count;

        int input_value = 0;
        pp.query("plt_radius",   input_value );
        write_real_comp[gap+SoArealData::radius] = input_value;

        input_value = 0;
        pp.query("plt_volume",   input_value );
        write_real_comp[gap+SoArealData::volume] = input_value;

        input_value = 0;
        pp.query("plt_mass",     input_value );
        write_real_comp[gap+SoArealData::mass] = input_value;

        input_value = 0;
        pp.query("plt_ro_p",     input_value );
        write_real_comp[gap+SoArealData::density] = input_value;

        input_value = 0;
        pp.query("plt_omoi"   ,  input_value );
        write_real_comp[gap+SoArealData::oneOverI] = input_value;

        input_value = 1;
        pp.query("plt_vel_p"   ,  input_value );
        write_real_comp[gap+SoArealData::velx] = input_value;
        write_real_comp[gap+SoArealData::vely] = input_value;
        write_real_comp[gap+SoArealData::velz] = input_value;

        input_value = 0;
        pp.query("plt_omega_p",   input_value );
        write_real_comp[gap+SoArealData::omegax] = input_value;
        write_real_comp[gap+SoArealData::omegay] = input_value;
        write_real_comp[gap+SoArealData::omegaz] = input_value;

        input_value = 0;
        pp.query("plt_statwt",   input_value );
        write_real_comp[gap+SoArealData::statwt] = input_value;

        input_value = 0;
        pp.query("plt_drag_p",   input_value );
        write_real_comp[gap+SoArealData::dragcoeff] = input_value;  // drag coeff
        write_real_comp[gap+SoArealData::dragx] = input_value;  // dragx
        write_real_comp[gap+SoArealData::dragy] = input_value;  // dragy
        write_real_comp[gap+SoArealData::dragz] = input_value;  // dragz

        input_value = 0;
        pp.query("plt_cp_s",  input_value);
        write_real_comp[gap+SoArealData::c_ps] = input_value;  // specific heat

        input_value = 0;
        pp.query("plt_T_p",   input_value );
        write_real_comp[gap+SoArealData::temperature] = input_value;  // temperature

        input_value = 0;
        pp.query("plt_convection", input_value );
        write_real_comp[gap+SoArealData::convection] = input_value;  // heat transfer coefficient

        gap = AoSrealData::count + SoArealData::count;

        input_value = 0;
        pp.query("plt_species_p",   input_value );
        if (SOLIDS::solve_species)
        {
          for(int n(0); n < SOLIDS::nspecies; ++n)
            write_real_comp[gap+n] = input_value;
        }

        gap = AoSrealData::count + SoArealData::count + SoAspeciesData::count*SOLIDS::nspecies;

        input_value = 0;
        pp.query("plt_ro_rates_p", input_value );
        if (SOLIDS::solve_species and REACTIONS::solve)
        {
          for(int n_s(0); n_s < SOLIDS::nspecies; ++n_s)
            for(int q(0); q < REACTIONS::nreactions; ++q) {
              const int comp = gap + n_s*REACTIONS::nreactions + q;
              write_real_comp[comp] = input_value;
            }
        }

        gap = AoSintData::count;

        input_value = 0;
        pp.query("plt_phase",   input_value );
        write_int_comp[gap+SoAintData::phase] = input_value;

        input_value = 0;
        pp.query("plt_state",   input_value );
        write_int_comp[gap+SoAintData::state] = input_value;

      }

    }

}

void
mfix::WritePlotFile (std::string& plot_file, int nstep, Real time )
{
    // If we've already written this plotfile, don't do it again!
    if (nstep == last_plt) return;

    // Now set last_plt to nstep ...
    last_plt = nstep;

    BL_PROFILE("mfix::WritePlotFile()");

    const std::string& plotfilename = amrex::Concatenate(plot_file,nstep);

    amrex::Print() << "  Writing plotfile " << plotfilename <<  " at time " << time << std::endl;


    if (pltVarCount > 0) {

      const int ngrow = 0;

      Vector<std::string> pltFldNames;
      Vector< std::unique_ptr<MultiFab> > mf(nlev);

      // Velocity components
      if( plt_vel_g   == 1) {
        pltFldNames.push_back("u_g");
        pltFldNames.push_back("v_g");
        pltFldNames.push_back("w_g");
      }

      // Pressure gradient
      if( plt_gradp_g == 1) {
        pltFldNames.push_back("gpx");
        pltFldNames.push_back("gpy");
        pltFldNames.push_back("gpz");
      }

      // Fluid volume fraction
      if( plt_ep_g == 1 )
        pltFldNames.push_back("ep_g");

      // Fluid pressure
      if( plt_p_g == 1 )
        pltFldNames.push_back("p_g");

      // Fluid density
      if( plt_ro_g == 1 )
        pltFldNames.push_back("ro_g");

      // Fluid molecular weight
      if( plt_MW_g == 1 )
        pltFldNames.push_back("MW_g");

      // Fluid enthalpy
      if( advect_enthalpy and plt_h_g == 1 )
        pltFldNames.push_back("h_g");

      // Temperature in fluid
      if( advect_enthalpy and plt_T_g == 1 )
        pltFldNames.push_back("T_g");

      // Tracer in fluid
      if( plt_trac == 1 )
        pltFldNames.push_back("trac");

      // Specific heat
      if( advect_enthalpy and plt_cp_g == 1 )
        pltFldNames.push_back("cp_g");

      // Thermal conductivity
      if( advect_enthalpy and plt_k_g == 1 )
        pltFldNames.push_back("k_g");

      // Fluid viscosity
      if( plt_mu_g == 1 )
        pltFldNames.push_back("mu_g");

      // vorticity
      if( plt_vort == 1 )
        pltFldNames.push_back("vort");

      // div(ep_g.u)
      if( plt_diveu == 1 )
        pltFldNames.push_back("diveu");

      // EB cell volume fraction
      if( plt_volfrac == 1 )
        pltFldNames.push_back("volfrac");

      // Fluid species mass fractions
      if(FLUID::solve_species and plt_X_gk == 1)
        for(std::string specie: FLUID::species)
          pltFldNames.push_back("X_"+specie+"_g");

      // Fluid species mass diffusivities
      if(FLUID::solve_species and plt_D_gk == 1)
        for(std::string specie: FLUID::species)
          pltFldNames.push_back("D_"+specie+"_g");

      // Fluid species specific heat
      if(FLUID::solve_species and advect_enthalpy and plt_cp_gk == 1)
        for(std::string specie: FLUID::species)
          pltFldNames.push_back("cp_"+specie+"_g");

      // Fluid species enthalpy
      if(FLUID::solve_species and advect_enthalpy and plt_h_gk == 1)
        for(std::string specie: FLUID::species)
          pltFldNames.push_back("h_"+specie+"_g");

      // Fluid species density reaction rates
      if (FLUID::solve_species and REACTIONS::solve and plt_ro_txfr == 1)
        for(std::string specie: FLUID::species)
          pltFldNames.push_back("ro_txfr_"+specie+"_g");

      for (int lev = 0; lev < nlev; ++lev)
      {
        // Multifab to hold all the variables -- there can be only one!!!!
        const int ncomp = pltVarCount;
        mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow,  MFInfo(), *ebfactory[lev]));

        int lc=0;

        // Velocity components
        if( plt_vel_g == 1 ) {
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->vel_g), 0, lc  , 1, 0);
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->vel_g), 1, lc+1, 1, 0);
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->vel_g), 2, lc+2, 1, 0);
          lc += 3;
        }

        // Pressure gradient
        if( plt_gradp_g == 1 ) {
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->gp, 0, lc  , 1, 0);
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->gp, 1, lc+1, 1, 0);
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->gp, 2, lc+2, 1, 0);
          lc += 3;
        }

        // Fluid volume fraction
        if( plt_ep_g == 1 ) {
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->ep_g, 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid pressure
        if( plt_p_g == 1 ) {
          MultiFab p_nd(m_leveldata[lev]->p_g->boxArray(), dmap[lev], 1, 0);
          p_nd.setVal(0.);
          MultiFab::Copy(p_nd, *m_leveldata[lev]->p_g, 0, 0, 1, 0);
          MultiFab::Add (p_nd, *m_leveldata[lev]->p0_g, 0, 0, 1, 0);
          amrex::average_node_to_cellcenter(*mf[lev], lc, p_nd, 0, 1);
          lc += 1;
        }

        // Fluid density
        if( plt_ro_g == 1 ) {
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->ro_g), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid molecular weight
        if( plt_MW_g == 1 ) {
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->MW_g), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid enthalpy
        if( advect_enthalpy and plt_h_g == 1 ) {
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->h_g), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid temperature
        if( advect_enthalpy and plt_T_g == 1 ) {
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->T_g), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid tracer
        if( plt_trac == 1 ) {
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->trac), 0, lc, 1, 0);
          lc += 1;
        }

        // Specific heat
        if( advect_enthalpy and plt_cp_g == 1 ) {
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->cp_g, 0, lc, 1, 0);
          lc += 1;
        }

        // Specific heat
        if( advect_enthalpy and plt_k_g == 1 ) {
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->k_g, 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid viscosity
        if( plt_mu_g == 1 ) {
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->mu_g, 0, lc, 1, 0);
          lc += 1;
        }

        // vorticity
        if( plt_vort == 1 ) {
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->vort, 0, lc, 1, 0);
          lc += 1;
        }

        // div(ep_g.u)
        if( plt_diveu == 1 ) {
          amrex::average_node_to_cellcenter(*mf[lev], lc, *m_leveldata[lev]->diveu, 0, 1);
          lc += 1;
        }

        // EB cell volume fraction
        if( plt_volfrac == 1 ) {
          if (ebfactory[lev]) {
            MultiFab::Copy(*mf[lev], ebfactory[lev]->getVolFrac(), 0, lc, 1, 0);
          } else {
            mf[lev]->setVal(1.0,lc,1,0);
          }
          lc += 1;
        }

        // Fluid species mass fractions
        if(FLUID::solve_species and plt_X_gk == 1 ) {
          for(int n(0); n < FLUID::nspecies; n++) {
            MultiFab::Copy(*mf[lev], *m_leveldata[lev]->X_gk, n, lc+n, 1, 0);
          }
          lc += FLUID::nspecies;
        }

        // Species mass fraction
        if(FLUID::solve_species and plt_D_gk == 1 ) {
          for(int n(0); n < FLUID::nspecies; n++) {
            MultiFab::Copy(*mf[lev], *m_leveldata[lev]->D_gk, n, lc+n, 1, 0);
          }
          lc += FLUID::nspecies;
        }

        // Fluid species specific heat
        if(FLUID::solve_species and advect_enthalpy and plt_cp_gk == 1 ) {
          for(int n(0); n < FLUID::nspecies; n++) {
            MultiFab::Copy(*mf[lev], *m_leveldata[lev]->cp_gk, n, lc+n, 1, 0);
          }
          lc += FLUID::nspecies;
        }

        // Fluid species enthalpy
        if(FLUID::solve_species and plt_h_gk == 1 ) {
          for(int n(0); n < FLUID::nspecies; n++) {
            MultiFab::Copy(*mf[lev], *m_leveldata[lev]->h_gk, n, lc+n, 1, 0);
          }
          lc += FLUID::nspecies;
        }

        if (FLUID::solve_species and REACTIONS::solve and plt_ro_txfr == 1) {
          for(int n(0); n < FLUID::nspecies; n++) {
            MultiFab::Copy(*mf[lev], *m_leveldata[lev]->ro_gk_txfr, n, lc+n, 1, 0);
          }
          lc += FLUID::nspecies;
        }
      }

      // Cleanup places where we have no data.
      Vector<const MultiFab*> mf2(nlev);
      for (int lev = 0; lev < nlev; ++lev) {
        EB_set_covered(*mf[lev], 0.0);
        mf2[lev] = mf[lev].get();
      }

      Vector<int> istep;
      istep.resize(nlev,nstep);
      amrex::WriteMultiLevelPlotfile(plotfilename, nlev, mf2, pltFldNames,
                                     Geom(), time, istep, refRatio());


      // no fluid
    } else {

      // Some post-processing tools (such as yt) might still need some basic
      // MultiFab header information to function. We provide this here by
      // creating an "empty" plotfile header (which essentially only contains
      // the BoxArray information). Particle data is saved elsewhere.

      Vector< std::unique_ptr<MultiFab> > mf(finest_level+1);
      Vector<std::string>  names;
      // NOTE: leave names vector empty => header should reflect nComp = 0
      //names.insert(names.end(), "placeholder");

      // Create empty MultiFab containing the right BoxArray (NOTE: setting
      // nComp = 1 here to avoid assertion fail in debug build).
      for (int lev = 0; lev <= finest_level; ++lev)
        mf[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0));

      Vector<const MultiFab*> mf2(finest_level+1);

      for (int lev = 0; lev <= finest_level; ++lev)
        mf2[lev] = mf[lev].get();

      // Write only the Headers corresponding to the "empty" mf/mf2 MultiFabs
      Vector<int> istep;
      istep.resize(nlev,nstep);
      amrex::WriteMultiLevelPlotfileHeaders(plotfilename, finest_level+1, mf2, names,
                                            Geom(), time, istep, refRatio());

    }

    WriteJobInfo(plotfilename);

    if ( DEM::solve or PIC::solve )
    {
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

        if (SOLIDS::solve_species)
          for(auto species: SOLIDS::species)
            real_comp_names.push_back("X_"+species+"_s");

        if (SOLIDS::solve_species and REACTIONS::solve)
          for(auto species: SOLIDS::species)
            for(auto reaction: REACTIONS::reactions)
              real_comp_names.push_back(reaction+"_"+species+"_txfr");

        int_comp_names.push_back("phase");
        int_comp_names.push_back("state");

        pc->WritePlotFile(plotfilename, "particles", write_real_comp,
                          write_int_comp, real_comp_names, int_comp_names);

    }

}

void mfix::WriteStaticPlotFile (const std::string & plotfilename) const
{
    BL_PROFILE("mfix::WriteStaticPlotFile()");

    Print() << "  Writing static quantities " << plotfilename << std::endl;

    /****************************************************************************
     *                                                                          *
     * Static (un-changing variables):                                          *
     *     1. level-set data                                                    *
     *     2. volfrac (from EB) data                                            *
     *                                                                          *
     ***************************************************************************/

    Vector<std::string> static_names = {"level_sets", "volfrac"};
    Vector< Vector< MultiFab const* > > static_vars = { amrex::GetVecOfConstPtrs(level_sets) };

    const int ngrow = 0;
    const int ncomp = static_names.size();


    /****************************************************************************
     *                                                                          *
     * Collect variables together into a single multi-component MultiFab        *
     *                                                                          *
     ***************************************************************************/

    Vector<std::unique_ptr<MultiFab>> mf(nlev);
    Vector<const MultiFab *>          mf_ptr(nlev);

    for (int lev = 0; lev < nlev; lev++)
    {
        mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow, MFInfo(),
                                   * particle_ebfactory[lev]));

        // Don't iterate over all ncomp => last component is for volfrac
        for (int dcomp = 0; dcomp < ncomp - 1; dcomp++)
        {
            const BoxArray nd_ba = amrex::convert(grids[lev], IntVect::TheNodeVector());
            MultiFab mf_loc = MFUtil::regrid(nd_ba, dmap[lev], *static_vars[dcomp][lev], true);
            amrex::average_node_to_cellcenter(* mf[lev], dcomp, mf_loc, 0, 1, ngrow);
        }

        if (ebfactory[lev]) {
            EBFArrayBoxFactory ebf(* eb_levels[lev], geom[lev], grids[lev], dmap[lev],
                                   {nghost_eb_basic(), nghost_eb_volume(),
                                    nghost_eb_full()}, m_eb_support_level);

            MultiFab::Copy(* mf[lev], ebf.getVolFrac(), 0, ncomp - 1, 1, ngrow);

        } else {
            // setVal (value_type val, int comp, int num_comp, int nghost=0)
            mf[lev]->setVal(1.0, ncomp - 1, 1, ngrow);
        }
    }

    for (int lev = 0; lev < nlev; ++lev)
    {
        // Don't do this (below) as it zeros out the covered cells...
        // EB_set_covered(* mf[lev], 0.0);
        mf_ptr[lev] = mf[lev].get();
    }

    Real time = 0.;
    Vector<int> istep;
    istep.resize(nlev,0);
    amrex::WriteMultiLevelPlotfile(plotfilename, nlev, mf_ptr, static_names,
                                   Geom(), time, istep, refRatio());

    WriteJobInfo(plotfilename);

    Print() << "  Done writing static quantities " << plotfilename << std::endl;
}
