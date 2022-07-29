#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_ParmParse.H>
#include <AMReX_EBFArrayBox.H>

#include <mfix_rw.H>
#include <mfix_pc.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>

using namespace amrex;

//namespace
//{
//    const std::string level_prefix {"Level_"};
//}

namespace MfixIO {

void
MfixRW::InitIOPltData ()
{
  if (ooo_debug) amrex::Print() << "InitIOPltData" << std::endl;

  // Variables to simplify checkpoint IO

  pltVarCount = 0;

  ParmParse pp("amr");

  if (fluid.solve()) {

      pp.query("plt_vel_g",     plt_vel_g    );
      pp.query("plt_ep_g",      plt_ep_g     );
      pp.query("plt_p_g",       plt_p_g      );
      pp.query("plt_ro_g",      plt_ro_g     );
      pp.query("plt_MW_g",      plt_MW_g     );
      pp.query("plt_h_g",       plt_h_g      );
      pp.query("plt_T_g",       plt_T_g      );
      pp.query("plt_trac",      plt_trac     );
      pp.query("plt_cp_g",      plt_cp_g     );
      pp.query("plt_k_g",       plt_k_g      );
      pp.query("plt_mu_g",      plt_mu_g     );
      pp.query("plt_diveu",     plt_diveu    );
      pp.query("plt_vort",      plt_vort     );
      pp.query("plt_volfrac",   plt_volfrac  );
      pp.query("plt_gradp_g",   plt_gradp_g  );
      pp.query("plt_X_g",       plt_X_gk     );
      pp.query("plt_D_g",       plt_D_gk     );
      pp.query("plt_cp_gk",     plt_cp_gk    );
      pp.query("plt_h_gk",      plt_h_gk     );
      pp.query("plt_txfr",      plt_txfr     );
      pp.query("plt_chem_txfr", plt_chem_txfr);
      pp.query("plt_proc",      plt_proc     );
      pp.query("plt_proc_p",    plt_proc_p   );
      pp.query("plt_cost_p",    plt_cost_p   );

      // Special test for CCSE regression test. Override all individual
      // flags and save all data to plot file.

      int plt_ccse_regtest = 0;
      pp.query("plt_regtest", plt_ccse_regtest);

      if (plt_ccse_regtest != 0) {
        plt_vel_g     = 1;
        plt_ep_g      = 1;
        plt_p_g       = 0;
        plt_ro_g      = 1;
        plt_MW_g      = reactions.solve();
        plt_h_g       = 1;
        plt_T_g       = 1;
        plt_trac      = 1;
        plt_cp_g      = 1;
        plt_k_g       = 1;
        plt_mu_g      = 1;
        plt_vort      = 1;
        plt_diveu     = 1;
        plt_volfrac   = 1;
        plt_gradp_g   = 1;
        plt_X_gk      = fluid.solve_species();
        plt_D_gk      = fluid.solve_species();
        plt_cp_gk     = 0; //fluid.solve_species() )&&)&& fluid.solve_enthalpy();
        plt_h_gk      = 0; //fluid.solve_species() )&&)&& fluid.solve_enthalpy();
        plt_txfr      = 0;
        plt_chem_txfr = fluid.solve_species() && reactions.solve();
        plt_proc      = 0;
        plt_proc_p    = 0;
        plt_cost_p    = 0;
      }

      // Count the number of variables to save.
      if (plt_vel_g   == 1) pltVarCount += 3;
      if (plt_gradp_g == 1) pltVarCount += 3;
      if (plt_ep_g    == 1) pltVarCount += 1;
      if (plt_p_g     == 1) pltVarCount += 1;
      if (plt_ro_g    == 1) pltVarCount += 1;
      if (plt_MW_g    == 1) pltVarCount += 1;
      if (plt_trac    == 1) pltVarCount += 1;
      if (plt_mu_g    == 1) pltVarCount += 1;
      if (plt_vort    == 1) pltVarCount += 1;
      if (plt_diveu   == 1) pltVarCount += 1;
      if (plt_volfrac == 1) pltVarCount += 1;
      if (plt_proc    == 1) pltVarCount += 1;
      if (plt_proc_p  == 1) pltVarCount += 1;
      if (plt_cost_p  == 1) pltVarCount += 1;

      if (fluid.solve_enthalpy()) {
        if (plt_T_g  == 1) pltVarCount += 1;
        if (plt_cp_g == 1) pltVarCount += 1;
        if (plt_k_g  == 1) pltVarCount += 1;
        if (plt_h_g  == 1) pltVarCount += 1;
      }

      if (fluid.solve_species()) {
        if (plt_X_gk == 1)  pltVarCount += fluid.nspecies();
        if (plt_D_gk == 1)  pltVarCount += fluid.nspecies();

        if (fluid.solve_enthalpy()) {
          if (plt_cp_gk == 1) pltVarCount += fluid.nspecies();
          if (plt_h_gk == 1)  pltVarCount += fluid.nspecies();
        }
      }

      Transfer txfr_idxs(fluid.nspecies(), reactions.nreactions());

      if (m_dem.solve() || m_pic.solve()) {
        if (plt_txfr == 1) pltVarCount += txfr_idxs.countNoChem;
      }

      if (fluid.solve_species() && reactions.solve()) {
        ChemTransfer chem_txfr_idxs(fluid.nspecies(), reactions.nreactions());
        if (plt_chem_txfr == 1) pltVarCount += chem_txfr_idxs.count;
      }
    }

    if (m_dem.solve() || m_pic.solve()) {

      GetSolidsIOPltFlags(pp, write_real_comp, write_int_comp);

    }
}


void
MfixRW::GetSolidsIOPltFlags (ParmParse& pp,
                             Vector<int>& write_real_comp_out,
                             Vector<int>& write_int_comp_out)
{
  int plt_ccse_regtest = 0;
  pp.query("plt_regtest", plt_ccse_regtest);

  const runtimeRealData rtData(solids.nspecies()*solids.solve_species(),
                               fluid.nspecies()*fluid.solve_species(),
                               reactions.nreactions()*reactions.solve());

  // Runtime-added variables
  const int size = SoArealData::count + rtData.count;
  write_real_comp_out.resize(size, 1);

  // All flags are true by default so we only need to turn off the
  // variables we don't want if not doing CCSE regression tests.
  if (plt_ccse_regtest == 0) {

    int input_value = 0;
    pp.query("plt_radius", input_value);
    write_real_comp_out[SoArealData::radius] = input_value;

    input_value = 0;
    pp.query("plt_volume", input_value);
    write_real_comp_out[SoArealData::volume] = input_value;

    input_value = 0;
    pp.query("plt_mass", input_value);
    write_real_comp_out[SoArealData::mass] = input_value;

    input_value = 0;
    pp.query("plt_ro_p", input_value);
    write_real_comp_out[SoArealData::density] = input_value;

    input_value = 0;
    pp.query("plt_omoi", input_value);
    write_real_comp_out[SoArealData::oneOverI] = input_value;

    input_value = 1;
    pp.query("plt_vel_p", input_value);
    write_real_comp_out[SoArealData::velx] = input_value;
    write_real_comp_out[SoArealData::vely] = input_value;
    write_real_comp_out[SoArealData::velz] = input_value;

    input_value = 0;
    pp.query("plt_omega_p", input_value);
    write_real_comp_out[SoArealData::omegax] = input_value;
    write_real_comp_out[SoArealData::omegay] = input_value;
    write_real_comp_out[SoArealData::omegaz] = input_value;

    input_value = 0;
    pp.query("plt_statwt", input_value);
    write_real_comp_out[SoArealData::statwt] = input_value;

    input_value = 0;
    pp.query("plt_drag_p", input_value);
    write_real_comp_out[SoArealData::dragcoeff] = input_value;  // drag coeff
    write_real_comp_out[SoArealData::dragx] = input_value;  // dragx
    write_real_comp_out[SoArealData::dragy] = input_value;  // dragy
    write_real_comp_out[SoArealData::dragz] = input_value;  // dragz

    input_value = 0;
    pp.query("plt_cp_s", input_value);
    write_real_comp_out[SoArealData::cp_s] = input_value;  // specific heat

    input_value = 0;
    pp.query("plt_T_p", input_value);
    write_real_comp_out[SoArealData::temperature] = input_value;  // temperature

    input_value = 0;
    pp.query("plt_convection", input_value);
    write_real_comp_out[SoArealData::convection] = input_value;  // heat transfer coefficient

    int gap = SoArealData::count;

    if (solids.solve_species())
    {
      input_value = 0;
      pp.query("plt_X_s", input_value);

      const int start = gap + rtData.X_sn;
      for(int n(0); n < solids.nspecies(); ++n)
        write_real_comp_out[start+n] = input_value;
    }

    if (reactions.solve())
    {
      input_value = 0;
      pp.query("plt_vel_s_txfr", input_value);

      const int start = gap + rtData.vel_txfr;
      for(int n(0); n < 3; ++n)
        write_real_comp_out[start+n] = input_value;
    }

    if (reactions.solve())
    {
      input_value = 0;
      pp.query("plt_h_s_txfr", input_value);

      const int start = gap + rtData.h_txfr;
      write_real_comp_out[start] = input_value;
    }

    if (reactions.solve())
    {
      input_value = 0;
      pp.query("plt_mass_sn_txfr", input_value);

      const int start = gap + rtData.mass_txfr;
      for(int n(0); n < amrex::max(solids.nspecies(), fluid.nspecies()); ++n)
        // Since we allocated a number of components for mass transfer that
        // is equal to th max nb of species between solids and fluid, here
        // we check that n is smaller than solids species nb
        if (n < solids.nspecies())
          write_real_comp_out[start+n] = input_value;
        else // if n is larger, then we always set it to false
          write_real_comp_out[start+n] = 0;
    }

    // Int data
    input_value = 0;
    pp.query("plt_phase", input_value);
    write_int_comp_out[SoAintData::phase] = input_value;

    input_value = 0;
    pp.query("plt_state", input_value);
    write_int_comp_out[SoAintData::state] = input_value;
  }
}


void
MfixRW::WritePlotFile (std::string& plot_file_in, int nstep, Real time)
{
    // If we've already written this plotfile, don't do it again!
    if (nstep == last_plt) return;

    // Now set last_plt to nstep ...
    last_plt = nstep;

    BL_PROFILE("mfix::WritePlotFile()");

    const std::string& plotfilename = amrex::Concatenate(plot_file_in,nstep);

    amrex::Print() << "  Writing plotfile " << plotfilename <<  " at time " << time << std::endl;

    if (pltVarCount > 0) {

      const int ngrow = 0;

      Vector<std::string> pltFldNames;
      Vector< std::unique_ptr<MultiFab> > mf(nlev);

      // Velocity components
      if (plt_vel_g   == 1) {
        pltFldNames.push_back("u_g");
        pltFldNames.push_back("v_g");
        pltFldNames.push_back("w_g");
      }

      // Pressure gradient
      if (plt_gradp_g == 1) {
        pltFldNames.push_back("gpx");
        pltFldNames.push_back("gpy");
        pltFldNames.push_back("gpz");
      }

      // Fluid volume fraction
      if (plt_ep_g == 1)
        pltFldNames.push_back("ep_g");

      // Fluid pressure
      if (plt_p_g == 1)
        pltFldNames.push_back("p_g");

      // Fluid density
      if (plt_ro_g == 1)
        pltFldNames.push_back("ro_g");

      // Fluid molecular weight
      if (plt_MW_g == 1)
        pltFldNames.push_back("MW_g");

      // Fluid enthalpy
      if (fluid.solve_enthalpy() && plt_h_g == 1)
        pltFldNames.push_back("h_g");

      // Temperature in fluid
      if (fluid.solve_enthalpy() && plt_T_g == 1)
        pltFldNames.push_back("T_g");

      // Tracer in fluid
      if (plt_trac == 1)
        pltFldNames.push_back("trac");

      // Specific heat
      if (fluid.solve_enthalpy() && plt_cp_g == 1)
        pltFldNames.push_back("cp_g");

      // Thermal conductivity
      if (fluid.solve_enthalpy() && plt_k_g == 1)
        pltFldNames.push_back("k_g");

      // Fluid viscosity
      if (plt_mu_g == 1)
        pltFldNames.push_back("mu_g");

      // vorticity
      if (plt_vort == 1)
        pltFldNames.push_back("vort");

      // div(ep_g.u)
      if (plt_diveu == 1)
        pltFldNames.push_back("diveu");

      // EB cell volume fraction
      if (plt_volfrac == 1)
        pltFldNames.push_back("volfrac");

      // rank of fluid grids
      if (plt_proc == 1)
        pltFldNames.push_back("proc");

      // rank of particle grids
      if (plt_proc_p == 1)
        pltFldNames.push_back("proc_p");

      // cost of particle cell
      if (plt_cost_p == 1)
        pltFldNames.push_back("cost_p");

      // Fluid species mass fractions
      if (fluid.solve_species() && plt_X_gk == 1)
        for (std::string specie: fluid.species_names())
          pltFldNames.push_back("X_"+specie+"_g");

      // Fluid species mass diffusivities
      if (fluid.solve_species() && plt_D_gk == 1)
        for (std::string specie: fluid.species_names())
          pltFldNames.push_back("D_"+specie+"_g");

      // Fluid species specific heat
      if (fluid.solve_species() && fluid.solve_enthalpy() && plt_cp_gk == 1)
        for (std::string specie: fluid.species_names())
          pltFldNames.push_back("cp_"+specie+"_g");

      // Fluid species enthalpy
      if (fluid.solve_species() && fluid.solve_enthalpy() && plt_h_gk == 1)
        for (std::string specie: fluid.species_names())
          pltFldNames.push_back("h_"+specie+"_g");

      // Fluid species density reaction rates
      if (plt_txfr == 1) {
        pltFldNames.push_back("drag_x");
        pltFldNames.push_back("drag_y");
        pltFldNames.push_back("drag_z");
        pltFldNames.push_back("beta");
        pltFldNames.push_back("gammaTp");
        pltFldNames.push_back("gamma");
      }

      // Fluid species density reaction rates
      if (fluid.solve_species() && reactions.solve() && plt_chem_txfr == 1) {
        for(std::string specie: fluid.species_names())
          pltFldNames.push_back("chem_ro_txfr_"+specie);

        pltFldNames.push_back("chem_velx_txfr");
        pltFldNames.push_back("chem_vely_txfr");
        pltFldNames.push_back("chem_velz_txfr");
        pltFldNames.push_back("chem_h_txfr");
      }


      for (int lev = 0; lev < nlev; ++lev)
      {
        // Multifab to hold all the variables -- there can be only one!!!!
        const int ncomp = pltVarCount;
        mf[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], ncomp, ngrow,  MFInfo(), *ebfactory[lev]);

        int lc=0;

        // Velocity components
        if (plt_vel_g == 1) {
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->vel_g), 0, lc, 3, 0);
          lc += 3;
        }

        // Pressure gradient
        if (plt_gradp_g == 1) {
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->gp, 0, lc, 3, 0);
          lc += 3;
        }

        // Fluid volume fraction
        if (plt_ep_g == 1) {
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->ep_g, 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid pressure
        if (plt_p_g == 1) {
          MultiFab p_nd(m_leveldata[lev]->p_g->boxArray(), dmap[lev], 1, 0);
          p_nd.setVal(0.);
          MultiFab::Copy(p_nd, *m_leveldata[lev]->p_g, 0, 0, 1, 0);
          MultiFab::Add (p_nd, *m_leveldata[lev]->p0_g, 0, 0, 1, 0);
          amrex::average_node_to_cellcenter(*mf[lev], lc, p_nd, 0, 1);
          lc += 1;
        }

        // Fluid density
        if (plt_ro_g == 1) {
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->ro_g), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid molecular weight
        if (plt_MW_g == 1) {

          const int nspecies_g = fluid.nspecies();
          const auto& fluid_parms = fluid.parameters();
          const int fluid_is_a_mixture = fluid.isMixture();

          MultiFab& T_g = *(m_leveldata[lev]->T_g);

          MultiFab MW_g(T_g.boxArray(), T_g.DistributionMap(), T_g.nComp(),
                        T_g.nGrow(), MFInfo(), T_g.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(T_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real      > const& MW_g_array = MW_g.array(mfi);
            Array4<Real const> const& X_gk_array = fluid.isMixture() ?
              m_leveldata[lev]->X_gk->const_array(mfi) : Array4<const Real>();

            ParallelFor(bx, [MW_g_array,X_gk_array,nspecies_g,fluid_is_a_mixture,
                fluid_parms]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              if (fluid_is_a_mixture) {
                Real MW_g_loc(0);

                for (int n(0); n < nspecies_g; ++n) {
                  const Real MW_gk = fluid_parms.get_MW_gk<run_on>(n);

                  MW_g_loc += X_gk_array(i,j,k,n) / MW_gk;
                }

                MW_g_array(i,j,k) = 1. / MW_g_loc;
              } else {
                MW_g_array(i,j,k) = fluid_parms.get_MW_g<run_on>();
              }
            });
          }

          EB_set_covered(MW_g, 0, MW_g.nComp(), MW_g.nGrow(), covered_val);

          MultiFab::Copy(*mf[lev], MW_g, 0, lc, 1, 0);

          lc += 1;
        }

        // Fluid enthalpy
        if (fluid.solve_enthalpy() && plt_h_g == 1) {
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->h_g), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid temperature
        if (fluid.solve_enthalpy() && plt_T_g == 1) {
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->T_g), 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid tracer
        if (plt_trac == 1) {
          MultiFab::Copy(*mf[lev], (*m_leveldata[lev]->trac), 0, lc, 1, 0);
          lc += 1;
        }

        // Specific heat
        if (fluid.solve_enthalpy() && plt_cp_g == 1) {

          const auto& fluid_parms = fluid.parameters();
          int fluid_is_a_mixture = fluid.isMixture();
          int nspecies_g = fluid.nspecies();

          MultiFab& T_g = *(m_leveldata[lev]->T_g);

          MultiFab cp_g(T_g.boxArray(), T_g.DistributionMap(), T_g.nComp(),
                        T_g.nGrow(), MFInfo(), T_g.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(T_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real const> dummy_arr;

            Array4<Real      > const& cp_g_array = cp_g.array(mfi);
            Array4<Real const> const& T_g_array  = T_g.const_array(mfi);

            Array4<Real const> const& X_gk_array = fluid.isMixture() ?
              m_leveldata[lev]->X_gk->const_array(mfi) : dummy_arr;

            ParallelFor(bx, [cp_g_array,T_g_array,X_gk_array,fluid_parms,
                fluid_is_a_mixture,nspecies_g]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              const Real Tg = T_g_array(i,j,k);

              if (!fluid_is_a_mixture) {
                cp_g_array(i,j,k) = fluid_parms.calc_cp_g<run_on>(Tg);

              } else {
                Real cp_g_loc = 0;

                for (int n_g(0); n_g < nspecies_g; ++n_g) {
                  const Real cp_gk = fluid_parms.calc_cp_gk<run_on>(Tg,n_g);

                  cp_g_loc += X_gk_array(i,j,k,n_g)*cp_gk;
                }

                cp_g_array(i,j,k) = cp_g_loc;
              }
            });
          }

          EB_set_covered(cp_g, 0, cp_g.nComp(), cp_g.nGrow(), covered_val);

          MultiFab::Copy(*mf[lev], cp_g, 0, lc, 1, 0);
          lc += 1;
        }

        // Thermal conductivity
        if (fluid.solve_enthalpy() && plt_k_g == 1) {

          const auto& fluid_parms = fluid.parameters();

          MultiFab& T_g = *(m_leveldata[lev]->T_g);

          MultiFab k_g(T_g.boxArray(), T_g.DistributionMap(), T_g.nComp(),
                       T_g.nGrow(), MFInfo(), T_g.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(T_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real      > const& k_g_array = k_g.array(mfi);
            Array4<Real const> const& T_g_array = T_g.const_array(mfi);

            ParallelFor(bx, [k_g_array,T_g_array,fluid_parms]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              k_g_array(i,j,k) = fluid_parms.calc_k_g(T_g_array(i,j,k));
            });
          }

          EB_set_covered(k_g, 0, k_g.nComp(), k_g.nGrow(), covered_val);

          MultiFab::Copy(*mf[lev], k_g, 0, lc, 1, 0);
          lc += 1;
        }

        // Fluid viscosity
        if (plt_mu_g == 1) {

          MultiFab& ep_g = *(m_leveldata[lev]->ep_g);

          MultiFab mu_g(ep_g.boxArray(), ep_g.DistributionMap(), ep_g.nComp(),
                        ep_g.nGrow(), MFInfo(), ep_g.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(ep_g,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real      > const& mu_g_array = mu_g.array(mfi);
            Array4<Real const> const& T_g_array  = fluid.solve_enthalpy()?
              m_leveldata[lev]->T_g->const_array(mfi) : Array4<Real const>();

            const int solve_enthalpy = fluid.solve_enthalpy();
            const auto& fluid_parms = fluid.parameters();

            const Real mu_g0 = fluid.mu_g();

            ParallelFor(bx, [mu_g_array,T_g_array,solve_enthalpy,fluid_parms,mu_g0]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              if (solve_enthalpy)
                mu_g_array(i,j,k) = fluid_parms.calc_mu_g(T_g_array(i,j,k));
              else
                mu_g_array(i,j,k) = mu_g0;
            });
          }

          EB_set_covered(mu_g, 0, mu_g.nComp(), mu_g.nGrow(), covered_val);

          MultiFab::Copy(*mf[lev], mu_g, 0, lc, 1, 0);
          lc += 1;
        }

        // vorticity
        if (plt_vort == 1) {
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->vort, 0, lc, 1, 0);
          lc += 1;
        }

        // div(ep_g.u)
        if (plt_diveu == 1) {
          amrex::average_node_to_cellcenter(*mf[lev], lc, *m_leveldata[lev]->diveu, 0, 1);
          lc += 1;
        }

        // EB cell volume fraction
        if (plt_volfrac == 1) {
          if (ebfactory[lev]) {
            MultiFab::Copy(*mf[lev], ebfactory[lev]->getVolFrac(), 0, lc, 1, 0);
          } else {
            mf[lev]->setVal(1.0, lc, 1, 0);
          }
          lc += 1;
        }

        // rank of fluid grids
        if( plt_proc == 1 ) {
          MultiFab::Copy(*mf[lev], *fluid_proc[lev], 0, lc, 1, 0);
          lc += 1;
        }

        // rank of particle grids
        if ( plt_proc_p == 1 ) {
          mf[lev]->ParallelCopy(*particle_proc[lev], 0, lc, 1, 0, 0);
          lc += 1;
        }

        // cost of particle cell
        if (plt_cost_p == 1) {
          mf[lev]->ParallelCopy(*particle_cost[lev], 0, lc, 1, 0, 0);
          lc += 1;
        }

        // Fluid species mass fractions
        if (fluid.solve_species() && plt_X_gk == 1) {
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->X_gk, 0, lc, fluid.nspecies(), 0);
          lc += fluid.nspecies();
        }

        // Species mass fraction
        if (fluid.solve_species() && plt_D_gk == 1) {

          const auto& fluid_parms = fluid.parameters();
          const int nspecies_g = fluid.nspecies();

          MultiFab& X_gk = *(m_leveldata[lev]->X_gk);

          MultiFab D_gk(X_gk.boxArray(), X_gk.DistributionMap(), nspecies_g,
                        X_gk.nGrow(), MFInfo(), X_gk.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(X_gk,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real      > const& D_gk_array = D_gk.array(mfi);

            ParallelFor(bx, [D_gk_array,nspecies_g,fluid_parms]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              for (int n(0); n < nspecies_g; ++n) {
                D_gk_array(i,j,k,n) = fluid_parms.get_D_g();
              }
            });
          }

          EB_set_covered(D_gk, 0, D_gk.nComp(), D_gk.nGrow(), covered_val);

          MultiFab::Copy(*mf[lev], D_gk, 0, lc, D_gk.nComp(), 0);

          lc += D_gk.nComp();
        }

        // Fluid species specific heat
        if (fluid.solve_species() && fluid.solve_enthalpy() && plt_cp_gk == 1) {

          const auto& fluid_parms = fluid.parameters();
          const int nspecies_g = fluid.nspecies();

          MultiFab& X_gk = *(m_leveldata[lev]->X_gk);

          MultiFab cp_gk(X_gk.boxArray(), X_gk.DistributionMap(), X_gk.nComp(),
                         X_gk.nGrow(), MFInfo(), X_gk.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(X_gk,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            Box const& bx = mfi.tilebox();

            Array4<Real      > const& cp_gk_array = cp_gk.array(mfi);
            Array4<Real const> const& T_g_array  = fluid.solve_enthalpy()?
              m_leveldata[lev]->T_g->const_array(mfi) : Array4<const Real>();

            ParallelFor(bx, [cp_gk_array,T_g_array,nspecies_g,fluid_parms]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              for (int n(0); n < nspecies_g; ++n) {
                const Real T_g = T_g_array(i,j,k);

                cp_gk_array(i,j,k,n) = fluid_parms.calc_cp_gk<run_on>(T_g,n);
              }
            });
          }

          EB_set_covered(cp_gk, 0, cp_gk.nComp(), cp_gk.nGrow(), covered_val);

          MultiFab::Copy(*mf[lev], cp_gk, 0, lc, nspecies_g, 0);

          lc += nspecies_g;
        }

        // Fluid species enthalpy
        if (fluid.solve_species() && plt_h_gk == 1) {

          const auto& fluid_parms = fluid.parameters();
          const int nspecies_g = fluid.nspecies();

          auto& ld = *m_leveldata[lev];

          MultiFab& X_gk = *(ld.X_gk);

          MultiFab h_gk(X_gk.boxArray(), X_gk.DistributionMap(), X_gk.nComp(),
                        X_gk.nGrow(), MFInfo(), X_gk.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(X_gk,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
            const EBFArrayBox& Xgk_fab = static_cast<EBFArrayBox const&>(X_gk[mfi]);
            const EBCellFlagFab& flags = Xgk_fab.getEBCellFlagFab();

            Box const& bx = mfi.tilebox();

            Array4<const Real> dummy_arr;

            Array4<      Real> const& h_gk_array = h_gk.array(mfi);
            Array4<const Real> const& T_g_array  = fluid.solve_enthalpy() ?
              (ld.T_g)->const_array(mfi) : dummy_arr;

            auto const& flags_arr = flags.const_array();

            ParallelFor(bx, [h_gk_array,T_g_array,nspecies_g,fluid_parms,flags_arr]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

              const Real Tg_loc = T_g_array(i,j,k);

              for (int n(0); n < nspecies_g; ++n)
                h_gk_array(i,j,k,n) = fluid_parms.calc_h_gk<run_on>(Tg_loc, n, cell_is_covered);
            });
          }

          EB_set_covered(h_gk, 0, h_gk.nComp(), h_gk.nGrow(), covered_val);

          MultiFab::Copy(*mf[lev], h_gk, 0, lc, nspecies_g, 0);

          lc += nspecies_g;
        }

        Transfer txfr_idxs(fluid.nspecies(), reactions.nreactions());

        if (plt_txfr == 1) {
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->txfr, 0, lc, txfr_idxs.countNoChem, 0);
          lc += txfr_idxs.countNoChem;
        }

        if (fluid.solve_species() && reactions.solve() && plt_chem_txfr == 1) {
          ChemTransfer chem_txfr_idxs(fluid.nspecies(), reactions.nreactions());
          MultiFab::Copy(*mf[lev], *m_leveldata[lev]->chem_txfr, 0, lc, chem_txfr_idxs.count, 0);
          lc += chem_txfr_idxs.count;
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
                                     geom, time, istep, ref_ratio);


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
        mf[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], 1, 0);

      Vector<const MultiFab*> mf2(finest_level+1);

      for (int lev = 0; lev <= finest_level; ++lev)
        mf2[lev] = mf[lev].get();

      // Write only the Headers corresponding to the "empty" mf/mf2 MultiFabs
      Vector<int> istep;
      istep.resize(nlev,nstep);
      amrex::WriteMultiLevelPlotfileHeaders(plotfilename, finest_level+1, mf2, names,
                                            geom, time, istep, ref_ratio);

    }

    WriteJobInfo(plotfilename);

    if (m_dem.solve() || m_pic.solve())
    {
        pc->WritePlotFile(plotfilename, "particles", write_real_comp,
                          write_int_comp, real_comp_names, int_comp_names);
    }
}


void
MfixRW::WriteSolidsPlotFile (SolidsPlotRegion& plot_region,
                             std::string& plot_file_in,
                             int nstep,
                             Real time)
{
  if ((m_dem.solve() || m_pic.solve()) && (solids_plot_regions() == true)) {

    // If we've already written this plotfile, don't do it again!
    if (nstep == plot_region.m_last_solids_plt) return;

    // Now set last_solids_plt to nstep ...
    plot_region.m_last_solids_plt = nstep;

    const RealBox region_extents = plot_region.m_region_extents;
    const std::string& region_name = plot_region.m_region_name;

    std::string solidsfilename = amrex::Concatenate(plot_file_in,nstep);
    solidsfilename += "_" + region_name;
    amrex::Print() << "  Writing solids plotfile " << solidsfilename <<  " at time " << time << std::endl;

    const int plot_types_nb = plot_region.m_h_plot_types.size();
    int* plot_types_ptr = plot_region.m_d_plot_types.dataPtr();

    auto F = [region_extents,plot_types_ptr,plot_types_nb]
      AMREX_GPU_DEVICE (const MFIXParticleContainer::SuperParticleType& p,
                        const amrex::RandomEngine&) noexcept -> bool
    {
      int particle_in_region = static_cast<int>(region_extents.contains(p.pos()));

      int type_found = static_cast<int>(plot_types_nb == 0);
      for (int n(0); n < plot_types_nb; ++n) {
        if (p.idata(SoAintData::phase) == plot_types_ptr[n]) {
          type_found = 1;
          break;
        }
      }

      return particle_in_region && type_found;
    };

    // no fluid
    {
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
        mf[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], 1, 0);

      Vector<const MultiFab*> mf2(finest_level+1);

      for (int lev = 0; lev <= finest_level; ++lev)
        mf2[lev] = mf[lev].get();

      // Write only the Headers corresponding to the "empty" mf/mf2 MultiFabs
      Vector<int> istep;
      istep.resize(nlev,nstep);
      amrex::WriteMultiLevelPlotfileHeaders(solidsfilename, finest_level+1, mf2, names,
                                            geom, time, istep, ref_ratio);
    }

    WriteJobInfo(solidsfilename);

    pc->WritePlotFile(solidsfilename, "particles", plot_region.m_write_real_comp,
                      plot_region.m_write_int_comp, real_comp_names, int_comp_names, F);
  }
}


void
MfixRW::WriteStaticPlotFile (const std::string & plotfilename) const
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
        mf[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], ncomp, ngrow,
                                             MFInfo(), *particle_ebfactory[lev]);

        // Don't iterate over all ncomp => last component is for volfrac
        for (int dcomp = 0; dcomp < ncomp - 1; dcomp++)
        {
            const BoxArray nd_ba = amrex::convert(grids[lev], IntVect::TheNodeVector());
            MultiFab mf_loc = MFUtil::regrid(nd_ba, dmap[lev], *static_vars[dcomp][lev], true);
            amrex::average_node_to_cellcenter(* mf[lev], dcomp, mf_loc, 0, 1, ngrow);
        }

        if (ebfactory[lev]) {
            EBFArrayBoxFactory ebf(* eb_levels[lev], geom[lev], grids[lev], dmap[lev],
                                   {nghost_eb_basic, nghost_eb_volume,
                                    nghost_eb_full}, m_eb_support_level);

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
                                   geom, time, istep, ref_ratio);

    WriteJobInfo(plotfilename);

    Print() << "  Done writing static quantities " << plotfilename << std::endl;
}

} // end namespace MfixIO
