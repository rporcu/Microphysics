#include <mfix.H>
#include <mfix_fluid.H>
#include <mfix_fillpatch_bc.H>

#include <AMReX_FillPatchUtil.H>

using namespace amrex;

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
mfix::FillPatchVel (int lev,
                    Real time,
                    MultiFab& mf,
                    int icomp,
                    int ncomp,
                    const Vector<BCRec>& bcrec)
{
    // Hack so that ghost cells are not undefined
    mf.setVal(covered_val);

    set_velocity_bc_values(time);

    const int minf = BCList::minf;

    if (lev == 0) {

      Vector<MultiFab*> smf;
      Vector<Real> stime;
      GetDataVel(0, time, smf, stime);

      PhysBCFunct<GpuBndryFuncFab<MFIXVelFill> > physbc
        (geom[lev], bcrec, MFIXVelFill{minf,
            m_boundary_conditions.bc_u_g().data(),
            m_boundary_conditions.bc_v_g().data(),
            m_boundary_conditions.bc_w_g().data(),
            bc_list.bc_ilo[lev]->array(), bc_list.bc_ihi[lev]->array(),
            bc_list.bc_jlo[lev]->array(), bc_list.bc_jhi[lev]->array(),
            bc_list.bc_klo[lev]->array(), bc_list.bc_khi[lev]->array()
            });

      amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                  geom[lev], physbc, 0);
    }
    else
    {
      Vector<MultiFab*> cmf, fmf;
      Vector<Real> ctime, ftime;
      GetDataVel(lev-1, time, cmf, ctime);
      GetDataVel(lev  , time, fmf, ftime);

      PhysBCFunct<GpuBndryFuncFab<MFIXVelFill> > cphysbc
        (geom[lev-1], bcrec, MFIXVelFill{minf,
            m_boundary_conditions.bc_u_g().data(),
            m_boundary_conditions.bc_v_g().data(),
            m_boundary_conditions.bc_w_g().data(),
            bc_list.bc_ilo[lev-1]->array(), bc_list.bc_ihi[lev-1]->array(),
            bc_list.bc_jlo[lev-1]->array(), bc_list.bc_jhi[lev-1]->array(),
            bc_list.bc_klo[lev-1]->array(), bc_list.bc_khi[lev-1]->array()
            });

      PhysBCFunct<GpuBndryFuncFab<MFIXVelFill> > fphysbc
        (geom[lev], bcrec, MFIXVelFill{minf,
            m_boundary_conditions.bc_u_g().data(),
            m_boundary_conditions.bc_v_g().data(),
            m_boundary_conditions.bc_w_g().data(),
            bc_list.bc_ilo[lev]->array(), bc_list.bc_ihi[lev]->array(),
            bc_list.bc_jlo[lev]->array(), bc_list.bc_jhi[lev]->array(),
            bc_list.bc_klo[lev]->array(), bc_list.bc_khi[lev]->array()
            });

      Interpolater* mapper = &cell_cons_interp;

      amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                0, icomp, ncomp, geom[lev-1], geom[lev],
                                cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcrec, 0);
    }
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
mfix::FillPatchScalar (int lev, Real time, MultiFab& mf,
                       ScalarToFill scalar_id,
                       const Real* bc_scalar,
                       const Vector<BCRec>& bcrec)
{
    // Hack so that ghost cells are not undefined
    mf.setVal(covered_val);

    const int minf = BCList::minf;

    const int icomp = 0;
    const int ncomp = 1;

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;

        GetDataScalar(0, time, smf, scalar_id, stime);

        PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > physbc
          (geom[lev], bcrec, MFIXScalarFill{minf, bc_scalar,
              bc_list.bc_ilo[lev]->array(), bc_list.bc_ihi[lev]->array(),
              bc_list.bc_jlo[lev]->array(), bc_list.bc_jhi[lev]->array(),
              bc_list.bc_klo[lev]->array(), bc_list.bc_khi[lev]->array()
              });

        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp,
                                    geom[lev], physbc, icomp);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataScalar(lev-1, time, cmf, scalar_id, ctime);
        GetDataScalar(lev  , time, fmf, scalar_id, ftime);

        PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > cphysbc
          (geom[lev-1], bcrec, MFIXScalarFill{minf, bc_scalar,
              bc_list.bc_ilo[lev-1]->array(), bc_list.bc_ihi[lev-1]->array(),
              bc_list.bc_jlo[lev-1]->array(), bc_list.bc_jhi[lev-1]->array(),
              bc_list.bc_klo[lev-1]->array(), bc_list.bc_khi[lev-1]->array()
              });

        PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > fphysbc
          (geom[lev], bcrec, MFIXScalarFill{minf, bc_scalar,
              bc_list.bc_ilo[lev]->array(), bc_list.bc_ihi[lev]->array(),
              bc_list.bc_jlo[lev]->array(), bc_list.bc_jhi[lev]->array(),
              bc_list.bc_klo[lev]->array(), bc_list.bc_khi[lev]->array()
              });
        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, 0, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcrec, icomp);
    }
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
// NOTE: icomp here refers to whether we are filling 0: fluid species
void
mfix::FillPatchSpecies (int lev,
                        Real time,
                        MultiFab& mf,
                        int icomp,
                        int ncomp,
                        const Vector<BCRec>& bcrec)
{
    // Hack so that ghost cells are not undefined
    mf.setVal(covered_val);

    const int minf = BCList::minf;

    // icomp tells us which scalar we are fill-patching
    // But we send "0) into FillPatch since each scalar is stored in its own array

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;

        GetDataSpecies(0, time, smf, icomp, stime);

        PhysBCFunct<GpuBndryFuncFab<MFIXSpeciesFill> > physbc
          (geom[lev], bcrec, MFIXSpeciesFill{minf,
              m_boundary_conditions.bc_X_gk_ptr().data(),
              bc_list.bc_ilo[lev]->array(), bc_list.bc_ihi[lev]->array(),
              bc_list.bc_jlo[lev]->array(), bc_list.bc_jhi[lev]->array(),
              bc_list.bc_klo[lev]->array(), bc_list.bc_khi[lev]->array()
              });

        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp,
                                    geom[lev], physbc, icomp);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataSpecies(lev-1, time, cmf, icomp, ctime);
        GetDataSpecies(lev  , time, fmf, icomp, ftime);

        PhysBCFunct<GpuBndryFuncFab<MFIXSpeciesFill> > cphysbc
          (geom[lev-1], bcrec, MFIXSpeciesFill{minf,
              m_boundary_conditions.bc_X_gk_ptr().data(),
              bc_list.bc_ilo[lev-1]->array(), bc_list.bc_ihi[lev-1]->array(),
              bc_list.bc_jlo[lev-1]->array(), bc_list.bc_jhi[lev-1]->array(),
              bc_list.bc_klo[lev-1]->array(), bc_list.bc_khi[lev-1]->array()
              });

        PhysBCFunct<GpuBndryFuncFab<MFIXSpeciesFill> > fphysbc
          (geom[lev], bcrec, MFIXSpeciesFill{minf,
              m_boundary_conditions.bc_X_gk_ptr().data(),
              bc_list.bc_ilo[lev]->array(), bc_list.bc_ihi[lev]->array(),
              bc_list.bc_jlo[lev]->array(), bc_list.bc_jhi[lev]->array(),
              bc_list.bc_klo[lev]->array(), bc_list.bc_khi[lev]->array()
              });

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, 0, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcrec, icomp);
    }
}

// Utility to copy in data from phi_old and/or phi_new into another multifab
void
mfix::GetDataVel (int lev,
                  Real time,
                  Vector<MultiFab*>& data,
                  Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        data.push_back(m_leveldata[lev]->vel_g);
        datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        data.push_back(m_leveldata[lev]->vel_go);
        datatime.push_back(t_old[lev]);
    }
    else
    {
        data.push_back(m_leveldata[lev]->vel_go);
        data.push_back(m_leveldata[lev]->vel_g);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

// Utility to copy in data from phi_old and/or phi_new into another multifab
void
mfix::GetDataScalar (int lev,
                     Real time,
                     Vector<MultiFab*>& data,
                     ScalarToFill scalar_id,
                     Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps) {
      switch (scalar_id) {
        case (ScalarToFill::Density ): {
          data.push_back(m_leveldata[lev]->ro_g);
          break;
        }
        case (ScalarToFill::Tracer): {
          data.push_back(m_leveldata[lev]->trac);
          break;
        }
        case (ScalarToFill::Enthalpy): {
          data.push_back(m_leveldata[lev]->h_g);
          break;
        }
        case (ScalarToFill::VolFrac): {
          data.push_back(m_leveldata[lev]->ep_g);
          break;
        }
        case (ScalarToFill::Temperature): {
          data.push_back(m_leveldata[lev]->T_g);
          break;
        }
      }
      datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps) {

      switch (scalar_id) {
        case (ScalarToFill::Density ): {
          data.push_back(m_leveldata[lev]->ro_go);
          break;
        }
        case (ScalarToFill::Tracer): {
          data.push_back(m_leveldata[lev]->trac_o);
          break;
        }
        case (ScalarToFill::Enthalpy): {
          data.push_back(m_leveldata[lev]->h_go);
          break;
        }
        case (ScalarToFill::VolFrac): {
          data.push_back(m_leveldata[lev]->ep_g);
          break;
        }
        case (ScalarToFill::Temperature): {
          data.push_back(m_leveldata[lev]->T_go);
          break;
        }
      }
      datatime.push_back(t_old[lev]);
    }
    else
    {
      switch (scalar_id) {
        case (ScalarToFill::Density ): {
          data.push_back(m_leveldata[lev]->ro_go);
          data.push_back(m_leveldata[lev]->ro_g);
          break;
        }
        case (ScalarToFill::Tracer): {
          data.push_back(m_leveldata[lev]->trac_o);
          data.push_back(m_leveldata[lev]->trac);
          break;
        }
        case (ScalarToFill::Enthalpy): {
          data.push_back(m_leveldata[lev]->h_go);
          data.push_back(m_leveldata[lev]->h_g);
          break;
        }
        case (ScalarToFill::VolFrac): {
          data.push_back(m_leveldata[lev]->ep_g);
          data.push_back(m_leveldata[lev]->ep_g);
          break;
        }
        case (ScalarToFill::Temperature): {
          data.push_back(m_leveldata[lev]->T_go);
          data.push_back(m_leveldata[lev]->T_g);
          break;
        }
      }
      datatime.push_back(t_old[lev]);
      datatime.push_back(t_new[lev]);
    }
}

// Utility to copy in data from phi_old and/or phi_new into another multifab
void
mfix::GetDataSpecies (int lev,
                      Real time,
                      Vector<MultiFab*>& data,
                      int icomp,
                      Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        if (icomp == 0) {
           data.push_back(m_leveldata[lev]->X_gk);
        }
        datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        if (icomp == 0) {
           data.push_back(m_leveldata[lev]->X_gko);
        }
        datatime.push_back(t_old[lev]);
    }
    else
    {
        if (icomp == 0) {
           data.push_back(m_leveldata[lev]->X_gko);
           data.push_back(m_leveldata[lev]->X_gk);
        }
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}


void mfix::fillpatch_force (Real time, Vector<MultiFab*> const& force, int ng)
{

  /* JMusser: This was some serious copy-and-paste straight out of      *
   * incflow. If this doesn't do what it is supposed to do, it's not my *
   * fault, because I have no idea what any of this means.              *
   * You're Welcome.                                                    */
  const int ncomp = force[0]->nComp();
  const auto& bcrec = get_force_bcrec();

  int l_probtype = -1;
  int lev = 0;

  {
    PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > physbc
          (geom[lev], bcrec, MFIXForFill{l_probtype});

    FillPatchSingleLevel(*force[lev], IntVect(ng), time,
                         {force[lev]}, {time},
                         0, 0, ncomp, geom[lev],
                         physbc, 0);
    }
    for (lev = 1; lev <= finest_level; ++lev)
    {
        PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > cphysbc
            (geom[lev-1], bcrec, MFIXForFill{l_probtype});
        PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > fphysbc
            (geom[lev  ], bcrec, MFIXForFill{l_probtype});
        Interpolater* mapper = &pc_interp;
        FillPatchTwoLevels(*force[lev], IntVect(ng), time,
                           {force[lev-1]}, {time},
                           {force[lev  ]}, {time},
                           0, 0, ncomp, geom[lev-1], geom[lev],
                           cphysbc, 0, fphysbc, 0,
                           refRatio(lev-1), mapper, bcrec, 0);
    }
}

void
mfix::fillpatch_all (Vector< MultiFab* > const& vel_in,
                     Vector< MultiFab* > const& ro_g_in,
                     Vector< MultiFab* > const& h_g_in,
                     Vector< MultiFab* > const& trac_in,
                     Vector< MultiFab* > const& X_gk_in,
                     Real time)
{

  // First do FillPatch of {velocity, density, tracer, enthalpy} so we know
  // the ghost cells of these arrays are all filled
  for (int lev = 0; lev < nlev; lev++) {

    // State with ghost cells
    MultiFab Sborder_u(grids[lev], dmap[lev], vel_in[lev]->nComp(), nghost_state(),
                       MFInfo(), *ebfactory[lev]);
    FillPatchVel(lev, time, Sborder_u, 0, Sborder_u.nComp(), get_velocity_bcrec());

    // Copy each FAB back from Sborder_u into the vel array, complete with filled ghost cells
    MultiFab::Copy(*vel_in[lev], Sborder_u, 0, 0, vel_in[lev]->nComp(), vel_in[lev]->nGrow());


    MultiFab Sborder_s(grids[lev], dmap[lev], 1, nghost_state(), MFInfo(), *ebfactory[lev]);


    // We FillPatch density even if not advecting it because we need it in the projections
    FillPatchScalar(lev, time, Sborder_s, ScalarToFill::Density,
                    m_boundary_conditions.bc_ro_g().data(), get_density_bcrec());
    MultiFab::Copy(*ro_g_in[lev], Sborder_s, 0, 0, 1, ro_g_in[lev]->nGrow());


    if (fluid.solve_tracer()) {
      FillPatchScalar(lev, time, Sborder_s, ScalarToFill::Tracer,
                      m_boundary_conditions.bc_tracer().data(), get_tracer_bcrec());
      MultiFab::Copy(*trac_in[lev], Sborder_s, 0, 0, 1, trac_in[lev]->nGrow());
    }

    if (fluid.solve_enthalpy()) {
      FillPatchScalar(lev, time, Sborder_s, ScalarToFill::Enthalpy,
                      m_boundary_conditions.bc_h_g().data(), get_enthalpy_bcrec());
      MultiFab::Copy(*h_g_in[lev], Sborder_s, 0, 0, 1, h_g_in[lev]->nGrow());
    }

    if (fluid.solve_species()) {
      MultiFab Sborder_X(grids[lev], dmap[lev], X_gk_in[lev]->nComp(),
                         nghost_state(), MFInfo(), *ebfactory[lev]);
      Sborder_X.setVal(0);

      FillPatchSpecies(lev, time, Sborder_X, 0, Sborder_X.nComp(), get_species_bcrec());
      MultiFab::Copy(*X_gk_in[lev], Sborder_X, 0, 0,
                     X_gk_in[lev]->nComp(), X_gk_in[lev]->nGrow());
    }
  }
}
