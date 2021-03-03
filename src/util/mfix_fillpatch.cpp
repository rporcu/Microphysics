#include <mfix.H>
#include <mfix_fluid_parms.H>
#include <mfix_fillpatch_bc.H>

#include <AMReX_FillPatchUtil.H>

namespace
{
  mfix* mfix_for_fillpatching;
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
void set_ptr_to_mfix (mfix& mfix_for_fillpatching_in)
{
   mfix_for_fillpatching = &mfix_for_fillpatching_in;
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
inline
void VelFillBox (Box const& /*bx*/,
                 Array4<amrex::Real> const& dest,
                 const int dcomp,
                 const int numcomp,
                 GeometryData const& geom,
                 const Real time_in,
                 const BCRec* /*bcr*/,
                 const int /*bcomp*/,
                 const int /*orig_comp*/)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in VelFillBox");
    if (numcomp != 3)
         amrex::Abort("Must have numcomp = 3 in VelFillBox");

    const Box& domain = geom.Domain();

    // This is a bit hack-y but does get us the right level
    int lev = 0;
    for (int ilev = 0; ilev < 10; ilev++)
    {
//     const Geometry& lev_geom = mfix_for_fillpatching->GetParGDB()->Geom(ilev);
       const Geometry& lev_geom = mfix_for_fillpatching->get_geom_ref(ilev);
       if (domain.length()[0] == (lev_geom.Domain()).length()[0])
       {
         lev = ilev;
         break;
       }
    }

    // We are hard-wiring this fillpatch routine to define the Dirichlet values
    //    at the faces (not the ghost cell center)
    int extrap_dir_bcs = 0;

    // We only do this to make it not const
    Real time = time_in;

    FArrayBox dest_fab(dest);
    Elixir eli_dest_fab = dest_fab.elixir();

    mfix_for_fillpatching->set_velocity_bcs(time, lev, dest_fab, domain, &extrap_dir_bcs);
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
inline
void ScalarFillBox (Box const& /*bx*/,
                    Array4<amrex::Real> const& dest,
                    const int dcomp,
                    const int numcomp,
                    GeometryData const& geom,
                    const Real time_in,
                    const BCRec* /*bcr*/,
                    const int /*bcomp*/,
                    const int orig_comp)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in ScalarFillBox");
    if (numcomp != 1)
         amrex::Abort("Must have numcomp = 1 in ScalarFillBox");

    const Box& domain = geom.Domain();

    // This is a bit hack-y but does get us the right level
    int lev = 0;
    for (int ilev = 0; ilev < 10; ilev++)
    {
       const Geometry& lev_geom = mfix_for_fillpatching->GetParGDB()->Geom(ilev);
       if (domain.length()[0] == (lev_geom.Domain()).length()[0])
       {
         lev = ilev;
         break;
       }
    }

    // We only do this to make it not const
    Real time = time_in;

    FArrayBox dest_fab(dest);
    Elixir eli_dest_fab = dest_fab.elixir();

   if( orig_comp == 0 )
      mfix_for_fillpatching->set_density_bcs(time, lev, dest_fab, domain);
   else if(orig_comp == 1)
      mfix_for_fillpatching->set_tracer_bcs(time, lev, dest_fab, domain);
   else if(orig_comp == 5 and FLUID::solve_enthalpy)
      mfix_for_fillpatching->set_enthalpy_bcs(time, lev, dest_fab, domain);
   else
      amrex::Abort("Unknown component in ScalarFillBox!");

}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
inline
void SpeciesFillBox (Box const& /*bx*/,
                     Array4<amrex::Real> const& dest,
                     const int dcomp,
                     const int numcomp,
                     GeometryData const& geom,
                     const Real time_in,
                     const BCRec* /*bcr*/,
                     const int /*bcomp*/,
                     const int orig_comp)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in SpeciesFillBox");
    if (numcomp != FLUID::nspecies)
         amrex::Abort("Must have numcomp = nspecies_g in SpeciesFillBox");

    const Box& domain = geom.Domain();

    // This is a bit hack-y but does get us the right level
    int lev = 0;
    for (int ilev = 0; ilev < 10; ilev++)
    {
       const Geometry& lev_geom = mfix_for_fillpatching->GetParGDB()->Geom(ilev);
       if (domain.length()[0] == (lev_geom.Domain()).length()[0])
       {
         lev = ilev;
         break;
       }
    }

    // We only do this to make it not const
    Real time = time_in;

    FArrayBox dest_fab(dest);
    Elixir eli_dest_fab = dest_fab.elixir();

   if( orig_comp == 0 )
      mfix_for_fillpatching->set_mass_fractions_g_bcs(time, lev, dest_fab, domain);
   else if(orig_comp == 1)
      mfix_for_fillpatching->set_species_diffusivities_g_bcs(time, lev, dest_fab, domain);
   else
      amrex::Abort("Unknown component in ScalarFillBox!");

}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
mfix::FillPatchVel (int lev,
                    Real time,
                    MultiFab& mf,
                    int icomp,
                    int ncomp,
                    const Vector<BCRec>& bcs)
{
    // Hack so that ghost cells are not undefined
    mf.setVal(covered_val);

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetDataVel(0, time, smf, stime);

        CpuBndryFuncFab bfunc(VelFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                    geom[lev], physbc, 0);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataVel(lev-1, time, cmf, ctime);
        GetDataVel(lev  , time, fmf, ftime);

        CpuBndryFuncFab bfunc(VelFillBox);
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, 0);
    }
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
// NOTE: icomp here refers to whether we are filling 0: density, 1: tracer, 2: ep_g, 3: mu_g, 4: temperature, 5: enthalpy
void
mfix::FillPatchScalar (int lev,
                       Real time,
                       MultiFab& mf,
                       int icomp,
                       int ncomp,
                       const Vector<BCRec>& bcs)
{
    // Hack so that ghost cells are not undefined
    mf.setVal(covered_val);

    // icomp tells us which scalar we are fill-patching
    // But we send "0) into FillPatch since each scalar is stored in its own array

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;

        GetDataScalar(0, time, smf, icomp, stime);

        CpuBndryFuncFab bfunc(ScalarFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp,
                                    geom[lev], physbc, icomp);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataScalar(lev-1, time, cmf, icomp, ctime);
        GetDataScalar(lev  , time, fmf, icomp, ftime);

        CpuBndryFuncFab bfunc(ScalarFillBox);
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, 0, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, icomp);
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
                        const Vector<BCRec>& bcs)
{
    // Hack so that ghost cells are not undefined
    mf.setVal(covered_val);

    // icomp tells us which scalar we are fill-patching
    // But we send "0) into FillPatch since each scalar is stored in its own array

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;

        GetDataSpecies(0, time, smf, icomp, stime);

        CpuBndryFuncFab bfunc(SpeciesFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp,
                                    geom[lev], physbc, icomp);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataSpecies(lev-1, time, cmf, icomp, ctime);
        GetDataSpecies(lev  , time, fmf, icomp, ftime);

        CpuBndryFuncFab bfunc(SpeciesFillBox);
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, 0, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, icomp);
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
                     int icomp,
                     Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (icomp == 3)
       data.push_back(m_leveldata[lev]->mu_g);
    // TODO cp_g, k_g

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        if (icomp == 0) {
           data.push_back(m_leveldata[lev]->ro_g);
        } else if (icomp == 1) {
           data.push_back(m_leveldata[lev]->trac);
        } else if (icomp == 2) {
           data.push_back(m_leveldata[lev]->ep_g);
        } else if (icomp == 4) {
           data.push_back(m_leveldata[lev]->T_g);
        } else if (icomp == 5) {
           data.push_back(m_leveldata[lev]->h_g);
        }
        datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        if (icomp == 0) {
           data.push_back(m_leveldata[lev]->ro_go);
        } else if (icomp == 1) {
           data.push_back(m_leveldata[lev]->trac_o);
        } else if (icomp == 2) {
           data.push_back(m_leveldata[lev]->ep_go);
        } else if (icomp == 4) {
           data.push_back(m_leveldata[lev]->T_g);
        } else if (icomp == 5) {
           data.push_back(m_leveldata[lev]->h_go);
        }
        datatime.push_back(t_old[lev]);
    }
    else
    {
        if (icomp == 0) {
           data.push_back(m_leveldata[lev]->ro_go);
           data.push_back(m_leveldata[lev]->ro_g);
        } else if (icomp == 1) {
           data.push_back(m_leveldata[lev]->trac_o);
           data.push_back(m_leveldata[lev]->trac);
        } else if (icomp == 2) {
           data.push_back(m_leveldata[lev]->ep_go);
           data.push_back(m_leveldata[lev]->ep_g);
        } else if (icomp == 4) {
           data.push_back(m_leveldata[lev]->T_go);
           data.push_back(m_leveldata[lev]->T_g);
        } else if (icomp == 5) {
           data.push_back(m_leveldata[lev]->h_go);
           data.push_back(m_leveldata[lev]->h_g);
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

  amrex::Vector<amrex::BCRec> bcrec_force;
  bcrec_force.resize(ncomp);

  int l_probtype = -1;
  int lev = 0;

  {
    PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > physbc
          (geom[lev], bcrec_force, MFIXForFill{l_probtype});
        FillPatchSingleLevel(*force[lev], IntVect(ng), time,
                             {force[lev]}, {time},
                             0, 0, ncomp, geom[lev],
                             physbc, 0);
    }
    for (lev = 1; lev <= finest_level; ++lev)
    {
        PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > cphysbc
            (geom[lev-1], bcrec_force, MFIXForFill{l_probtype});
        PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > fphysbc
            (geom[lev  ], bcrec_force, MFIXForFill{l_probtype});
        Interpolater* mapper = &pc_interp;
        FillPatchTwoLevels(*force[lev], IntVect(ng), time,
                           {force[lev-1]}, {time},
                           {force[lev  ]}, {time},
                           0, 0, ncomp, geom[lev-1], geom[lev],
                           cphysbc, 0, fphysbc, 0,
                           refRatio(lev-1), mapper, bcrec_force, 0);
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

  const int l_nspecies = FLUID::nspecies;

  // First do FillPatch of {velocity, density, tracer, enthalpy} so we know
  // the ghost cells of these arrays are all filled
  for (int lev = 0; lev < nlev; lev++) {

    int state_comp, num_comp;

    // State with ghost cells
    MultiFab Sborder_u(grids[lev], dmap[lev], vel_in[lev]->nComp(), nghost_state(),
                       MFInfo(), *ebfactory[lev]);
    FillPatchVel(lev, time, Sborder_u, 0, Sborder_u.nComp(), bcs_u);

    // Copy each FAB back from Sborder_u into the vel array, complete with filled ghost cells
    MultiFab::Copy(*vel_in[lev], Sborder_u, 0, 0, vel_in[lev]->nComp(), vel_in[lev]->nGrow());

    MultiFab Sborder_s(grids[lev], dmap[lev], 1, nghost_state(), MFInfo(), *ebfactory[lev]);

    // We FillPatch density even if not advecting it because we need it in the projections
    state_comp =  0; // comp = 0 --> density
    num_comp = 1;
    FillPatchScalar(lev, time, Sborder_s, state_comp, num_comp, bcs_s);
    MultiFab::Copy(*ro_g_in[lev], Sborder_s, 0, 0, num_comp, ro_g_in[lev]->nGrow());

    if (advect_tracer) {
      state_comp =  1; // comp = 1 --> tracer
      num_comp = 1;
      FillPatchScalar(lev, time, Sborder_s, state_comp, num_comp, bcs_s);
      MultiFab::Copy(*trac_in[lev], Sborder_s, 0, 0, num_comp, trac_in[lev]->nGrow());
    }

    if (advect_enthalpy) {
      state_comp =  5; // comp = 1 --> enthalpy
      num_comp = 1;
      FillPatchScalar(lev, time, Sborder_s, state_comp, num_comp, bcs_s);
      MultiFab::Copy(*h_g_in[lev], Sborder_s, 0, 0, num_comp, h_g_in[lev]->nGrow());
    }

    if (advect_fluid_species) {
      MultiFab Sborder_X(grids[lev], dmap[lev], FLUID::nspecies, nghost_state(),
                         MFInfo(), *ebfactory[lev]);
      Sborder_X.setVal(0);
      state_comp = 0;
      num_comp = l_nspecies;
      FillPatchSpecies(lev, time, Sborder_X, state_comp, num_comp, bcs_X);
      MultiFab::Copy(*X_gk_in[lev], Sborder_X, 0, 0, num_comp,
                     X_gk_in[lev]->nGrow());
    }
  }
}
