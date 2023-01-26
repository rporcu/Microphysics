#include <mfix_rw.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_monitors.H>
#include <mfix_utils.H>

#include <AMReX_FabArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiCutFab.H>

using namespace amrex;


namespace MfixIO {

void
MfixRW::WriteMassBalanceReport (const Real new_time)
{

  if (not report_mass_balance) {
    return;
  }

  // Compute current mass in system
  ComputeMassAccum(1);

  const int offset = MFIXSpecies::NMAX;
  const int nspecies_g = fluid.nspecies();
  const int nspecies_s = solids.nspecies();

  if(ParallelDescriptor::IOProcessor()) {
    printf("\n**********");
    for (int col=0; col < 6; ++col) {
      printf("**************");
    }
    printf("****************\n");

    printf("  Species mass balance for interval:  %12.6f  to %12.6f\n",
           mass_balance_report_time, new_time);

    mass_balance_report_time = new_time;

    printf("\n  %-8s%14s%14s%14s%14s%14s%14s%14s\n", "Species",
           "mass(t+dt)", "mass(t) ", "production",
           "mass(in)", "mass(out)","net accu","flux   ");

    amrex::Array<amrex::Real,7> totals{0.};

    if (fluid.solve()) {
      for (int n=0; n < nspecies_g; ++n) {
        Real delta_accum( m_mass_accum[n+offset] - m_mass_accum[n] - m_mass_prod[n] );
        Real bc_flux( m_mass_inflow[n] - m_mass_outflow[n]);

        //net_acc  += delta_accum;
        //net_flux += bc_flux;

        printf("  %-8s%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e\n", fluid.species_names(n).c_str(),
               m_mass_accum[n+offset], m_mass_accum[n],m_mass_prod[n],
               m_mass_inflow[n], m_mass_outflow[n], delta_accum, bc_flux);

        totals[0] += m_mass_accum[n+offset];
        totals[1] += m_mass_accum[n];
        totals[2] += m_mass_prod[n];
        totals[3] += m_mass_inflow[n];
        totals[4] += m_mass_outflow[n];
        totals[5] += delta_accum;
        totals[6] += bc_flux;

        m_mass_accum[n] = m_mass_accum[n+offset];
        m_mass_inflow[n] = 0.;
        m_mass_outflow[n] = 0.;
        m_mass_prod[n] = 0.;
      }
    }

    if (m_dem.solve() || m_pic.solve()) {
      for (int n=0; n < nspecies_s; ++n) {
        Real delta_accum(pc->get_mass_accum(n+offset) - pc->get_mass_accum(n) - pc->get_mass_prod(n));
        Real bc_flux(pc->get_mass_inflow(n) - pc->get_mass_outflow(n));

        //net_acc += delta_accum;
        //net_flux += bc_flux;

        printf("  %-8s%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e\n", solids.species_names(n).c_str(),
               pc->get_mass_accum(n+offset), pc->get_mass_accum(n), pc->get_mass_prod(n),
               pc->get_mass_inflow(n), pc->get_mass_outflow(n),
               delta_accum, bc_flux);

        totals[0] += pc->get_mass_accum(n+offset);
        totals[1] += pc->get_mass_accum(n);
        totals[2] += pc->get_mass_prod(n);
        totals[3] += pc->get_mass_inflow(n);
        totals[4] += pc->get_mass_outflow(n);
        totals[5] += delta_accum;
        totals[6] += bc_flux;

        pc->ResetMassBalance(n);
      }
    }

    printf("----------");
    for (int col=0; col < 6; ++col) {
      printf("--------------");
    }
    printf("----------------\n");

    printf("  %-8s%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e\n", "Total",
           totals[0], totals[1], totals[2], totals[3], totals[4], totals[5], totals[6]);

    //printf("\n  total accu := %14.4e %77s\n", net_acc,"sum(mass(t+dt) - mass(t) - production)");
    //printf("    net flux := %14.4e %77s\n", net_flux,"sum(mass(in) - mass(out))");
    //printf("  difference := %14.4e %77s\n", totals[5]-totals[6],"(total acc) - (net flux))");

    if (m_dem.solve() || m_pic.solve()) {
      long tnp = pc->getTotalNumParticles();
      if (m_dem.solve() ) { printf("\n  total number of particles = "); }
      if (m_pic.solve() ) { printf("\n  total number of parcels =   "); }
      printf("%16s\n", MfixIO::FormatWithCommas(tnp).c_str() ) ;
    }

    Real diff = std::abs(totals[5]-totals[6]);
    printf("  |(total acc) - (net flux)| =  %14.4e\n", diff);

    printf("**********");
    for (int col=0; col < 6; ++col) {
      printf("**************");
    }
    printf("****************\n\n");
  }
}


void
MfixRW::ComputeMassAccum (const int offset)
{
  BL_PROFILE("mfix::ComputeMassAccum()");

  if (not report_mass_balance) {
    return;
  }


  if (fluid.solve()) {

    GpuArray<Real, MFIXSpecies::NMAX> accum = {0.};

    const int nspecies_g = fluid.nspecies();
    for (int lev = 0; lev < nlev; lev++) {

      const GpuArray<Real,3> dx = geom[lev].CellSizeArray();
      const Real vol = dx[0]*dx[1]*dx[2];

      for (MFIter mfi(*(m_leveldata[lev]->X_gk), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        EBCellFlagFab const& flagfab = ebfactory[lev]->getMultiEBCellFlagFab()[mfi];

        if (flagfab.getType(bx) != FabType::covered) {

          Array4<Real const> const& ep_g = m_leveldata[lev]->ep_g->const_array(mfi);
          Array4<Real const> const& ro_g = m_leveldata[lev]->ro_g->const_array(mfi);
          Array4<Real const> const& X_gk = m_leveldata[lev]->X_gk->const_array(mfi);

          if (flagfab.getType(bx) == FabType::singlevalued ) {

            Array4<Real const> const& volfrac = (ebfactory[lev]->getVolFrac()).const_array(mfi);

            for (int n=0; n < nspecies_g; ++n){

              ReduceOps<ReduceOpSum> reduce_op;
              ReduceData<Real> reduce_data(reduce_op);
              using ReduceTuple = typename decltype(reduce_data)::Type;

              reduce_op.eval(bx, reduce_data, [n, ep_g, ro_g, X_gk, volfrac]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                return { volfrac(i,j,k)*ep_g(i,j,k)*ro_g(i,j,k)*X_gk(i,j,k,n) };
              });

              ReduceTuple host_tuple = reduce_data.value(reduce_op);
              accum[n] += amrex::get<0>(host_tuple)*vol;

            } /* End loop over species */

          } else {

            for (int n=0; n < nspecies_g; ++n){

              ReduceOps<ReduceOpSum> reduce_op;
              ReduceData<Real> reduce_data(reduce_op);
              using ReduceTuple = typename decltype(reduce_data)::Type;

              reduce_op.eval(bx, reduce_data, [n, ep_g, ro_g, X_gk]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                return { ep_g(i,j,k)*ro_g(i,j,k)*X_gk(i,j,k,n) };
              });

              ReduceTuple host_tuple = reduce_data.value(reduce_op);
              accum[n] += amrex::get<0>(host_tuple)*vol;

            } /* End loop over species */
          } // FabType
        } // Not covered Fab
      } // loop over MFIter

    } // nlev

    // Global sum and copy to global variable with offset
    ParallelDescriptor::ReduceRealSum(accum.data(), nspecies_g);
    for (int n=0; n < nspecies_g; ++n) {
      m_mass_accum[n + offset*MFIXSpecies::NMAX] = accum[n];
    }

  }// solve fluid

  if (m_dem.solve() || m_pic.solve()) {
    pc->ComputeMassAccum(offset);
  }
}


void
MfixRW::ComputeMassProduction (const Real dt,
                               Vector< MultiFab const*> const& txfr)
{
  BL_PROFILE("mfix::ComputeMassProduction()");

  const int nspecies_g = fluid.nspecies();
  std::vector<Real> prod(nspecies_g, 0.);

  InterphaseTxfrIndexes txfr_idxs(nspecies_g, reactions.nreactions());

  const int scomp = txfr_idxs.chem_ro_gk;

  for (int lev = 0; lev < nlev; lev++) {

    const GpuArray<Real,3> dx = geom[lev].CellSizeArray();
    const Real vol = dx[0]*dx[1]*dx[2];

    for (int n=0; n < nspecies_g; ++n){

      ReduceOps<ReduceOpSum> reduce_op;
      ReduceData<Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      for (MFIter mfi(*txfr[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        EBCellFlagFab const& flagfab = ebfactory[lev]->getMultiEBCellFlagFab()[mfi];

        if (flagfab.getType(bx) != FabType::covered) {

          Array4<Real const> const& ep_g = m_leveldata[lev]->ep_g->const_array(mfi);
          Array4<Real const> const& ro_gk_txfr = txfr[lev]->const_array(mfi);

          if (flagfab.getType(bx) == FabType::singlevalued ) {

            Array4<Real const> const& volfrac = (ebfactory[lev]->getVolFrac()).const_array(mfi);

            reduce_op.eval(bx, reduce_data, [n, ep_g, ro_gk_txfr, scomp, volfrac]
            AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            { return { volfrac(i,j,k)*ep_g(i,j,k)*ro_gk_txfr(i,j,k,scomp+n) }; });

          } else {

            reduce_op.eval(bx, reduce_data, [n, ep_g, ro_gk_txfr, scomp]
            AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            { return { ep_g(i,j,k)*ro_gk_txfr(i,j,k,scomp+n) }; });

          } // check on FabType
        } // FabType is covered
      } // MFIter loop

      ReduceTuple host_tuple = reduce_data.value(reduce_op);
      prod[n] = amrex::get<0>(host_tuple)*vol;

    } /* End loop over species */

  } // nlev

  // Global sum and copy to global variable
  ParallelDescriptor::ReduceRealSum(prod.data(), nspecies_g);
  for (int n=0; n < nspecies_g; ++n) {
    m_mass_prod[n] += dt*prod[n];
  }
}


void
MfixRW::ComputeMassFlux (Vector< MultiFab const*> const& flux_x,
                         Vector< MultiFab const*> const& flux_y,
                         Vector< MultiFab const*> const& flux_z,
                         const int scomp,
                         const int ncomp,
                         const bool fluxes_are_area_weighted,
                         const int eb_has_flow,
                         Vector< MultiFab const*> const& eb_vel_in,
                         Vector< MultiFab const*> const& eb_species_in,
                         const Real dt)
{
  amrex::ignore_unused(ncomp, fluxes_are_area_weighted);

  const int nspecies_g = fluid.nspecies();

  std::vector<Real> mass_flow(2*nspecies_g, 0.);

  for (int lev = 0; lev < nlev; lev++) {

    Array<const MultiFab*,3> flux = {flux_x[lev], flux_y[lev], flux_z[lev]};

    const GpuArray<Real,3> dx = geom[lev].CellSizeArray();

    // This should be caught elsewhere but just in case...
    AMREX_ASSERT(dx[0] == dx[1] && dx[1] == dx[2]);
    Real const da( dx[0]*dx[0] );

    amrex::EBFArrayBoxFactory const& fact =
      static_cast<amrex::EBFArrayBoxFactory const&>(*ebfactory[lev]);

    Box domain(geom[lev].Domain());

    Array<Array4<int>,3> bct_lo = {bc_list.bc_ilo[lev]->array(),
                                   bc_list.bc_jlo[lev]->array(),
                                   bc_list.bc_klo[lev]->array()};

    Array<Array4<int>,3> bct_hi = {bc_list.bc_ihi[lev]->array(),
                                   bc_list.bc_jhi[lev]->array(),
                                   bc_list.bc_khi[lev]->array()};

    Array<Array<Array4<int>,3>,2> bct_extents = {bct_lo, bct_hi};

    auto areafrac = fact.getAreaFrac();

    for (int n=0; n < nspecies_g; ++n) {

      ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
      ReduceData<Real,Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      // Loop over MFIter
      for (MFIter mfi(*m_leveldata[lev]->ep_g, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];

        Array<IntVect,2> domain_extents = {domain.smallEnd(), domain.bigEnd()};


        // Loop over directions
        for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {

            // Get flux over direction dir
          Array4<Real const> const& flux_arr = flux[dir]->const_array(mfi);

          const Box& sbox = (*flux[dir])[mfi].box();

          // Save local box extents in a specific container
          Array<IntVect,2> sbox_extents = {sbox.smallEnd(), sbox.bigEnd()};

          // Loop over face normal
          for (int normal(-1); normal <= 1; normal += 2) {

            // index to map normal from {-1,1} to {0,1}
            const int idx = (normal + 1) / 2;

            // On domain face
            if (sbox_extents[idx][dir] == (domain_extents[idx][dir]+idx)) {

              // Create a copy of local box extents
              Array<IntVect,2> operative_box_extents(sbox_extents);

              // Modify operative box extents to make a 2D Box
              operative_box_extents[(idx+1)%2][dir] = operative_box_extents[idx][dir];

              // Create the 2D Box
              const Box operative_box(operative_box_extents[0],
                                      operative_box_extents[1]);

              // Get the boundary condition type on current face
              Array4<int>& bc_type = bct_extents[idx][dir];

              // Get domain index on this face
              const int domain_idx = domain_extents[idx][dir];

              if (flagfab.getType(bx) == FabType::regular ) {

                reduce_op.eval(operative_box, reduce_data, [dir,bc_type,
                domain_idx,normal,flux_arr,n,scomp]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                  Real m_in  = 0.;
                  Real m_out = 0.;

                  IntVect bct_cell(i,j,k);
                  bct_cell[dir] = domain_idx + normal;
                  const int bct = bc_type(bct_cell,0);

                  if(bct == BCList::pout) {
                    //m_out =  1.;
                    m_out =  normal*flux_arr(i,j,k,n+scomp);
                  }
                  else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                    //m_in  = 1.;
                    m_in  = -normal*flux_arr(i,j,k,n+scomp);
                  }

                  return {m_in, m_out};
                });


              } else if (flagfab.getType(bx) == FabType::singlevalued ) {

                Array4<Real const> const& areafrac_arr = areafrac[dir]->const_array(mfi);

                reduce_op.eval(operative_box, reduce_data, [dir,bc_type,
                domain_idx,normal,flux_arr,areafrac_arr,n,scomp]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                  Real m_in  = 0.;
                  Real m_out = 0.;

                  IntVect bct_cell(i,j,k);
                  bct_cell[dir] = domain_idx + normal;
                  const int bct = bc_type(bct_cell,0);

                  if(bct == BCList::pout) {
                    m_out =  normal*areafrac_arr(i,j,k)*flux_arr(i,j,k,n+scomp);
                    //m_out =  areafrac_arr(i,j,k);
                  }
                  else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                    m_in  = -normal*areafrac_arr(i,j,k)*flux_arr(i,j,k,n+scomp);
                    //m_in  = areafrac_arr(i,j,k);
                  }

                  return {m_in, m_out};
                });

              } // single valued fab
            } // on domain face
          } // loop over normal direction
        } // loop over direction

        // Add in mass flow from EB
        if( eb_has_flow && flagfab.getType(bx) == FabType::singlevalued) {

          // Area of eb face
          Array4<Real const> const& barea  = fact.getBndryArea().const_array(mfi);

          Array4<Real const> const& eb_vel     = eb_vel_in[lev]->const_array(mfi);
          Array4<Real const> const& eb_species = eb_species_in[lev]->const_array(mfi);

          reduce_op.eval(bx, reduce_data, [n, eb_vel, eb_species, barea]
            AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
          {

            Real vel_mag = std::sqrt(eb_vel(i,j,k,0)*eb_vel(i,j,k,0) +
                                     eb_vel(i,j,k,1)*eb_vel(i,j,k,1) +
                                     eb_vel(i,j,k,2)*eb_vel(i,j,k,2));

            return {vel_mag*eb_species(i,j,k,n)*barea(i,j,k), 0.};

          });

        } // eb_has_flow

      } // MFIter loop

      ReduceTuple host_tuple = reduce_data.value(reduce_op);

      // mass in
      mass_flow[n] += amrex::get<0>(host_tuple)*da;
      // mass out
      mass_flow[n+nspecies_g] += amrex::get<1>(host_tuple)*da;

    } // Loop over species

  } // loop over levels
#if 0
#endif


  ParallelDescriptor::ReduceRealSum(mass_flow.data(), 2*nspecies_g);

  // Copy into global variables.
  for (int n=0; n < nspecies_g; ++n) {
    m_mass_inflow[n]  += dt*mass_flow[n];
    m_mass_outflow[n] += dt*mass_flow[n+nspecies_g];
  }
}

void
MfixRW::InitMassBalance ()
{
  if (report_mass_balance) {
    for(int n(0); n<MFIXSpecies::NMAX; n++) {
        m_mass_accum[n] = m_mass_accum[n];
        m_mass_accum[n] = m_mass_accum[n+MFIXSpecies::NMAX];
        m_mass_inflow[n] = 0.;
        m_mass_outflow[n] = 0.;
        m_mass_prod[n] = 0.;
    }
  }
}

} // end namespace MfixIO
