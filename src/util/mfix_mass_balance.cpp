#include <mfix_rw.H>
#include <mfix_fluid.H>
#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_monitors.H>

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


  if(ParallelDescriptor::IOProcessor()) {
    printf("\n**********");
    for (int col=0; col < 7; ++col) {
      printf("**************");
    }
    printf("****************\n");

    printf("  Species mass balance for interval:  %12.6f  to %12.6f\n",
           mass_balance_report_time, new_time);

    mass_balance_report_time = new_time;

    Real tot_flux_in = 0.;
    std::vector<Real> error(nspecies_g, 100.);
    std::vector<Real> delta_accum(nspecies_g, 0.);
    std::vector<Real> bc_flux(nspecies_g, 0.);
    for (int n=0; n < nspecies_g; ++n) {
      delta_accum[n] = mass_accum[n+offset] - mass_accum[n] - mass_prod[n];
      bc_flux[n] = mass_inflow[n] - mass_outflow[n];
      tot_flux_in += mass_inflow[n];
    }

    printf("\n  %-8s%14s%14s%14s%14s%14s%14s%14s%14s\n", "Species",
           "mass(t+dt)", "mass(t) ", "production",
           "mass(in)", "mass(out)","net accu","flux   ","% error  ");

    for (int n=0; n < nspecies_g; ++n) {
      if( tot_flux_in > 0.) {
        error[n] = Math::abs((bc_flux[n] - delta_accum[n]) / tot_flux_in )*100.0;
      }
      printf("  %-8s%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e\n", fluid.species_names(n).c_str(),
             mass_accum[n+offset], mass_accum[n],mass_prod[n],
             mass_inflow[n], mass_outflow[n],
             delta_accum[n], bc_flux[n], error[n]);


      mass_accum[n] = mass_accum[n+offset];
      mass_inflow[n] = 0.;
      mass_outflow[n] = 0.;
      mass_prod[n] = 0.;
    }

    printf("\n  net accu := mass(t+dt) - mass(t) - production\n");
    printf("      flux := mass(in) - mass(out)\n");
    printf("   %% error := 100*(flux - net accu) / (total flux)\n");

    printf("**********");
    for (int col=0; col < 7; ++col) {
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

  GpuArray<Real, MFIXSpecies::NMAX> accum;

  const int nspecies_g = fluid.nspecies();
  for (int lev = 0; lev < nlev; lev++) {

    const GpuArray<Real,3> dx = geom[lev].CellSizeArray();
    const Real vol = dx[0]*dx[1]*dx[2];

    const MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

    MultiFab      & ep_g = *(m_leveldata[lev]->ep_g);
    MultiFab const& ro_g = *(m_leveldata[lev]->ro_g);

    // Convert "ep_g" into (rho * ep_g) -- no ghosts
    MultiFab::Multiply(ep_g, ro_g, 0, 0, 1, 0);

    for (int n=0; n < nspecies_g; ++n){

      accum[n] = 0.;

      MultiFab const& X_gk = *(m_leveldata[lev]->X_gk);

      accum[n] = amrex::ReduceSum(*volfrac, ep_g, X_gk, 0,
      [=] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                 Array4<const Real> const & vfrc,
                                 Array4<const Real> const & ep_ro,
                                 Array4<const Real> const & Xgk) -> amrex::Real
      {
        Real dm = 0.;
        amrex::Loop(bx, [n,vfrc,ep_ro,Xgk,&dm] (int i, int j, int k) noexcept
          {if(vfrc(i,j,k) > 0.0) dm += vfrc(i,j,k)*ep_ro(i,j,k)*Xgk(i,j,k,n);});
        return dm;
      });
      accum[n] *= vol;

    } /* End loop over species */

      // Convert (rho * ep_g) back into ep_g
    MultiFab::Divide(ep_g, ro_g, 0, 0, 1, 0);

  } // nlev


  // Global sum and copy to global variable with offset
  ParallelDescriptor::ReduceRealSum(accum.data(), nspecies_g);
  for (int n=0; n < nspecies_g; ++n) {
    mass_accum[n + offset*MFIXSpecies::NMAX] = accum[n];
  }
}


void
MfixRW::ComputeMassProduction (const Real /*dt*/,
                               Vector< MultiFab const*> const& txfr)
{
  BL_PROFILE("mfix::ComputeMassProduction()");

  const int nspecies_g = fluid.nspecies();
  std::vector<Real> prod(nspecies_g, 0.);

  Transfer txfr_idxs(nspecies_g, reactions.nreactions());

  const int scomp = txfr_idxs.ro_gk_txfr;

  for (int lev = 0; lev < nlev; lev++) {

    const GpuArray<Real,3> dx = geom[lev].CellSizeArray();
    const Real vol = dx[0]*dx[1]*dx[2];

    const MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

     MultiFab const& ro_gk_txfr_fab = *(txfr[lev]);

    for (int n=0; n < nspecies_g; ++n){

      prod[n] = amrex::ReduceSum(*volfrac, ro_gk_txfr_fab, 0,
      [=] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                 Array4<const Real> const & vfrc,
                                 Array4<const Real> const & ro_gk_txfr) -> amrex::Real
      {
        Real dm = 0.;
        amrex::Loop(bx, [n,vfrc,ro_gk_txfr,scomp,&dm] (int i, int j, int k) noexcept
        {if(vfrc(i,j,k) > 0.0) dm += vfrc(i,j,k)*ro_gk_txfr(i,j,k,n+scomp);});
        return dm;
      });
      prod[n] *= vol;

    } /* End loop over species */
  } // nlev

  // Global sum and copy to global variable
  ParallelDescriptor::ReduceRealSum(prod.data(), nspecies_g);
  for (int n=0; n < nspecies_g; ++n) {
//    mass_prod[n] += dt*prod[n];
  }
}


void
MfixRW::ComputeMassFlux (Vector< MultiFab const*> const& flux_x,
                         Vector< MultiFab const*> const& flux_y,
                         Vector< MultiFab const*> const& flux_z,
                         const int scomp,
                         const int ncomp,
                         const bool fluxes_are_area_weighted,
                         const Real dt)
{
  amrex::ignore_unused(ncomp, fluxes_are_area_weighted);

  const int nspecies_g = fluid.nspecies();

  std::vector<Real> mass_flow(2*nspecies_g, 0.);

  for (int lev = 0; lev < nlev; lev++) {

    Array<const MultiFab*,3> flux = {flux_x[lev], flux_y[lev], flux_z[lev]};

    const GpuArray<Real,3> dx = geom[lev].CellSizeArray();
    GpuArray<Real,3> da = {dx[1]*dx[2], dx[0]*dx[2], dx[0]*dx[1]};

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

    // Loop over MFIter
    for (MFIter mfi(*flux_x[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

      EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];

      // Only if fab is not covered
      if (flagfab.getType(amrex::grow(bx,1)) != FabType::covered) {

        Array<IntVect,2> domain_extents = {domain.smallEnd(), domain.bigEnd()};

        auto areafrac = fact.getAreaFrac();

        // Loop over directions
        for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {

          // Get flux over direction dir
          Array4<Real const> const& flux_arr = flux[dir]->const_array(mfi);
          Array4<Real const> const& areafrac_arr = areafrac[dir]->const_array(mfi);

          const Box& local_box = (*flux[dir])[mfi].box();

          // Save local box extents in a specific container
          Array<IntVect,2> local_box_extents = {local_box.smallEnd(), local_box.bigEnd()};

          // Save a copy of facet area along this direction
          const Real dA = da[dir];

          // Loop over face normal
          for (int normal(-1); normal <= 1; normal += 2) {

            // index to map normal from {-1,1} to {0,1}
            const int idx = (normal + 1) / 2;

            // On low face
            if (local_box_extents[idx][dir] == (domain_extents[idx][dir]+idx)) {

              for (int n=0; n < nspecies_g; ++n) {

                ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
                ReduceData<Real,Real> reduce_data(reduce_op);
                using ReduceTuple = typename decltype(reduce_data)::Type;

                // Create a copy of local box extents
                Array<IntVect,2> operative_box_extents(local_box_extents);

                // Modify operative box extents to make a 2D Box
                operative_box_extents[(idx+1)%2][dir] = operative_box_extents[idx][dir];

                // Create the 2D Box
                const Box operative_box(operative_box_extents[0],
                                        operative_box_extents[1]);

                // Get the boundary condition type on current face
                Array4<int>& bc_type = bct_extents[idx][dir];

                // Get domain index on this face
                const int domain_idx = domain_extents[idx][dir];

                reduce_op.eval(operative_box, reduce_data, [dir,bc_type,
                    domain_idx,dA,normal,flux_arr,areafrac_arr,n,scomp]
                  AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                  Real m_in  = 0.;
                  Real m_out = 0.;

                  IntVect bct_cell(i,j,k);
                  bct_cell[dir] = domain_idx + normal;
                  const int bct = bc_type(bct_cell,0);

                  if(bct == BCList::pout) {
                    m_out += normal*areafrac_arr(i,j,k)*flux_arr(i,j,k,n+scomp)*dA;
                  }
                  else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                    m_in  += normal*areafrac_arr(i,j,k)*flux_arr(i,j,k,n+scomp)*dA;
                  }

                  return {m_in, m_out};
                });

                ReduceTuple host_tuple = reduce_data.value(reduce_op);

                Real mass_in = amrex::get<0>(host_tuple);
                Real mass_out = amrex::get<1>(host_tuple);

                amrex::ParallelDescriptor::ReduceRealSum(&mass_in, 1);
                amrex::ParallelDescriptor::ReduceRealSum(&mass_out, 1);

                mass_flow[n] += mass_in;
                mass_flow[n+nspecies_g] += mass_out;
              }

            } // Loop over species
          } // loop over normal direction
        } // loop over dir
      } // Not covered Fab
    } // loop over MFIter
  } // Loop over levels

  ParallelDescriptor::ReduceRealSum(mass_flow.data(), 2*nspecies_g);

  // Copy into global variables.
  for (int n=0; n < nspecies_g; ++n) {
    mass_inflow[n]  += dt*mass_flow[n];
    mass_outflow[n] += dt*mass_flow[n+nspecies_g];
  }
}

} // end namespace MfixIO
