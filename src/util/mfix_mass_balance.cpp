#include <mfix_rw.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_reactions_parms.H>

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

  const int offset = SPECIES::NMAX;
  const int nspecies_g = fluid.nspecies;


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
      printf("  %-8s%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e\n", fluid.species[n].c_str(),
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

  GpuArray<Real,SPECIES::NMAX> accum;

  const int nspecies_g = fluid.nspecies;
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
    mass_accum[n + offset*SPECIES::NMAX] = accum[n];
  }
}


void
MfixRW::ComputeMassProduction (const Real /*dt*/,
                               Vector< MultiFab const*> const& chem_txfr)
{
  BL_PROFILE("mfix::ComputeMassProduction()");

  const int nspecies_g = fluid.nspecies;
  std::vector<Real> prod(nspecies_g, 0.);

  ChemTransfer chem_txfr_idxs(nspecies_g, reactions.nreactions);

  const int scomp = chem_txfr_idxs.ro_gk_txfr;

  for (int lev = 0; lev < nlev; lev++) {

    const GpuArray<Real,3> dx = geom[lev].CellSizeArray();
    const Real vol = dx[0]*dx[1]*dx[2];

    const MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

     MultiFab const& ro_gk_txfr_fab = *(chem_txfr[lev]);

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

  const int nspecies_g = fluid.nspecies;

  std::vector<Real> mass_flow(2*nspecies_g, 0.);

  for (int lev = 0; lev < nlev; lev++) {

    const GpuArray<Real,3> dx = geom[lev].CellSizeArray();
    const Real dxdy = dx[0]*dx[1];
    const Real dxdz = dx[0]*dx[2];
    const Real dydz = dx[1]*dx[2];

    amrex::EBFArrayBoxFactory const& fact =
      static_cast<amrex::EBFArrayBoxFactory const&>(*ebfactory[lev]);

    Box domain(geom[lev].Domain());

    Array4<int> const& bct_ilo = bc_list.bc_ilo[lev]->array();
    Array4<int> const& bct_ihi = bc_list.bc_ihi[lev]->array();
    Array4<int> const& bct_jlo = bc_list.bc_jlo[lev]->array();
    Array4<int> const& bct_jhi = bc_list.bc_jhi[lev]->array();
    Array4<int> const& bct_klo = bc_list.bc_klo[lev]->array();
    Array4<int> const& bct_khi = bc_list.bc_khi[lev]->array();

    for (int n=0; n < nspecies_g; ++n) {

      ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
      ReduceData<Real, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      for (MFIter mfi(*(m_leveldata[lev]->ep_g), false); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Box& ubx = (*flux_x[lev])[mfi].box();
        IntVect ubx_lo(ubx.loVect());
        IntVect ubx_hi(ubx.hiVect());

        const Box& vbx = (*flux_y[lev])[mfi].box();
        IntVect vbx_lo(vbx.loVect());
        IntVect vbx_hi(vbx.hiVect());

        const Box& wbx = (*flux_z[lev])[mfi].box();
        IntVect wbx_lo(wbx.loVect());
        IntVect wbx_hi(wbx.hiVect());

        IntVect dom_lo(domain.loVect());
        IntVect dom_hi(domain.hiVect());

        Array4<Real const> const& flux_x_arr = flux_x[lev]->const_array(mfi);
        Array4<Real const> const& flux_y_arr = flux_y[lev]->const_array(mfi);
        Array4<Real const> const& flux_z_arr = flux_z[lev]->const_array(mfi);

        EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];

        if (flagfab.getType() == FabType::regular ) {

          if (dom_lo[0] == ubx_lo[0]) {

            // Create InVects for following Box
            IntVect ulo_bx_yz_lo(ubx_lo);
            IntVect ulo_bx_yz_hi(ubx_hi);

            // Fix lo and hi limits
            ulo_bx_yz_lo[0] = dom_lo[0];
            ulo_bx_yz_hi[0] = dom_lo[0];

            const Box ulo_bx_yz(ulo_bx_yz_lo, ulo_bx_yz_hi);

            reduce_op.eval(ulo_bx_yz, reduce_data,
              [bct_ilo, dom_lo, dydz,  flux_x_arr, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

                if(bct == BCList::pout) {
                  m_out -= flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                  m_in  += flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                return {m_in, m_out};
              });
          }

          if (ubx_hi[0] == (dom_hi[0] + 1)) {

            // Create InVects for following Box
            IntVect uhi_bx_yz_lo(ubx_lo);
            IntVect uhi_bx_yz_hi(ubx_hi);

            // Fix lo and hi limits
            uhi_bx_yz_lo[0] = dom_hi[0]+1;
            uhi_bx_yz_hi[0] = dom_hi[0]+1;

            const Box uhi_bx_yz(uhi_bx_yz_lo, uhi_bx_yz_hi);

            reduce_op.eval(uhi_bx_yz, reduce_data,
              [bct_ihi, dom_hi, dydz,  flux_x_arr, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

                if(bct == BCList::pout) {
                  m_out += flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                  m_in  -= flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                return {m_in, m_out};
              });

          }


          if (dom_lo[1] == vbx_lo[1]) {

              // Create InVects for following Box
              IntVect vlo_bx_xz_lo(vbx_lo);
              IntVect vlo_bx_xz_hi(vbx_hi);

              // Fix lo and hi limits
              vlo_bx_xz_lo[1] = dom_lo[1];
              vlo_bx_xz_hi[1] = dom_lo[1];

              const Box vlo_bx_xz(vlo_bx_xz_lo, vlo_bx_xz_hi);
              reduce_op.eval(vlo_bx_xz, reduce_data,
                [bct_jlo, dom_lo, dxdz,  flux_y_arr, n, scomp]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {

                  Real m_in  = 0.;
                  Real m_out = 0.;

                  const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

                  if(bct == BCList::pout) {
                    m_out -= flux_y_arr(i,j,k,n+scomp)*dxdz;
                  }
                  else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                    m_in  += flux_y_arr(i,j,k,n+scomp)*dxdz;
                  }
                  return {m_in, m_out};
                });

          }

          if (vbx_hi[1] == (dom_hi[1] + 1)) {

            // Create InVects for following Box
            IntVect vhi_bx_xz_lo(vbx_lo);
            IntVect vhi_bx_xz_hi(vbx_hi);

            // Fix lo and hi limits
            vhi_bx_xz_lo[1] = dom_hi[1]+1;
            vhi_bx_xz_hi[1] = dom_hi[1]+1;

            const Box vhi_bx_xz(vhi_bx_xz_lo, vhi_bx_xz_hi);

            reduce_op.eval(vhi_bx_xz, reduce_data,
              [bct_jhi, dom_hi, dxdz,  flux_y_arr, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {

                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

                if(bct == BCList::pout) {
                  m_out += flux_y_arr(i,j,k,n+scomp)*dxdz;
                }
                else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                  m_in  -= flux_y_arr(i,j,k,n+scomp)*dxdz;
                }
                return {m_in, m_out};
              });

          }

          if (dom_lo[2] == wbx_lo[2]) {

            // Create InVects for following Boxes
            IntVect wlo_bx_xy_lo(wbx_lo);
            IntVect wlo_bx_xy_hi(wbx_hi);

            // Fix lo and hi limits
            wlo_bx_xy_lo[2] = dom_lo[2];
            wlo_bx_xy_hi[2] = dom_lo[2];

            const Box wlo_bx_xy(wlo_bx_xy_lo, wlo_bx_xy_hi);

            reduce_op.eval(wlo_bx_xy, reduce_data,
              [bct_klo, dom_lo, dxdy,  flux_z_arr, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_klo(i,j,dom_lo[2]-1,0);

                if(bct == BCList::pout) {
                  m_out -= flux_z_arr(i,j,k,n+scomp)*dxdy;
                }
                else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                  m_in  += flux_z_arr(i,j,k,n+scomp)*dxdy;
                }
                return {m_in, m_out};
              });
          }

          if (wbx_hi[2] == (dom_hi[2] + 1)) {

            // Create InVects for following Boxes
            IntVect whi_bx_xy_lo(wbx_lo);
            IntVect whi_bx_xy_hi(wbx_hi);

            // Fix lo and hi limits
            whi_bx_xy_lo[2] = dom_hi[2]+1;
            whi_bx_xy_hi[2] = dom_hi[2]+1;

            const Box whi_bx_xy(whi_bx_xy_lo, whi_bx_xy_hi);

            reduce_op.eval(whi_bx_xy, reduce_data,
               [bct_khi, dom_hi, dxdy,  flux_z_arr, n, scomp]
               AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
               {
                 Real m_in  = 0.;
                 Real m_out = 0.;

                 const int bct = bct_khi(i,j,dom_hi[2]+1,0);

                 if(bct == BCList::pout) {
                   m_out += flux_z_arr(i,j,k,n+scomp)*dxdy;
                 }
                 else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                   m_in  -= flux_z_arr(i,j,k,n+scomp)*dxdy;
                 }
                 return {m_in, m_out};
               });
          }

        } else if (flagfab.getType(amrex::grow(bx,1)) != FabType::covered ) {

          Array4<Real const> apx = fact.getAreaFrac()[0]->const_array(mfi);
          Array4<Real const> apy = fact.getAreaFrac()[1]->const_array(mfi);
          Array4<Real const> apz = fact.getAreaFrac()[2]->const_array(mfi);

          if (dom_lo[0] == ubx_lo[0]) {

            // Create InVects for following Box
            IntVect ulo_bx_yz_lo(ubx_lo);
            IntVect ulo_bx_yz_hi(ubx_hi);

            // Fix lo and hi limits
            ulo_bx_yz_lo[0] = dom_lo[0];
            ulo_bx_yz_hi[0] = dom_lo[0];

            const Box ulo_bx_yz(ulo_bx_yz_lo, ulo_bx_yz_hi);

            reduce_op.eval(ulo_bx_yz, reduce_data,
              [bct_ilo, dom_lo, dydz, apx,  flux_x_arr, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {

                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

                if(bct == BCList::pout) {
                  m_out -= apx(i,j,k)*flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                  m_in  += apx(i,j,k)*flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                return {m_in, m_out};
              });

          }

          if (ubx_hi[0] == (dom_hi[0] + 1)) {

            // Create InVects for following Box
            IntVect uhi_bx_yz_lo(ubx_lo);
            IntVect uhi_bx_yz_hi(ubx_hi);

            // Fix lo and hi limits
            uhi_bx_yz_lo[0] = dom_hi[0]+1;
            uhi_bx_yz_hi[0] = dom_hi[0]+1;

            const Box uhi_bx_yz(uhi_bx_yz_lo, uhi_bx_yz_hi);

            reduce_op.eval(uhi_bx_yz, reduce_data,
              [bct_ihi, dom_hi, dydz,  flux_x_arr, apx, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

                if(bct == BCList::pout) {
                  m_out += apx(i,j,k)*flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                  m_in  -= apx(i,j,k)*flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                return {m_in, m_out};
              });

          }


          if (dom_lo[1] == vbx_lo[1]) {

              // Create InVects for following Box
              IntVect vlo_bx_xz_lo(vbx_lo);
              IntVect vlo_bx_xz_hi(vbx_hi);

              // Fix lo and hi limits
              vlo_bx_xz_lo[1] = dom_lo[1];
              vlo_bx_xz_hi[1] = dom_lo[1];

              const Box vlo_bx_xz(vlo_bx_xz_lo, vlo_bx_xz_hi);
              reduce_op.eval(vlo_bx_xz, reduce_data,
                [bct_jlo, dom_lo, dxdz,  flux_y_arr, apy, n, scomp]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {

                  Real m_in  = 0.;
                  Real m_out = 0.;

                  const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

                  if(bct == BCList::pout) {
                    m_out -= apy(i,j,k)*flux_y_arr(i,j,k,n+scomp)*dxdz;
                  }
                  else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                    m_in  += apy(i,j,k)*flux_y_arr(i,j,k,n+scomp)*dxdz;
                  }
                  return {m_in, m_out};
                });

          }

          if (vbx_hi[1] == (dom_hi[1] + 1)) {

            // Create InVects for following Box
            IntVect vhi_bx_xz_lo(vbx_lo);
            IntVect vhi_bx_xz_hi(vbx_hi);

            // Fix lo and hi limits
            vhi_bx_xz_lo[1] = dom_hi[1]+1;
            vhi_bx_xz_hi[1] = dom_hi[1]+1;

            const Box vhi_bx_xz(vhi_bx_xz_lo, vhi_bx_xz_hi);

            reduce_op.eval(vhi_bx_xz, reduce_data,
              [bct_jhi, dom_hi, dxdz,  flux_y_arr, apy, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {

                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

                if(bct == BCList::pout) {
                  m_out += apy(i,j,k)*flux_y_arr(i,j,k,n+scomp)*dxdz;
                }
                else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                  m_in  -= apy(i,j,k)*flux_y_arr(i,j,k,n+scomp)*dxdz;
                }
                return {m_in, m_out};
              });

          }

          if (dom_lo[2] == wbx_lo[2]) {

            // Create InVects for following Boxes
            IntVect wlo_bx_xy_lo(wbx_lo);
            IntVect wlo_bx_xy_hi(wbx_hi);

            // Fix lo and hi limits
            wlo_bx_xy_lo[2] = dom_lo[2];
            wlo_bx_xy_hi[2] = dom_lo[2];

            const Box wlo_bx_xy(wlo_bx_xy_lo, wlo_bx_xy_hi);

            reduce_op.eval(wlo_bx_xy, reduce_data,
              [bct_klo, dom_lo, dxdy,  flux_z_arr, apz, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_klo(i,j,dom_lo[2]-1,0);

                if(bct == BCList::pout) {
                  m_out -= apz(i,j,k)*flux_z_arr(i,j,k,n+scomp)*dxdy;
                }
                else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                  m_in  += apz(i,j,k)*flux_z_arr(i,j,k,n+scomp)*dxdy;
                }
                return {m_in, m_out};
              });
          }

          if (wbx_hi[2] == (dom_hi[2] + 1)) {

            // Create InVects for following Boxes
            IntVect whi_bx_xy_lo(wbx_lo);
            IntVect whi_bx_xy_hi(wbx_hi);

            // Fix lo and hi limits
            whi_bx_xy_lo[2] = dom_hi[2]+1;
            whi_bx_xy_hi[2] = dom_hi[2]+1;

            const Box whi_bx_xy(whi_bx_xy_lo, whi_bx_xy_hi);

            reduce_op.eval(whi_bx_xy, reduce_data,
               [bct_khi, dom_hi, dxdy,  flux_z_arr, apz, n, scomp]
               AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
               {
                 Real m_in  = 0.;
                 Real m_out = 0.;

                 const int bct = bct_khi(i,j,dom_hi[2]+1,0);

                 if(bct == BCList::pout) {
                   m_out += apz(i,j,k)*flux_z_arr(i,j,k,n+scomp)*dxdy;
                 }
                 else if ((bct == BCList::minf) || (bct == BCList::pinf)) {
                   m_in  -= apz(i,j,k)*flux_z_arr(i,j,k,n+scomp)*dxdy;
                 }
                 return {m_in, m_out};
               });
          }
        } // Not regular Fab
      } // loop over MFIter

      ReduceTuple host_tuple = reduce_data.value(reduce_op);

      mass_flow[n] += amrex::get<0>(host_tuple);
      mass_flow[n+nspecies_g] += amrex::get<1>(host_tuple);

    } // Loop over species
  } // Loop over levels

  ParallelDescriptor::ReduceRealSum(mass_flow.data(), 2*nspecies_g);

  // Copy into global variables.
  for (int n=0; n < nspecies_g; ++n) {
    mass_inflow[n]  += dt*mass_flow[n];
    mass_outflow[n] += dt*mass_flow[n+nspecies_g];
  }
}

} // end namespace MfixIO
