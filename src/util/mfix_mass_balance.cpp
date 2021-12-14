#include <mfix.H>
#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>


void
mfix::WriteMassBalanceReport(const Real dt)
{

  if (not report_mass_balance) {
    return;
  }

  const int offset = SPECIES::NMAX;
  const int nspecies_g = fluid.nspecies;

  std::vector<Real> error(nspecies_g, 0.);


  amrex::Print () << "Fluid name is " << fluid.species[0] << std::endl;

  if(ParallelDescriptor::IOProcessor()) {
    printf("\n****************");
    for (int n=0; n < nspecies_g; ++n) {
      printf("**************");
    }

    printf("\n  Species mass balance:\n");

    Real total_dm = 0.;

    printf("\n%-14.14s", " ");
    // printf("\n              ");
    for (int n=0; n < nspecies_g; ++n) {
      printf("%12s  ", fluid.species[n].c_str());
      error[n] = mass_accum[n+offset];
    }

    printf("\n%15s:", "mass (t)");
    for (int n=0; n < nspecies_g; ++n) {
      printf("  %12.4e", mass_accum[n]);
      error[n] -= mass_accum[n];
    }

    printf("\n%15s:", "inflow");
    for (int n=0; n < nspecies_g; ++n) {
      printf("  %12.4e", dt*mass_inflow[n]);
      error[n] += dt*mass_inflow[n];
    }

    printf("\n%15s:", "outflow");
    for (int n=0; n < nspecies_g; ++n) {
      printf("  %12.4e", dt*mass_outflow[n]);
      error[n] -= dt*mass_outflow[n];
    }

    if (reactions.solve) {
    }

    printf("\n%15s:", "mass (t+dt)");
    for (int n=0; n < nspecies_g; ++n) {
      printf("  %12.4e", mass_accum[n+offset]);
    }

    printf("\n%15s:", "error");
    for (int n=0; n < nspecies_g; ++n) {
      printf("  %12.4e", error[n]);
      total_dm += error[n];
    }



    // for (int n=0; n < nspecies_g; ++n) {

    //   // Current mass - previous mass
    //   Real dm = (mass_accum[n+offset] - mass_accum[n]);
    //   dm += dt*mass_inflow[n];  // Add in flow from boundaries
    //   dm -= dt*mass_outflow[n]; // Mass loss from outflows

    //   //printf("      %2d   %12.4e\n", n, dm);
    //   printf("      %2d   %12.4e  =%12.4e  +%12.4e  -%12.4e\n", n,
    //          mass_accum[n+offset], mass_accum[n],
    //          mass_inflow[n], mass_outflow[n]);

    //   mass_accum[n] = mass_accum[n+offset];
    //   total_dm += dm;
    // }

    printf("\n\n   Total   %12.4e\n", total_dm);
    for (int n=0; n < nspecies_g; ++n) {
      printf("**************");
    }
    printf("****************\n");
  }

}



void
mfix::ComputeMassAccum ( const int offset )
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
    // MultiFab::Multiply(ep_g, ro_g, 0, 0, 1, 0);

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
        amrex::Loop(bx, [n,vol,vfrc,ep_ro,Xgk,&dm] (int i, int j, int k) noexcept
          {if(vfrc(i,j,k) > 0.0) dm += vfrc(i,j,k)*ep_ro(i,j,k)*Xgk(i,j,k,n);});
        return dm;
      });
      accum[n] *= vol;

    } /* End loop over species */

      // Convert (rho * ep_g) back into ep_g
    // MultiFab::Divide(ep_g, ro_g, 0, 0, 1, 0);

  } // nlev


  // Global sum and copy to global variable with offset
  ParallelDescriptor::ReduceRealSum(accum.data(), nspecies_g);
  for (int n=0; n < nspecies_g; ++n) {
    mass_accum[n + offset*SPECIES::NMAX] = accum[n];
  }

  if( !offset ) {
    amrex::Print() << "Computing initial mass for each species:\n";
    if(ParallelDescriptor::IOProcessor()) {
      printf("****************************************\n");
      printf("  Species mass report:\n");
      for (int n=0; n < nspecies_g; ++n){
        printf("   %2d   %12.6e\n", n, mass_accum[n]);
      }
      printf("****************************************\n\n");
    }
  }

}





void
mfix::ComputeMassFlux (Vector< MultiFab const*> const& flux_x,
                       Vector< MultiFab const*> const& flux_y,
                       Vector< MultiFab const*> const& flux_z,
                       const int scomp,
                       const int ncomp,
                       const bool fluxes_are_area_weighted,
                       const int offset)
{


  GpuArray<Real,SPECIES::NMAX> mass_in, mass_out;

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  const int nspecies_g = fluid.nspecies;

  for (int lev = 0; lev < nlev; lev++) {

    const GpuArray<Real,3> dx = geom[lev].CellSizeArray();
    const Real dxdy = dx[0]*dx[1];
    const Real dxdz = dx[0]*dx[2];
    const Real dydz = dx[1]*dx[2];

    auto const& fact = EBFactory(lev);

    Box domain(geom[lev].Domain());

    Array4<int> const& bct_ilo = bc_ilo[lev]->array();
    Array4<int> const& bct_ihi = bc_ihi[lev]->array();
    Array4<int> const& bct_jlo = bc_jlo[lev]->array();
    Array4<int> const& bct_jhi = bc_jhi[lev]->array();
    Array4<int> const& bct_klo = bc_klo[lev]->array();
    Array4<int> const& bct_khi = bc_khi[lev]->array();

    for (int n=0; n < nspecies_g; ++n) {

      mass_in[n]  = 0.;
      mass_out[n] = 0.;

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

        const int nlft = amrex::max(0, dom_lo[0]-ubx_lo[0]);
        const int nbot = amrex::max(0, dom_lo[1]-vbx_lo[1]);
        const int ndwn = amrex::max(0, dom_lo[2]-wbx_lo[2]);

        const int nrgt = amrex::max(0, ubx_hi[0]-dom_hi[0]);
        const int ntop = amrex::max(0, vbx_hi[1]-dom_hi[1]);
        const int nup  = amrex::max(0, wbx_hi[2]-dom_hi[2]);

        Array4<Real const> const& flux_x_arr = flux_x[lev]->const_array(mfi);
        Array4<Real const> const& flux_y_arr = flux_y[lev]->const_array(mfi);
        Array4<Real const> const& flux_z_arr = flux_z[lev]->const_array(mfi);

        EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
        Array4<EBCellFlag const> const& flag = flagfab.const_array();

        if (flagfab.getType() == FabType::regular ) {

          if (nlft == 0) { // Fluxes don't have ghost cells

            // Create InVects for following Box
            IntVect ulo_bx_yz_lo(ubx_lo);
            IntVect ulo_bx_yz_hi(ubx_hi);

            // Fix lo and hi limits
            ulo_bx_yz_lo[0] = dom_lo[0];
            ulo_bx_yz_hi[0] = dom_lo[0];

            const Box ulo_bx_yz(ulo_bx_yz_lo, ulo_bx_yz_hi);

            reduce_op.eval(ulo_bx_yz, reduce_data,
              [bct_ilo, dom_lo, dydz, minf, pinf, pout, flux_x_arr, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

                if(bct == pout) {
                  m_out -= flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                else if ((bct == minf) || (bct == pinf)) {
                  m_in  += flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                return {m_in, m_out};
              });
          } // nlft

          if (nrgt == 1) {

            // Create InVects for following Box
            IntVect uhi_bx_yz_lo(ubx_lo);
            IntVect uhi_bx_yz_hi(ubx_hi);

            // Fix lo and hi limits
            uhi_bx_yz_lo[0] = dom_hi[0]+1;
            uhi_bx_yz_hi[0] = dom_hi[0]+1;

            const Box uhi_bx_yz(uhi_bx_yz_lo, uhi_bx_yz_hi);

            reduce_op.eval(uhi_bx_yz, reduce_data,
              [bct_ihi, dom_hi, dydz, minf, pinf, pout, flux_x_arr, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

                if(bct == pout) {
                  m_out += flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                else if ((bct == minf) || (bct == pinf)) {
                  m_in  -= flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                return {m_in, m_out};
              });

          } // nrgt


          if (nbot == 0) { // dom_lo[1] == vbx_lo[1]

              // Create InVects for following Box
              IntVect vlo_bx_xz_lo(vbx_lo);
              IntVect vlo_bx_xz_hi(vbx_hi);

              // Fix lo and hi limits
              vlo_bx_xz_lo[1] = dom_lo[1];
              vlo_bx_xz_hi[1] = dom_lo[1];

              const Box vlo_bx_xz(vlo_bx_xz_lo, vlo_bx_xz_hi);
              reduce_op.eval(vlo_bx_xz, reduce_data,
                [bct_jlo, dom_lo, dxdz, minf, pinf, pout, flux_y_arr, n, scomp]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {

                  Real m_in  = 0.;
                  Real m_out = 0.;

                  const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

                  if(bct == pout) {
                    m_out -= flux_y_arr(i,j,k,n+scomp)*dxdz;
                  }
                  else if ((bct == minf) || (bct == pinf)) {
                    m_in  += flux_y_arr(i,j,k,n+scomp)*dxdz;
                  }
                  return {m_in, m_out};
                });

          } // nbot

          if (ntop == 1) {

            // Create InVects for following Box
            IntVect vhi_bx_xz_lo(vbx_lo);
            IntVect vhi_bx_xz_hi(vbx_hi);

            // Fix lo and hi limits
            vhi_bx_xz_lo[1] = dom_hi[1]+1;
            vhi_bx_xz_hi[1] = dom_hi[1]+1;

            const Box vhi_bx_xz(vhi_bx_xz_lo, vhi_bx_xz_hi);

            reduce_op.eval(vhi_bx_xz, reduce_data,
              [bct_jhi, dom_hi, dxdz, minf, pinf, pout, flux_y_arr, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {

                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

                if(bct == pout) {
                  m_out += flux_y_arr(i,j,k,n+scomp)*dxdz;
                }
                else if ((bct == minf) || (bct == pinf)) {
                  m_in  -= flux_y_arr(i,j,k,n+scomp)*dxdz;
                }
                return {m_in, m_out};
              });

          } // ntop

          if (ndwn == 0) {

            // Create InVects for following Boxes
            IntVect wlo_bx_xy_lo(wbx_lo);
            IntVect wlo_bx_xy_hi(wbx_hi);

            // Fix lo and hi limits
            wlo_bx_xy_lo[2] = dom_lo[2];
            wlo_bx_xy_hi[2] = dom_lo[2];

            const Box wlo_bx_xy(wlo_bx_xy_lo, wlo_bx_xy_hi);

            reduce_op.eval(wlo_bx_xy, reduce_data,
              [bct_klo, dom_lo, dxdy, minf, pinf, pout, flux_z_arr, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_klo(i,j,dom_lo[2]-1,0);

                if(bct == pout) {
                  m_out -= flux_z_arr(i,j,k,n+scomp)*dxdy;
                }
                else if ((bct == minf) || (bct == pinf)) {
                  m_in  += flux_z_arr(i,j,k,n+scomp)*dxdy;
                }
                return {m_in, m_out};
              });
          } // ndwn

          if (nup == 0) {

            // Create InVects for following Boxes
            IntVect whi_bx_xy_lo(wbx_lo);
            IntVect whi_bx_xy_hi(wbx_hi);

            // Fix lo and hi limits
            whi_bx_xy_lo[2] = dom_hi[2]+1;
            whi_bx_xy_hi[2] = dom_hi[2]+1;

            const Box whi_bx_xy(whi_bx_xy_lo, whi_bx_xy_hi);

            reduce_op.eval(whi_bx_xy, reduce_data,
               [bct_khi, dom_hi, dxdy, minf, pinf, pout, flux_z_arr, n, scomp]
               AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
               {
                 Real m_in  = 0.;
                 Real m_out = 0.;

                 const int bct = bct_khi(i,j,dom_hi[2]+1,0);

                 if(bct == pout) {
                   m_out += flux_z_arr(i,j,k,n+scomp)*dxdy;
                 }
                 else if ((bct == minf) || (bct == pinf)) {
                   m_in  -= flux_z_arr(i,j,k,n+scomp)*dxdy;
                 }
                 return {m_in, m_out};
               });
          } // nup

        } else if (flagfab.getType(amrex::grow(bx,1)) != FabType::covered ) {

          Array4<Real const> apx = fact.getAreaFrac()[0]->const_array(mfi);
          Array4<Real const> apy = fact.getAreaFrac()[1]->const_array(mfi);
          Array4<Real const> apz = fact.getAreaFrac()[2]->const_array(mfi);

          if (nlft == 0) {

            // Create InVects for following Box
            IntVect ulo_bx_yz_lo(ubx_lo);
            IntVect ulo_bx_yz_hi(ubx_hi);

            // Fix lo and hi limits
            ulo_bx_yz_lo[0] = dom_lo[0];
            ulo_bx_yz_hi[0] = dom_lo[0];

            const Box ulo_bx_yz(ulo_bx_yz_lo, ulo_bx_yz_hi);

            reduce_op.eval(ulo_bx_yz, reduce_data,
              [bct_ilo, dom_lo, dydz, apx, minf, pinf, pout, flux_x_arr, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {

                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

                if(bct == pout) {
                  m_out -= apx(i,j,k)*flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                else if ((bct == minf) || (bct == pinf)) {
                  m_in  += apx(i,j,k)*flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                return {m_in, m_out};
              });

          }//nlft

          if (nrgt == 1) {

            // Create InVects for following Box
            IntVect uhi_bx_yz_lo(ubx_lo);
            IntVect uhi_bx_yz_hi(ubx_hi);

            // Fix lo and hi limits
            uhi_bx_yz_lo[0] = dom_hi[0]+1;
            uhi_bx_yz_hi[0] = dom_hi[0]+1;

            const Box uhi_bx_yz(uhi_bx_yz_lo, uhi_bx_yz_hi);

            reduce_op.eval(uhi_bx_yz, reduce_data,
              [bct_ihi, dom_hi, dydz, minf, pinf, pout, flux_x_arr, apx, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

                if(bct == pout) {
                  m_out += apx(i,j,k)*flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                else if ((bct == minf) || (bct == pinf)) {
                  m_in  -= apx(i,j,k)*flux_x_arr(i,j,k,n+scomp)*dydz;
                }
                return {m_in, m_out};
              });

          }


          if (nbot == 0) { // dom_lo[1] == vbx_lo[1]

              // Create InVects for following Box
              IntVect vlo_bx_xz_lo(vbx_lo);
              IntVect vlo_bx_xz_hi(vbx_hi);

              // Fix lo and hi limits
              vlo_bx_xz_lo[1] = dom_lo[1];
              vlo_bx_xz_hi[1] = dom_lo[1];

              const Box vlo_bx_xz(vlo_bx_xz_lo, vlo_bx_xz_hi);
              reduce_op.eval(vlo_bx_xz, reduce_data,
                [bct_jlo, dom_lo, dxdz, minf, pinf, pout, flux_y_arr, apy, n, scomp]
                AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {

                  Real m_in  = 0.;
                  Real m_out = 0.;

                  const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

                  if(bct == pout) {
                    m_out -= apy(i,j,k)*flux_y_arr(i,j,k,n+scomp)*dxdz;
                  }
                  else if ((bct == minf) || (bct == pinf)) {
                    m_in  += apy(i,j,k)*flux_y_arr(i,j,k,n+scomp)*dxdz;
                  }
                  return {m_in, m_out};
                });

          } // nbot

          if (ntop == 1) {

            // Create InVects for following Box
            IntVect vhi_bx_xz_lo(vbx_lo);
            IntVect vhi_bx_xz_hi(vbx_hi);

            // Fix lo and hi limits
            vhi_bx_xz_lo[1] = dom_hi[1]+1;
            vhi_bx_xz_hi[1] = dom_hi[1]+1;

            const Box vhi_bx_xz(vhi_bx_xz_lo, vhi_bx_xz_hi);

            reduce_op.eval(vhi_bx_xz, reduce_data,
              [bct_jhi, dom_hi, dxdz, minf, pinf, pout, flux_y_arr, apy, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {

                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

                if(bct == pout) {
                  m_out += apy(i,j,k)*flux_y_arr(i,j,k,n+scomp)*dxdz;
                }
                else if ((bct == minf) || (bct == pinf)) {
                  m_in  -= apy(i,j,k)*flux_y_arr(i,j,k,n+scomp)*dxdz;
                }
                return {m_in, m_out};
              });

          } // ntop

          if (ndwn == 0) {

            // Create InVects for following Boxes
            IntVect wlo_bx_xy_lo(wbx_lo);
            IntVect wlo_bx_xy_hi(wbx_hi);

            // Fix lo and hi limits
            wlo_bx_xy_lo[2] = dom_lo[2];
            wlo_bx_xy_hi[2] = dom_lo[2];

            const Box wlo_bx_xy(wlo_bx_xy_lo, wlo_bx_xy_hi);

            reduce_op.eval(wlo_bx_xy, reduce_data,
              [bct_klo, dom_lo, dxdy, minf, pinf, pout, flux_z_arr, apz, n, scomp]
              AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
              {
                Real m_in  = 0.;
                Real m_out = 0.;

                const int bct = bct_klo(i,j,dom_lo[2]-1,0);

                if(bct == pout) {
                  m_out -= apz(i,j,k)*flux_z_arr(i,j,k,n+scomp)*dxdy;
                }
                else if ((bct == minf) || (bct == pinf)) {
                  m_in  += apz(i,j,k)*flux_z_arr(i,j,k,n+scomp)*dxdy;
                }
                return {m_in, m_out};
              });
          } // ndwn

          if (nup == 0) {

            // Create InVects for following Boxes
            IntVect whi_bx_xy_lo(wbx_lo);
            IntVect whi_bx_xy_hi(wbx_hi);

            // Fix lo and hi limits
            whi_bx_xy_lo[2] = dom_hi[2]+1;
            whi_bx_xy_hi[2] = dom_hi[2]+1;

            const Box whi_bx_xy(whi_bx_xy_lo, whi_bx_xy_hi);

            reduce_op.eval(whi_bx_xy, reduce_data,
               [bct_khi, dom_hi, dxdy, minf, pinf, pout, flux_z_arr, apz, n, scomp]
               AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
               {
                 Real m_in  = 0.;
                 Real m_out = 0.;

                 const int bct = bct_khi(i,j,dom_hi[2]+1,0);

                 if(bct == pout) {
                   m_out += apz(i,j,k)*flux_z_arr(i,j,k,n+scomp)*dxdy;
                 }
                 else if ((bct == minf) || (bct == pinf)) {
                   m_in  -= apz(i,j,k)*flux_z_arr(i,j,k,n+scomp)*dxdy;
                 }
                 return {m_in, m_out};
               });
          } // nup


        }


      }

      ReduceTuple host_tuple = reduce_data.value(reduce_op);

      mass_in[n]  += amrex::get<0>(host_tuple);
      mass_out[n] += amrex::get<1>(host_tuple);


    } // Loop over species
  } // Loop over levels

  // Copy to a single array for one reduce call.
  std::vector<Real> mass_flow(2*nspecies_g, 0.);
  for (int n=0; n < nspecies_g; ++n) {
    mass_flow[n] = mass_in[n];
    mass_flow[n+nspecies_g] = mass_out[n];
  }

  ParallelDescriptor::ReduceRealSum(mass_flow.data(), 2*nspecies_g);

  // Copy into global variables.
  for (int n=0; n < nspecies_g; ++n) {
    mass_inflow[n] = mass_flow[n];
    mass_outflow[n] = mass_flow[n+nspecies_g];
  }
}
