#include <mfix.H>

#include <mfix_calc_cell.H>
#include <mfix_bc.H>
#include <mfix_fluid.H>

#include <AMReX_ParallelReduce.H>

using namespace amrex;

void
MFIXBoundaryConditions::calc_eb_bc_areas (MFIXEmbeddedBoundaries& embedded_boundaries,
                                          Vector< MultiFab* > const& eb_mf)
{
  BL_PROFILE("MFIXBoundaryConditions::calc_eb_bc_areas()");

  const int nlev = eb_mf.size();

  embedded_boundaries.set_compute_area(0);

  // Number of boundary conditions
  const int nbc = m_bc.size();

  amrex::Vector<Real> eb_areas(nbc, Real(0.));

  for (int lev = 0; lev < nlev; lev++) {

    const GpuArray<Real,3> dx = m_geom[lev].CellSizeArray();
    const Real dx2 = dx[0]*dx[0];

    const auto& factory =
      dynamic_cast<EBFArrayBoxFactory const&>(eb_mf[lev]->Factory());

    const auto& flags = factory.getMultiEBCellFlagFab();

    for(int bcv(0); bcv < nbc; ++bcv) {

      if (m_bc[bcv].type != BCList::eb) {

        eb_areas[bcv] = Real(0.0);

      } else {

        const int  has_normal = m_bc[bcv].eb.has_normal;
        amrex::GpuArray<amrex::Real,3> normal{0.};
        if (has_normal) {
           normal[0] = m_bc[bcv].eb.normal[0];
           normal[1] = m_bc[bcv].eb.normal[1];
           normal[2] = m_bc[bcv].eb.normal[2];
        }

        const Real pad = std::numeric_limits<float>::epsilon();
        const Real normal_tol = m_bc[bcv].eb.normal_tol;

        const Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
        const Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        const Box ic_bx = calc_ic_box(m_geom[lev], m_bc[bcv].region);

        for (MFIter mfi(*eb_mf[lev], false); mfi.isValid(); ++mfi) {

          const Box& bx = mfi.tilebox();
          FabType t = flags[mfi].getType(bx);

          if (t == FabType::singlevalued && ic_bx.intersects(bx)) {

            Array4<EBCellFlag const> const& flagsfab = flags.const_array(mfi);

            // Area fractions
            Array4<Real const> const& apx = factory.getAreaFrac()[0]->const_array(mfi);
            Array4<Real const> const& apy = factory.getAreaFrac()[1]->const_array(mfi);
            Array4<Real const> const& apz = factory.getAreaFrac()[2]->const_array(mfi);

            // Area of eb face
            Array4<Real const> const& barea = factory.getBndryArea().const_array(mfi);

            const Box bx_int = bx&(ic_bx);

            reduce_op.eval(bx_int, reduce_data, [flagsfab, barea, has_normal, normal,
              norm_tol_lo, norm_tol_hi, apx, apy, apz]
            AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {

              Real area(0.0);

              if(flagsfab(i,j,k).isSingleValued()) {

                area = barea(i,j,k);

                if(has_normal) {

                     Real apxm = apx(i  ,j  ,k  );
                     Real apxp = apx(i+1,j  ,k  );
                     Real apym = apy(i  ,j  ,k  );
                     Real apyp = apy(i  ,j+1,k  );
                     Real apzm = apz(i  ,j  ,k  );
                     Real apzp = apz(i  ,j  ,k+1);

                     Real dapx = apxm-apxp;
                     Real dapy = apym-apyp;
                     Real dapz = apzm-apzp;
                     Real anorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
                     Real anorminv = 1.0/anorm;

                     Real anrmx = dapx * anorminv;
                     Real anrmy = dapy * anorminv;
                     Real anrmz = dapz * anorminv;

                     const Real dotprod = anrmx*normal[0] + anrmy*normal[1] + anrmz*normal[2];

                     area *= (norm_tol_lo <= dotprod) ? Real(1.0) : Real(0.0);
                     area *= (dotprod <= norm_tol_hi) ? Real(1.0) : Real(0.0);
                }
              }
              return {area};
            });
          }
        } // loop over MFIter

        ReduceTuple host_tuple = reduce_data.value(reduce_op);
        eb_areas[bcv] = dx2*amrex::get<0>(host_tuple);

      }
    } // Loop over boundary conditions

    ParallelAllReduce::Sum(eb_areas.data(), eb_areas.size(), ParallelContext::CommunicatorAll());

    for(int bcv(0); bcv < nbc; ++bcv) {

      BC_t& bc = m_bc[bcv];

      if (bc.type == BCList::eb) {
        bc.eb.area = eb_areas[bcv];
      }
    }

  } //nlev
}
