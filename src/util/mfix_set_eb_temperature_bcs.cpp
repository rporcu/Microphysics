#include <mfix.H>

#include <mfix_calc_cell.H>
#include <mfix_bc.H>
#include <mfix_fluid.H>

using namespace amrex;

void
MFIXBoundaryConditions::set_eb_temperature_bcs (Vector< MultiFab* > const& eb_T_g_in)
{
  BL_PROFILE("MFIXBoundaryConditions::set_eb_temperature_bcs()");

  const int nlev = eb_T_g_in.size();

  for (int lev = 0; lev < nlev; lev++)
  {
    Box domain(m_geom[lev].Domain());

    // NOTE: setting T_g and k_g equal to zero should give us conditions
    // equivalent to homogeneous Neumann. To set inhomogeneous Dirichlet we
    // need to call set_eb_temperature_dirichlet_values
    eb_T_g_in[lev]->setVal(0.);

    const auto& factory =
      dynamic_cast<EBFArrayBoxFactory const&>(eb_T_g_in[lev]->Factory());

    const auto& flags = factory.getMultiEBCellFlagFab();

    Real dx = m_geom[lev].CellSize(0);
    Real dy = m_geom[lev].CellSize(1);
    Real dz = m_geom[lev].CellSize(2);

    const GpuArray<Real, 3> plo = m_geom[lev].ProbLoArray();    

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*eb_T_g_in[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();
      const Box& sbx = (*eb_T_g_in[lev])[mfi].box();

      // We update only cut-cells values
      if (flags[mfi].getType(bx) == FabType::singlevalued) {

        FArrayBox& eb_T_g_fab = (*eb_T_g_in[lev])[mfi];

        IntVect dom_lo(domain.loVect());
        IntVect dom_hi(domain.hiVect());

        // EB arrays
        Array4<Real> const& eb_T_g = eb_T_g_fab.array();

        // Flags
        Array4<const EBCellFlag> const& flags_arr = flags.const_array(mfi);

        IntVect eb_T_g_lo(eb_T_g_fab.loVect());
        IntVect eb_T_g_hi(eb_T_g_fab.hiVect());

        const IntVect slo(sbx.loVect());
        const IntVect shi(sbx.hiVect());

        const IntVect domlo(domain.loVect());
        const IntVect domhi(domain.hiVect());

        // Set the initial conditions.
        for(int bcv(0); bcv < m_bc.size(); ++bcv)
        {
          if (m_bc[bcv].type == BCList::eb) {
            int i_w(0), j_s(0), k_b(0);
            int i_e(0), j_n(0), k_t(0);

            calc_cell_ic(dx, dy, dz,
                         m_bc[bcv].region->lo(),
                         m_bc[bcv].region->hi(),
                         plo.data(),
                         i_w, i_e, j_s, j_n, k_b, k_t);

            // Use the volume fraction already calculated from particle data
            const Real eb_temperature = m_bc[bcv].eb.temperature;

            const int istart = amrex::max(slo[0], i_w);
            const int jstart = amrex::max(slo[1], j_s);
            const int kstart = amrex::max(slo[2], k_b);
            const int iend   = amrex::min(shi[0], i_e);
            const int jend   = amrex::min(shi[1], j_n);
            const int kend   = amrex::min(shi[2], k_t);

            {
              const int first = (slo[0] < dom_lo[0] && dom_lo[0] == istart) ? slo[0] : istart;
              const int last  = (shi[0] > dom_hi[0] && dom_hi[0] == iend)   ? shi[0] : iend;

              const Box local_bx(IntVect(first, jstart, kstart),
                                 IntVect(last,  jend,   kend));

              ParallelFor(local_bx, [flags_arr,eb_temperature,eb_T_g]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                if (flags_arr(i,j,k).isSingleValued()) {
                  eb_T_g(i,j,k) = eb_temperature;
                }
              });
            }

            {
              const int first = (slo[1] < dom_lo[1] && dom_lo[1] == jstart) ? slo[1] : jstart;
              const int last  = (shi[1] > dom_hi[1] && dom_hi[1] == jend)   ? shi[1] : jend;

              const Box local_bx(IntVect(istart, first, kstart),
                                 IntVect(iend,   last,  kend));

              ParallelFor(local_bx, [flags_arr,eb_temperature,eb_T_g]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                if (flags_arr(i,j,k).isSingleValued()) {
                  eb_T_g(i,j,k) = eb_temperature;
                }
              });
            }

            {
              const int first = (slo[2] < dom_lo[2] && dom_lo[2] == kstart) ? slo[2] : kstart;
              const int last  = (shi[2] > dom_hi[2] && dom_hi[2] == kend)   ? shi[2] : kend;

              const Box local_bx(IntVect(istart, jstart, first),
                                 IntVect(iend,   jend,   last));

              ParallelFor(local_bx, [flags_arr,eb_temperature,eb_T_g]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                if (flags_arr(i,j,k).isSingleValued()) {
                  eb_T_g(i,j,k) = eb_temperature;
                }
              });
            }
          }
        }
      }

    } // end MFIter loop

    // Do this after as well as before to pick up terms that got updated in the
    // call above
    eb_T_g_in[lev]->FillBoundary(m_geom[lev].periodicity());
  }
}
