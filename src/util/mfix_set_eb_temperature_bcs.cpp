#include <mfix.H>

#include <mfix_calc_cell.H>
#include <mfix_bc_parms.H>
#include <mfix_fluid_parms.H>
#include <mfix_bc_list.H>

using namespace amrex;

//
// Set the BCs for temperature only
//
void
mfix::mfix_set_eb_temperature_bcs (Vector< MultiFab* > const& eb_T_g_in,
                                   Vector< MultiFab* > const& eb_k_g_in)
{
  BL_PROFILE("mfix::mfix_set_eb_temperature_bcs()");

  Vector< MultiFab* > k_g = get_k_g();

  for (int lev = 0; lev < nlev; lev++)
  {
     Box domain(geom[lev].Domain());

     // NOTE: setting T_g and k_g equal to zero should give us conditions
     // equivalent to homogeneous Neumann. To set inhomogeneous Dirichlet we
     // need to call set_eb_temperature_dirichlet_values
     eb_T_g_in[lev]->setVal(0.);
     eb_k_g_in[lev]->setVal(0.);

     const auto& factory =
       dynamic_cast<EBFArrayBoxFactory const&>(eb_T_g_in[lev]->Factory());

     const auto& flags = factory.getMultiEBCellFlagFab();

     Real dx = geom[lev].CellSize(0);
     Real dy = geom[lev].CellSize(1);
     Real dz = geom[lev].CellSize(2);

     const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();    

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*eb_T_g_in[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       const Box& bx = mfi.tilebox();
       const Box& sbx = (*eb_T_g_in[lev])[mfi].box();

       FabType t = flags[mfi].getType(bx);

       // We update only cut-cells values
       if (t == FabType::singlevalued) {
         set_eb_temperature_bcs(sbx, bx, domain, dx, dy, dz, plo,
             (*eb_T_g_in[lev])[mfi], (*eb_k_g_in[lev])[mfi], (*k_g[lev])[mfi],
             flags[mfi]);
       }
     }

     // Do this after as well as before to pick up terms that got updated in the
     // call above
     eb_T_g_in[lev]->FillBoundary(geom[lev].periodicity());
     eb_k_g_in[lev]->FillBoundary(geom[lev].periodicity());
  }
}

void
mfix::set_eb_temperature_bcs (const Box& sbx,
                              const Box& bx,
                              const Box& domain,
                              const Real dx,
                              const Real dy,
                              const Real dz,
                              const GpuArray<Real, 3>& plo,
                              FArrayBox& eb_T_g_fab,
                              FArrayBox& eb_k_g_fab,
                              FArrayBox& k_g_fab,
                              const EBCellFlagFab& flags_fab)
{
  BL_PROFILE("mfix::set_eb_temperature_bcs()");

  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  // EB arrays
  Array4<Real> const& eb_T_g = eb_T_g_fab.array();
  Array4<Real> const& eb_k_g = eb_k_g_fab.array();

  //// Problem variable
  //Array4<const Real> const& k_g = k_g_fab.array();

  // TODO this has to be replaced by the above k_g
  const Real k_g = FLUID::k_g0;

  // Flags
  Array4<const EBCellFlag> const& flags = flags_fab.array();

  IntVect eb_T_g_lo(eb_T_g_fab.loVect());
  IntVect eb_T_g_hi(eb_T_g_fab.hiVect());

  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  const int eb = bc_list.get_eb();

  // Set the initial conditions.
  for(int bcv(0); bcv < BC::bc.size(); ++bcv)
  {
    if ( BC::bc[bcv].type == eb ) {
      int i_w(0), j_s(0), k_b(0);
      int i_e(0), j_n(0), k_t(0);

      calc_cell_ic(dx, dy, dz,
                   BC::bc[bcv].region->lo(),
                   BC::bc[bcv].region->hi(),
                   plo.data(),
                   i_w, i_e, j_s, j_n, k_b, k_t);

      // Use the volume fraction already calculated from particle data
      const Real eb_temperature = BC::bc[bcv].eb.temperature;

      const int istart = amrex::max(slo[0], i_w);
      const int jstart = amrex::max(slo[1], j_s);
      const int kstart = amrex::max(slo[2], k_b);
      const int iend   = amrex::min(shi[0], i_e);
      const int jend   = amrex::min(shi[1], j_n);
      const int kend   = amrex::min(shi[2], k_t);

      {
        const int first = (slo[0] < dom_lo[0] and dom_lo[0] == istart) ? slo[0] : istart;
        const int last  = (shi[0] > dom_hi[0] and dom_hi[0] == iend)   ? shi[0] : iend;

        const Box bx(IntVect(first, jstart, kstart),
                     IntVect(last,  jend,   kend));

        ParallelFor(bx, [flags,eb_temperature,eb_T_g,k_g,eb_k_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (flags(i,j,k).isSingleValued()) {
            eb_T_g(i,j,k) = eb_temperature;
            //eb_k_g(i,j,k) = k_g(i,j,k);
            eb_k_g(i,j,k) = k_g;
          }
        });
      }

      {
        const int first = (slo[1] < dom_lo[1] and dom_lo[1] == jstart) ? slo[1] : jstart;
        const int last  = (shi[1] > dom_hi[1] and dom_hi[1] == jend)   ? shi[1] : jend;

        const Box bx(IntVect(istart, first, kstart),
                     IntVect(iend,   last,  kend));

        ParallelFor(bx, [flags,eb_temperature,eb_T_g,k_g,eb_k_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (flags(i,j,k).isSingleValued()) {
            eb_T_g(i,j,k) = eb_temperature;
            //eb_k_g(i,j,k) = k_g(i,j,k);
            eb_k_g(i,j,k) = k_g;
          }
        });
      }

      {
        const int first = (slo[2] < dom_lo[2] and dom_lo[2] == kstart) ? slo[2] : kstart;
        const int last  = (shi[2] > dom_hi[2] and dom_hi[2] == kend)   ? shi[2] : kend;

        const Box bx(IntVect(istart, jstart, first),
                     IntVect(iend,   jend,   last));

        ParallelFor(bx, [flags,eb_temperature,eb_T_g,k_g,eb_k_g]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (flags(i,j,k).isSingleValued()) {
            eb_T_g(i,j,k) = eb_temperature;
            //eb_k_g(i,j,k) = k_g(i,j,k);
            eb_k_g(i,j,k) = k_g;
          }
        });
      }

    }
  }
}
