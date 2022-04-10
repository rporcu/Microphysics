#include <mfix.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>
#include <mfix_mf_helpers.H>

void
mfix::set_temperature_bc0 (const Box& sbx,
                           MFIter* mfi,
                           const int lev,
                           const Box& domain)
{
  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

  auto& ld = *m_leveldata[lev];

  const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*ld.ep_g)[*mfi]);
  const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

  const int nspecies_g = fluid.nspecies;

  // Flag to understand if fluid is a mixture
  const int fluid_is_a_mixture = fluid.is_a_mixture;

  Real* p_bc_t_g = m_bc_t_g.data();
  Real** p_bc_X_gk = fluid_is_a_mixture ? m_bc_X_gk_ptr.data() : nullptr;

  Array4<Real> const& a_T_g  = m_leveldata[lev]->T_g->array(*mfi);
  Array4<Real> const& a_h_g  = m_leveldata[lev]->h_g->array(*mfi);

  const IntVect sbx_lo(sbx.loVect());
  const IntVect sbx_hi(sbx.hiVect());

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  Array4<const int> const& a_bc_ilo = bc_list.bc_ilo[lev]->array();
  Array4<const int> const& a_bc_ihi = bc_list.bc_ihi[lev]->array();
  Array4<const int> const& a_bc_jlo = bc_list.bc_jlo[lev]->array();
  Array4<const int> const& a_bc_jhi = bc_list.bc_jhi[lev]->array();
  Array4<const int> const& a_bc_klo = bc_list.bc_klo[lev]->array();
  Array4<const int> const& a_bc_khi = bc_list.bc_khi[lev]->array();

  const int nlft = amrex::max(0,dom_lo[0]-sbx_lo[0]);
  const int nbot = amrex::max(0,dom_lo[1]-sbx_lo[1]);
  const int ndwn = amrex::max(0,dom_lo[2]-sbx_lo[2]);

  const int nrgt = amrex::max(0,sbx_hi[0]-dom_hi[0]);
  const int ntop = amrex::max(0,sbx_hi[1]-dom_hi[1]);
  const int nup  = amrex::max(0,sbx_hi[2]-dom_hi[2]);

  auto& fluid_parms = *fluid.parameters;

  auto const& flags_arr = flags.const_array();

  auto set_temperature_bc0_in_box = [a_T_g,p_bc_t_g,
      a_h_g,fluid_is_a_mixture,nspecies_g,p_bc_X_gk,fluid_parms,run_on_device,
      flags_arr]
  AMREX_GPU_DEVICE (int bct, int bcv, int i, int j, int k) noexcept
  {
    const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

    if((bct == BCList::pinf) || (bct == BCList::pout) || (bct == BCList::minf))
    {
      a_T_g(i,j,k)  = p_bc_t_g[bcv];

      if (!fluid_is_a_mixture) {
        a_h_g(i,j,k) = run_on_device ?
          fluid_parms.calc_h_g<RunOn::Device>(p_bc_t_g[bcv], cell_is_covered) :
          fluid_parms.calc_h_g<RunOn::Host>(p_bc_t_g[bcv], cell_is_covered);
      }
      else {
        Real h_g_sum(0);

        for (int n(0); n < nspecies_g; n++) {
          const Real h_gk = run_on_device ?
            fluid_parms.calc_h_gk<RunOn::Device>(p_bc_t_g[bcv], n, cell_is_covered) :
            fluid_parms.calc_h_gk<RunOn::Host>(p_bc_t_g[bcv], n, cell_is_covered);

          h_g_sum += p_bc_X_gk[n][bcv]*h_gk;
        }

        a_h_g(i,j,k)  = h_g_sum;
      }
    }
  };

  if (nlft > 0)
  {
    IntVect bx_yz_lo_hi_3D(sbx_hi);

    // Fix lo and hi limits
    bx_yz_lo_hi_3D[0] = dom_lo[0]-1;

    const Box bx_yz_lo_3D(sbx_lo, bx_yz_lo_hi_3D);

    ParallelFor(bx_yz_lo_3D, [dom_lo,a_bc_ilo,set_temperature_bc0_in_box]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_ilo(dom_lo[0]-1,j,k,1);
      const int bct = a_bc_ilo(dom_lo[0]-1,j,k,0);

      set_temperature_bc0_in_box(bct, bcv, i, j, k);
    });
  }

  if (nrgt > 0)
  {
    IntVect bx_yz_hi_lo_3D(sbx_lo);

    // Fix lo and hi limits
    bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

    const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, sbx_hi);

    ParallelFor(bx_yz_hi_3D, [dom_hi,a_bc_ihi,set_temperature_bc0_in_box]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_ihi(dom_hi[0]+1,j,k,1);
      const int bct = a_bc_ihi(dom_hi[0]+1,j,k,0);

      set_temperature_bc0_in_box(bct, bcv, i, j, k);
    });
  }

  if (nbot > 0)
  {
    IntVect bx_xz_lo_hi_3D(sbx_hi);

    // Fix lo and hi limits
    bx_xz_lo_hi_3D[1] = dom_lo[1]-1;

    const Box bx_xz_lo_3D(sbx_lo, bx_xz_lo_hi_3D);

    ParallelFor(bx_xz_lo_3D, [dom_lo,a_bc_jlo,set_temperature_bc0_in_box]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_jlo(i,dom_lo[1]-1,k,1);
      const int bct = a_bc_jlo(i,dom_lo[1]-1,k,0);

      set_temperature_bc0_in_box(bct, bcv, i, j, k);
    });
  }

  if (ntop > 0)
  {
    IntVect bx_xz_hi_lo_3D(sbx_lo);

    // Fix lo and hi limits
    bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

    const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, sbx_hi);

    ParallelFor(bx_xz_hi_3D, [dom_hi,a_bc_jhi,set_temperature_bc0_in_box]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_jhi(i,dom_hi[1]+1,k,1);
      const int bct = a_bc_jhi(i,dom_hi[1]+1,k,0);

      set_temperature_bc0_in_box(bct, bcv, i, j, k);
    });
  }

  if (ndwn > 0)
  {
    IntVect bx_xy_lo_hi_3D(sbx_hi);

    // Fix lo and hi limits
    bx_xy_lo_hi_3D[2] = dom_lo[2]-1;

    const Box bx_xy_lo_3D(sbx_lo, bx_xy_lo_hi_3D);

    ParallelFor(bx_xy_lo_3D, [dom_lo,a_bc_klo,set_temperature_bc0_in_box]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_klo(i,j,dom_lo[2]-1,1);
      const int bct = a_bc_klo(i,j,dom_lo[2]-1,0);

      set_temperature_bc0_in_box(bct, bcv, i, j, k);
    });
  }

  if (nup > 0)
  {
    IntVect bx_xy_hi_lo_3D(sbx_lo);

    // Fix lo and hi limits
    bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

    const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, sbx_hi);

    ParallelFor(bx_xy_hi_3D, [dom_hi,a_bc_khi,set_temperature_bc0_in_box]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bcv = a_bc_khi(i,j,dom_hi[2]+1,1);
      const int bct = a_bc_khi(i,j,dom_hi[2]+1,0);

      set_temperature_bc0_in_box(bct, bcv, i, j, k);
    });
  }

  //Gpu::synchronize();
}
