#include <mfix.H>

#include <mfix_fluid.H>

using namespace amrex;

//
// Set the BCs for density only
//
void
mfix::mfix_set_enthalpy_bcs (Real time,
                             Vector< MultiFab* > const& h_g_in)
{
  BL_PROFILE("mfix::mfix_set_enthalpy_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
    Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*h_g_in[lev], false); mfi.isValid(); ++mfi)
      set_enthalpy_bcs(time, lev, (*h_g_in[lev])[mfi], domain);

    h_g_in[lev]->FillBoundary(geom[lev].periodicity());

    EB_set_covered(*h_g_in[lev], 0, h_g_in[lev]->nComp(), h_g_in[lev]->nGrow(),
        covered_val);
  }
}

void
mfix::set_enthalpy_bcs (Real time,
                        const int lev,
                        FArrayBox& h_g_fab,
                        const Box& domain)
{
  BL_PROFILE("mfix::set_enthalpy_bcs()");

  const EBFArrayBox& hg_EB_fab = static_cast<EBFArrayBox const&>(h_g_fab);
  const EBCellFlagFab& flags = hg_EB_fab.getEBCellFlagFab();

  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<const int> const& bct_ilo = bc_list.bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_list.bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_list.bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_list.bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_list.bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_list.bc_khi[lev]->array();

  const int nspecies_g = fluid.nspecies();

  // Flag to understand if fluid is a mixture
  const int fluid_is_a_mixture = fluid.isMixture();

  Array4<Real> const& h_g = h_g_fab.array();

  IntVect h_g_lo(h_g_fab.loVect());
  IntVect h_g_hi(h_g_fab.hiVect());

  const int nlft = amrex::max(0, dom_lo[0]-h_g_lo[0]);
  const int nbot = amrex::max(0, dom_lo[1]-h_g_lo[1]);
  const int ndwn = amrex::max(0, dom_lo[2]-h_g_lo[2]);

  const int nrgt = amrex::max(0, h_g_hi[0]-dom_hi[0]);
  const int ntop = amrex::max(0, h_g_hi[1]-dom_hi[1]);
  const int nup  = amrex::max(0, h_g_hi[2]-dom_hi[2]);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(h_g_hi), bx_xz_lo_hi_3D(h_g_hi), bx_xy_lo_hi_3D(h_g_hi);
  IntVect bx_yz_hi_lo_3D(h_g_lo), bx_xz_hi_lo_3D(h_g_lo), bx_xy_hi_lo_3D(h_g_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for GPU loops
  const Box bx_yz_lo_3D(h_g_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, h_g_hi);

  const Box bx_xz_lo_3D(h_g_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, h_g_hi);

  const Box bx_xy_lo_3D(h_g_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, h_g_hi);

  // Update temperature before using to update enthalpy
  set_temperature_bc_values (time);
  Real* p_bc_t_g = m_bc_t_g.data();

  // Update species boundary values
  if(fluid_is_a_mixture)
    set_species_bc_values(time);

  Real** p_bc_X_gk = fluid_is_a_mixture ? m_bc_X_gk_ptr.data() : nullptr;

  const auto& fluid_parms = fluid.parameters();

  auto const& flags_arr = flags.const_array();

  auto set_enthalpy_bcs_in_box = [h_g,p_bc_t_g,fluid_is_a_mixture,
       p_bc_X_gk,nspecies_g,fluid_parms,flags_arr]
    AMREX_GPU_DEVICE (int bct, int bcv, IntVect dom_ijk, int i, int j, int k) noexcept
  {
    const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

    if(bct == BCList::pout) {
      h_g(i,j,k) = h_g(dom_ijk);
    }
    else if (bct == BCList::minf || bct == BCList::pinf) {
      if (!fluid_is_a_mixture) {
        h_g(i,j,k) = fluid_parms.calc_h_g<run_on>(p_bc_t_g[bcv], cell_is_covered);
      }
      else {
        Real h_g_sum(0);

        for (int n(0); n < nspecies_g; n++) {
          const Real h_gk = fluid_parms.calc_h_gk<run_on>(p_bc_t_g[bcv], n, cell_is_covered);

          h_g_sum += p_bc_X_gk[n][bcv]*h_gk;
        }

        h_g(i,j,k) = h_g_sum;
      }
    }
  };

  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D, [bct_ilo,dom_lo,set_enthalpy_bcs_in_box]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const IntVect dom_ijk(dom_lo[0],j,k);

      set_enthalpy_bcs_in_box(bct, bcv, dom_ijk, i, j, k);
    });
  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D, [bct_ihi,dom_hi,set_enthalpy_bcs_in_box]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const IntVect dom_ijk(dom_hi[0],j,k);

      set_enthalpy_bcs_in_box(bct, bcv, dom_ijk, i, j, k);
    });
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D, [bct_jlo,dom_lo,set_enthalpy_bcs_in_box]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const IntVect dom_ijk(i,dom_lo[1],k);

      set_enthalpy_bcs_in_box(bct, bcv, dom_ijk, i, j, k);
    });
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D, [bct_jhi,dom_hi,set_enthalpy_bcs_in_box]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const IntVect dom_ijk(i,dom_hi[1],k);

      set_enthalpy_bcs_in_box(bct, bcv, dom_ijk, i, j, k);
    });
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D, [bct_klo,dom_lo,set_enthalpy_bcs_in_box]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const IntVect dom_ijk(i,j,dom_lo[2]);

      set_enthalpy_bcs_in_box(bct, bcv, dom_ijk, i, j, k);
    });
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D, [bct_khi,dom_hi,set_enthalpy_bcs_in_box]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const IntVect dom_ijk(i,j,dom_hi[2]);

      set_enthalpy_bcs_in_box(bct, bcv, dom_ijk, i, j, k);
    });
  }

  //Gpu::synchronize();
}
