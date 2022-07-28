#include <mfix.H>

#include <mfix_fluid.H>

using namespace amrex;

//
// Set the BCs for density only
//
void
mfix::mfix_set_density_bcs (Real time,
                            Vector< MultiFab* > const& ro_g_in)
{
  BL_PROFILE("mfix::mfix_set_density_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*(m_leveldata[lev]->ep_g), false); mfi.isValid(); ++mfi)
     {
        set_density_bcs(time, lev, (*ro_g_in[lev])[mfi], domain);
     }

     ro_g_in[lev]->FillBoundary(geom[lev].periodicity());

     EB_set_covered(*ro_g_in[lev], 0, ro_g_in[lev]->nComp(), ro_g_in[lev]->nGrow(), covered_val);
  }
}

void
mfix::set_density_bcs (Real time,
                       const int lev,
                       FArrayBox& scal_fab,
                       const Box& domain)

{
  BL_PROFILE("mfix::set_density_bcs()");

  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<const int> const& bct_ilo = bc_list.bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_list.bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_list.bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_list.bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_list.bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_list.bc_khi[lev]->array();

  Array4<Real> const& scal_arr = scal_fab.array();

  IntVect scal_lo(scal_fab.loVect());
  IntVect scal_hi(scal_fab.hiVect());

  const int nlft = amrex::max(0, dom_lo[0]-scal_lo[0]);
  const int nbot = amrex::max(0, dom_lo[1]-scal_lo[1]);
  const int ndwn = amrex::max(0, dom_lo[2]-scal_lo[2]);

  const int nrgt = amrex::max(0, scal_hi[0]-dom_hi[0]);
  const int ntop = amrex::max(0, scal_hi[1]-dom_hi[1]);
  const int nup  = amrex::max(0, scal_hi[2]-dom_hi[2]);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(scal_hi), bx_xz_lo_hi_3D(scal_hi), bx_xy_lo_hi_3D(scal_hi);
  IntVect bx_yz_hi_lo_3D(scal_lo), bx_xz_hi_lo_3D(scal_lo), bx_xy_hi_lo_3D(scal_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for GPU loops
  const Box bx_yz_lo_3D(scal_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, scal_hi);

  const Box bx_xz_lo_3D(scal_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, scal_hi);

  const Box bx_xy_lo_3D(scal_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, scal_hi);

  // The temperature will be useful in the future
  // set_temperature_bc_values (time);
  // Real* p_bc_t_g  = m_bc_t_g.data();
  set_density_bc_values(time);
  Real* p_bc_ro_g  = m_boundary_conditions.bc_ro_g().data();

  auto set_density_bcs_in_box = [scal_arr,p_bc_ro_g]
  AMREX_GPU_DEVICE (int bct, int bcv, IntVect dom_ijk, int i, int j, int k) noexcept
  {
    if ((bct == BCList::pinf) || (bct == BCList::pout))
      scal_arr(i,j,k) = scal_arr(dom_ijk);
    else if (bct == BCList::minf)
      scal_arr(i,j,k) = p_bc_ro_g[bcv];
  };


  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D, [bct_ilo,dom_lo,set_density_bcs_in_box]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const IntVect dom_ijk(dom_lo[0],j,k);

      set_density_bcs_in_box(bct, bcv, dom_ijk, i, j, k);
    });
  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D, [bct_ihi,dom_hi,set_density_bcs_in_box]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const IntVect dom_ijk(dom_hi[0],j,k);

      set_density_bcs_in_box(bct, bcv, dom_ijk, i, j, k);
    });
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D, [bct_jlo,dom_lo,set_density_bcs_in_box]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const IntVect dom_ijk(i,dom_lo[1],k);

      set_density_bcs_in_box(bct, bcv, dom_ijk, i, j, k);
    });
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D, [bct_jhi,dom_hi,set_density_bcs_in_box]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const IntVect dom_ijk(i,dom_hi[1],k);

      set_density_bcs_in_box(bct, bcv, dom_ijk, i, j, k);
    });
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D, [bct_klo,dom_lo,set_density_bcs_in_box]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const IntVect dom_ijk(i,j,dom_lo[2]);

      set_density_bcs_in_box(bct, bcv, dom_ijk, i, j, k);
    });
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D, [bct_khi,dom_hi,set_density_bcs_in_box]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const IntVect dom_ijk(i,j,dom_hi[2]);

      set_density_bcs_in_box(bct, bcv, dom_ijk, i, j, k);

    });
  }


}
