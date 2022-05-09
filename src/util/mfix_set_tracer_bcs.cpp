#include <mfix.H>

#include <mfix_fluid_parms.H>

using namespace amrex;

//
// Set the BCs for all the variables EXCEPT pressure or velocity.
//
void
mfix::mfix_set_tracer_bcs (Real time,
                           Vector< MultiFab* > const& trac_in)
{
  BL_PROFILE("mfix::mfix_set_tracer_bcs()");

  if (fluid.solve_tracer)
  {
    for (int lev = 0; lev < nlev; lev++)
    {
      Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*(m_leveldata[lev]->ep_g), false); mfi.isValid(); ++mfi)
      {
        set_tracer_bcs(time, lev, (*trac_in[lev])[mfi], domain);
      }

      trac_in[lev] -> FillBoundary (geom[lev].periodicity());

      EB_set_covered(*trac_in[lev], 0, trac_in[lev]->nComp(), trac_in[lev]->nGrow(), covered_val);
    }
  }
}


void
mfix::set_tracer_bcs (Real /*time*/,
                      const int lev,
                      FArrayBox& scal_fab,
                      const Box& domain)

{
  BL_PROFILE("mfix::set_tracer_bcs()");

  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<const int> const& bct_ilo = bc_list.bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_list.bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_list.bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_list.bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_list.bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_list.bc_khi[lev]->array();

  Array4<Real> const& scal_arr = scal_fab.array();

  Real bc0 = fluid.trac_0;

  IntVect scal_lo(scal_fab.loVect());
  IntVect scal_hi(scal_fab.hiVect());

  const int nlft = amrex::max(0, dom_lo[0]-scal_lo[0]);
  const int nbot = amrex::max(0, dom_lo[1]-scal_lo[1]);
  const int ndwn = amrex::max(0, dom_lo[2]-scal_lo[2]);

  const int nrgt = amrex::max(0, scal_hi[0]-dom_hi[0]);
  const int ntop = amrex::max(0, scal_hi[1]-dom_hi[1]);
  const int nup  = amrex::max(0, scal_hi[2]-dom_hi[2]);

  // Create InVects for following 2D Boxes
  IntVect bx_yz_lo_lo_2D(scal_lo), bx_yz_lo_hi_2D(scal_hi);
  IntVect bx_yz_hi_lo_2D(scal_lo), bx_yz_hi_hi_2D(scal_hi);
  IntVect bx_xz_lo_lo_2D(scal_lo), bx_xz_lo_hi_2D(scal_hi);
  IntVect bx_xz_hi_lo_2D(scal_lo), bx_xz_hi_hi_2D(scal_hi);
  IntVect bx_xy_lo_lo_2D(scal_lo), bx_xy_lo_hi_2D(scal_hi);
  IntVect bx_xy_hi_lo_2D(scal_lo), bx_xy_hi_hi_2D(scal_hi);

  // Fix lo and hi limits
  bx_yz_lo_lo_2D[0] = dom_lo[0]-1;
  bx_yz_lo_hi_2D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_2D[0] = dom_hi[0]+1;
  bx_yz_hi_hi_2D[0] = dom_hi[0]+1;

  bx_xz_lo_lo_2D[1] = dom_lo[1]-1;
  bx_xz_lo_hi_2D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_2D[1] = dom_hi[1]+1;
  bx_xz_hi_hi_2D[1] = dom_hi[1]+1;

  bx_xy_lo_lo_2D[2] = dom_lo[2]-1;
  bx_xy_lo_hi_2D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_2D[2] = dom_hi[2]+1;
  bx_xy_hi_hi_2D[2] = dom_hi[2]+1;

  // Create 2D boxes for GPU loops
  const Box bx_yz_lo_2D(bx_yz_lo_lo_2D, bx_yz_lo_hi_2D);
  const Box bx_yz_hi_2D(bx_yz_hi_lo_2D, bx_yz_hi_hi_2D);

  const Box bx_xz_lo_2D(bx_xz_lo_lo_2D, bx_xz_lo_hi_2D);
  const Box bx_xz_hi_2D(bx_xz_hi_lo_2D, bx_xz_hi_hi_2D);

  const Box bx_xy_lo_2D(bx_xy_lo_lo_2D, bx_xy_lo_hi_2D);
  const Box bx_xy_hi_2D(bx_xy_hi_lo_2D, bx_xy_hi_hi_2D);

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

  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D,
      [bct_ilo,dom_lo,bc0,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      // const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if ((bct == BCList::pinf) || (bct == BCList::pout))
      {
        scal_arr(i,j,k) = scal_arr(dom_lo[0],j,k);
      }
      else if (bct == BCList::minf)
      {
        scal_arr(i,j,k) = bc0;
      }
    });

  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D,
      [bct_ihi,dom_hi,bc0,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      // const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
        scal_arr(i,j,k) = scal_arr(dom_hi[0],j,k);
      }
      else if(bct == BCList::minf)
      {
        scal_arr(i,j,k) = bc0;
      }
    });
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D,
      [bct_jlo,dom_lo,bc0,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      // const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
        scal_arr(i,j,k) = scal_arr(i,dom_lo[1],k);
      }
      else if(bct == BCList::minf)
      {
        scal_arr(i,j,k) = bc0;
      }
    });
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D,
      [bct_jhi,dom_hi,bc0,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      // const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
        scal_arr(i,j,k) = scal_arr(i,dom_hi[1],k);
      }
      else if(bct == BCList::minf)
      {
        scal_arr(i,j,k) = bc0;
      }
    });
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D,
      [bct_klo,dom_lo,bc0,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      // const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
        scal_arr(i,j,k) = scal_arr(i,j,dom_lo[2]);
      }
      else if(bct == BCList::minf)
      {
        scal_arr(i,j,k) = bc0;
      }
    });
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D,
      [bct_khi,dom_hi,bc0,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      // const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if ((bct == BCList::pinf) || (bct == BCList::pout))
      {
        scal_arr(i,j,k) = scal_arr(i,j,dom_hi[2]);
      }
      else if (bct == BCList::minf)
      {
        scal_arr(i,j,k) = bc0;
      }
    });
  }
}
