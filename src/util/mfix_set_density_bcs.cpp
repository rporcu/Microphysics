#include <mfix.H>
#include <bc_mod_F.H>
#include <fld_constants_mod_F.H>
#include <param_mod_F.H>

using namespace amrex;

//
// Set the BCs for density only
//
void
mfix::mfix_set_density_bcs (Real time,
                            Vector< std::unique_ptr<MultiFab> > & ro_g_in)
{
  BL_PROFILE("mfix::mfix_set_density_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*ep_g[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
     {
        set_density_bcs(time, lev, (*ro_g_in[lev])[mfi], domain);
     }

     ro_g_in[lev] -> FillBoundary (geom[lev].periodicity());

     EB_set_covered(*ro_g_in[lev], 0, ro_g_in[lev]->nComp(), ro_g_in[lev]->nGrow(), covered_val);
  }
}

void 
mfix::set_density_bcs(Real time,
                      const int lev,
                      FArrayBox& scal_fab,
                      const Box& domain)

{
  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  Array4<Real> const& scal_arr = scal_fab.array();

  Real bc0 = get_ro_g0();

  IntVect scal_lo(scal_fab.loVect());
  IntVect scal_hi(scal_fab.hiVect());

  const int nlft = std::max(0, dom_lo[0]-scal_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-scal_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-scal_lo[2]);

  const int nrgt = std::max(0, scal_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, scal_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, scal_hi[2]-dom_hi[2]);

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

  // Create 2D boxes for CUDA loops
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

  // Create 3D boxes for CUDA loops
  const Box bx_yz_lo_3D(scal_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, scal_hi);

  const Box bx_xz_lo_3D(scal_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, scal_hi);

  const Box bx_xy_lo_3D(scal_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, scal_hi);

  const Real undefined = get_undefined();

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  amrex::Real* p_bc_ep_g = m_bc_ep_g.data();
  amrex::Real* p_bc_t_g  = m_bc_t_g.data();

  if (nlft > 0)
  {
    int ilo = dom_lo[0];
    amrex::ParallelFor(bx_yz_lo_3D,
      [bct_ilo,ilo,bc0,pinf,pout,minf,undefined,p_bc_t_g,scal_arr] 
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ilo(ilo-1,j,k,0);

      if ((bct == pinf) or (bct == pout))
         scal_arr(i,j,k) = scal_arr(ilo,j,k);
      else if (bct == minf)
         scal_arr(i,j,k) = bc0;
    });
  }

  if (nrgt > 0)
  {
    int ihi = dom_hi[0];
    amrex::ParallelFor(bx_yz_hi_3D,
      [bct_ihi,ihi,bc0,pinf,pout,minf,undefined,p_bc_t_g,scal_arr] 
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ihi(ihi+1,j,k,0);

      if((bct == pinf) or (bct == pout))
         scal_arr(i,j,k) = scal_arr(ihi,j,k);
      else if(bct == minf)
         scal_arr(i,j,k) = bc0;
    });
  }

  if (nbot > 0)
  {
    int jlo = dom_lo[1];
    amrex::ParallelFor(bx_xz_lo_3D,
      [bct_jlo,jlo,bc0,pinf,pout,minf,undefined,p_bc_t_g,scal_arr] 
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jlo(i,jlo-1,k,0);

      if((bct == pinf) or (bct == pout))
         scal_arr(i,j,k) = scal_arr(i,jlo,k);
      else if(bct == minf)
         scal_arr(i,j,k) = bc0;
    });
  }

  if (ntop > 0)
  {
    int jhi = dom_hi[1];
    amrex::ParallelFor(bx_xz_hi_3D,
      [bct_jhi,jhi,dom_hi,bc0,pinf,pout,minf,undefined,p_bc_t_g,scal_arr] 
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_jhi(i,jhi+1,k,0);
  
        if((bct == pinf) or (bct == pout))
           scal_arr(i,j,k) = scal_arr(i,dom_hi[1],k);
        else if(bct == minf)
           scal_arr(i,j,k) = bc0;
      });
  }

  if (ndwn > 0)
  {
    int klo = dom_lo[2];
    amrex::ParallelFor(bx_xy_lo_3D,
      [bct_klo,klo,dom_lo,bc0,pinf,pout,minf,undefined,p_bc_t_g,scal_arr] 
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_klo(i,j,klo-1,0);

        if((bct == pinf) or (bct == pout))
           scal_arr(i,j,k) = scal_arr(i,j,klo);
        else if(bct == minf)
           scal_arr(i,j,k) = bc0;
      });
  }

  if (nup > 0)
  {
    int khi = dom_hi[2];
    amrex::ParallelFor(bx_xy_hi_3D, 
      [bct_khi,khi,dom_hi,bc0,pinf,pout,minf,undefined,p_bc_t_g,scal_arr] 
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bct = bct_khi(i,j,khi+1,0);
  
        if ((bct == pinf) or (bct == pout))
           scal_arr(i,j,k) = scal_arr(i,j,dom_hi[2]);
        else if (bct == minf)
           scal_arr(i,j,k) = bc0;
      });
  }
}
