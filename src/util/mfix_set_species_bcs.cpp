#include <mfix.H>

#include <mfix_fluid_parms.H>
#include <mfix_species_parms.H>

using namespace amrex;

//
// Set the BCs for density only
//
void
mfix::mfix_set_species_bcs (Real time,
                            Vector< MultiFab* > const& X_gk_in)
{
  BL_PROFILE("mfix::mfix_set_species_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*(m_leveldata[lev]->ep_g), false); mfi.isValid(); ++mfi)
     {
        set_mass_fractions_g_bcs(time, lev, (*X_gk_in[lev])[mfi], domain);
     }

     X_gk_in[lev]->FillBoundary(geom[lev].periodicity());

     EB_set_covered(*X_gk_in[lev], 0, X_gk_in[lev]->nComp(), X_gk_in[lev]->nGrow(),
         covered_val);
  }
}

void 
mfix::set_mass_fractions_g_bcs (Real time,
                                const int lev,
                                FArrayBox& X_gk_fab,
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

  Array4<Real> const& X_gk = X_gk_fab.array();

  IntVect X_gk_lo(X_gk_fab.loVect());
  IntVect X_gk_hi(X_gk_fab.hiVect());

  const int nlft = std::max(0, dom_lo[0]-X_gk_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-X_gk_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-X_gk_lo[2]);

  const int nrgt = std::max(0, X_gk_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, X_gk_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, X_gk_hi[2]-dom_hi[2]);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(X_gk_hi), bx_xz_lo_hi_3D(X_gk_hi), bx_xy_lo_hi_3D(X_gk_hi);
  IntVect bx_yz_hi_lo_3D(X_gk_lo), bx_xz_hi_lo_3D(X_gk_lo), bx_xy_hi_lo_3D(X_gk_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for GPU loops
  const Box bx_yz_lo_3D(X_gk_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, X_gk_hi);

  const Box bx_xz_lo_3D(X_gk_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, X_gk_hi);

  const Box bx_xy_lo_3D(X_gk_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, X_gk_hi);

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  const int nspecies_g = fluid.nspecies;

  set_species_bc_values(time);
  Real** p_bc_X_gk = m_bc_X_gk_ptr.data();

  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D, nspecies_g, [bct_ilo,dom_lo,pout,minf,pinf,
        X_gk,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);

      if((bct == pout)) {
        X_gk(i,j,k,n) = X_gk(dom_lo[0],j,k,n);
      }
      else if ((bct == minf) || (bct == pinf)) {
        X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
      }
    });
  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D, nspecies_g, [bct_ihi,dom_hi,pout,minf,pinf,
        X_gk,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);

      if(bct == pout) {
        X_gk(i,j,k,n) = X_gk(dom_hi[0],j,k,n);
      }
      else if ((bct == minf) || (bct == pinf)) {
        X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
      }
    });
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D, nspecies_g, [bct_jlo,dom_lo,pout,minf,pinf,
        X_gk,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);

      if(bct == pout) {
        X_gk(i,j,k,n) = X_gk(i,dom_lo[1],k,n);
        }
      else if ((bct == minf) || (bct == pinf))
        X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
    });
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D, nspecies_g, [bct_jhi,dom_hi,pout,minf,pinf,
        X_gk,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);

      if(bct == pout)
        X_gk(i,j,k,n) = X_gk(i,dom_hi[1],k,n);
      else if ((bct == minf) || (bct == pinf))
        X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
    });
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D, nspecies_g, [bct_klo,dom_lo,pout,minf,pinf,
        X_gk,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);

      if(bct == pout)
        X_gk(i,j,k,n) = X_gk(i,j,dom_lo[2],n);
      else if ((bct == minf) || (bct == pinf))
        X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
    });
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D, nspecies_g, [bct_khi,dom_hi,pout,minf,pinf,
        X_gk,p_bc_X_gk]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);

      if(bct == pout)
        X_gk(i,j,k,n) = X_gk(i,j,dom_hi[2],n);
      else if ((bct == minf) || (bct == pinf))
        X_gk(i,j,k,n) = p_bc_X_gk[n][bcv];
    });
  }
}
