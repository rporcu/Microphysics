#include <mfix.H>

#include <MFIX_FLUID_Parms.H>
#include <MFIX_SPECIES_Parms.H>

using namespace amrex;

//
// Set the BCs for density only
//
void
mfix::mfix_set_species_bcs (Real time,
                            Vector< MultiFab* > const& X_g_in,
                            Vector< MultiFab* > const& D_g_in)
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
        set_mass_fractions_g_bcs(time, lev, (*X_g_in[lev])[mfi], domain);
        set_species_diffusivities_g_bcs(time, lev, (*D_g_in[lev])[mfi], domain);
     }

     X_g_in[lev]->FillBoundary(geom[lev].periodicity());
     D_g_in[lev]->FillBoundary(geom[lev].periodicity());

     EB_set_covered(*X_g_in[lev], 0, X_g_in[lev]->nComp(), X_g_in[lev]->nGrow(),
         covered_val);
     EB_set_covered(*D_g_in[lev], 0, D_g_in[lev]->nComp(), D_g_in[lev]->nGrow(),
         covered_val);
  }
}

void 
mfix::set_mass_fractions_g_bcs (Real time,
                                const int lev,
                                FArrayBox& X_g_fab,
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

  Array4<Real> const& X_g = X_g_fab.array();

  IntVect X_g_lo(X_g_fab.loVect());
  IntVect X_g_hi(X_g_fab.hiVect());

  const int nlft = std::max(0, dom_lo[0]-X_g_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-X_g_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-X_g_lo[2]);

  const int nrgt = std::max(0, X_g_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, X_g_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, X_g_hi[2]-dom_hi[2]);

  // Create InVects for following 3D Boxes
  IntVect bx_yz_lo_hi_3D(X_g_hi), bx_xz_lo_hi_3D(X_g_hi), bx_xy_lo_hi_3D(X_g_hi);
  IntVect bx_yz_hi_lo_3D(X_g_lo), bx_xz_hi_lo_3D(X_g_lo), bx_xy_hi_lo_3D(X_g_lo);

  // Fix lo and hi limits
  bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
  bx_yz_hi_lo_3D[0] = dom_hi[0]+1;

  bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
  bx_xz_hi_lo_3D[1] = dom_hi[1]+1;

  bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
  bx_xy_hi_lo_3D[2] = dom_hi[2]+1;

  // Create 3D boxes for GPU loops
  const Box bx_yz_lo_3D(X_g_lo, bx_yz_lo_hi_3D);
  const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, X_g_hi);

  const Box bx_xz_lo_3D(X_g_lo, bx_xz_lo_hi_3D);
  const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, X_g_hi);

  const Box bx_xy_lo_3D(X_g_lo, bx_xy_lo_hi_3D);
  const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, X_g_hi);

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  const int nspecies_g = FLUID::nspecies_g;

  Gpu::ManagedVector< Real* > m_bc_X_g_managed(nspecies_g);

  for (int n(0); n < nspecies_g; n++)
    m_bc_X_g_managed[n] = m_bc_X_g[n].data();

  Real** p_bc_X_g = m_bc_X_g_managed.data();

  if (nlft > 0)
  {
    amrex::ParallelFor(bx_yz_lo_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
      const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);

      if((bct == pout)) {
        X_g(i,j,k,n) = X_g(dom_lo[0],j,k,n);
      }
      else if ((bct == minf) or (bct == pinf)) {
        X_g(i,j,k,n) = p_bc_X_g[n][bcv];
      }
    });
  }

  if (nrgt > 0)
  {
    amrex::ParallelFor(bx_yz_hi_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
      const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);

      if(bct == pout) {
        X_g(i,j,k,n) = X_g(dom_hi[0],j,k,n);
      }
      else if ((bct == minf) or (bct == pinf)) {
        X_g(i,j,k,n) = p_bc_X_g[n][bcv];
      }
    });
  }

  if (nbot > 0)
  {
    amrex::ParallelFor(bx_xz_lo_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
      const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);

      if(bct == pout) {
        X_g(i,j,k,n) = X_g(i,dom_lo[1],k,n);
        }
      else if ((bct == minf) or (bct == pinf))
        X_g(i,j,k,n) = p_bc_X_g[n][bcv];
    });
  }

  if (ntop > 0)
  {
    amrex::ParallelFor(bx_xz_hi_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
      const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);

      if(bct == pout)
        X_g(i,j,k,n) = X_g(i,dom_hi[1],k,n);
      else if ((bct == minf) or (bct == pinf))
        X_g(i,j,k,n) = p_bc_X_g[n][bcv];
    });
  }

  if (ndwn > 0)
  {
    amrex::ParallelFor(bx_xy_lo_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);
      const int bcv = bct_klo(i,j,dom_lo[2]-1,1);

      if(bct == pout)
        X_g(i,j,k,n) = X_g(i,j,dom_lo[2],n);
      else if ((bct == minf) or (bct == pinf))
        X_g(i,j,k,n) = p_bc_X_g[n][bcv];
    });
  }

  if (nup > 0)
  {
    amrex::ParallelFor(bx_xy_hi_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);
      const int bcv = bct_khi(i,j,dom_hi[2]+1,1);

      if(bct == pout)
        X_g(i,j,k,n) = X_g(i,j,dom_hi[2],n);
      else if ((bct == minf) or (bct == pinf))
        X_g(i,j,k,n) = p_bc_X_g[n][bcv];
    });
  }

}


void 
mfix::set_species_diffusivities_g_bcs (Real time,
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

  const int nspecies_g = FLUID::nspecies_g;

  Gpu::ManagedVector< Real > D_g0_managed(nspecies_g);

  for (int n(0); n < nspecies_g; n++)
    D_g0_managed[n] = FLUID::D_g0[n];

  Real* p_D_g0 = D_g0_managed.data();

  IntVect scal_lo(scal_fab.loVect());
  IntVect scal_hi(scal_fab.hiVect());

  const int nlft = std::max(0, dom_lo[0]-scal_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-scal_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-scal_lo[2]);

  const int nrgt = std::max(0, scal_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, scal_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, scal_hi[2]-dom_hi[2]);

  const int minf = bc_list.get_minf();
  const int pinf = bc_list.get_pinf();
  const int pout = bc_list.get_pout();

  if (nlft > 0)
  {
    IntVect bx_yz_lo_hi_3D(scal_hi);
    bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
    const Box bx_yz_lo_3D(scal_lo, bx_yz_lo_hi_3D);

    int ilo = dom_lo[0];

    amrex::ParallelFor(bx_yz_lo_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_ilo(ilo-1,j,k,0);

      if (bct == pout)
         scal_arr(i,j,k,n) = scal_arr(ilo,j,k,n);
      else if (bct == minf or bct == pinf)
         scal_arr(i,j,k,n) = p_D_g0[n];
    });
  }

  if (nrgt > 0)
  {
    IntVect bx_yz_hi_lo_3D(scal_lo);
    bx_yz_hi_lo_3D[0] = dom_hi[0]+1;
    const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, scal_hi);

    int ihi = dom_hi[0];

    amrex::ParallelFor(bx_yz_hi_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_ihi(ihi+1,j,k,0);

      if (bct == pout)
         scal_arr(i,j,k,n) = scal_arr(ihi,j,k,n);
      else if (bct == minf or bct == pinf)
         scal_arr(i,j,k,n) = p_D_g0[n];
    });
  }

  if (nbot > 0)
  {
    IntVect bx_xz_lo_hi_3D(scal_hi);
    bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
    const Box bx_xz_lo_3D(scal_lo, bx_xz_lo_hi_3D);

    int jlo = dom_lo[1];

    amrex::ParallelFor(bx_xz_lo_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_jlo(i,jlo-1,k,0);

      if (bct == pout)
         scal_arr(i,j,k,n) = scal_arr(i,jlo,k,n);
      else if (bct == minf or bct == pinf)
         scal_arr(i,j,k,n) = p_D_g0[n];
    });
  }

  if (ntop > 0)
  {
    IntVect bx_xz_hi_lo_3D(scal_lo);
    bx_xz_hi_lo_3D[1] = dom_hi[1]+1;
    const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, scal_hi);

    int jhi = dom_hi[1];

    amrex::ParallelFor(bx_xz_hi_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_jhi(i,jhi+1,k,0);

      if (bct == pout)
         scal_arr(i,j,k,n) = scal_arr(i,dom_hi[1],k,n);
      else if (bct == minf or bct == pinf)
         scal_arr(i,j,k,n) = p_D_g0[n];
    });
  }

  if (ndwn > 0)
  {
    IntVect bx_xy_lo_hi_3D(scal_hi);
    bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
    const Box bx_xy_lo_3D(scal_lo, bx_xy_lo_hi_3D);

    int klo = dom_lo[2];

    amrex::ParallelFor(bx_xy_lo_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_klo(i,j,klo-1,0);

      if (bct == pout)
         scal_arr(i,j,k,n) = scal_arr(i,j,klo,n);
      else if (bct == minf or bct == pinf)
         scal_arr(i,j,k,n) = p_D_g0[n];
    });
  }

  if (nup > 0)
  {
    IntVect bx_xy_hi_lo_3D(scal_lo);
    bx_xy_hi_lo_3D[2] = dom_hi[2]+1;
    const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, scal_hi);

    int khi = dom_hi[2];

    amrex::ParallelFor(bx_xy_hi_3D, nspecies_g,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_khi(i,j,khi+1,0);

      if (bct == pout)
         scal_arr(i,j,k,n) = scal_arr(i,j,dom_hi[2],n);
      else if (bct == minf or bct == pinf)
         scal_arr(i,j,k,n) = p_D_g0[n];
    });
  }
}
