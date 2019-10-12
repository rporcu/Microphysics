#include <mfix_redist_diff.hpp>
#include <param_mod_F.H>
#include <mfix_util_F.H>

#include <cmath>
#include <limits>

namespace redist_diff_aux {

void
compute_delta_mass(const Box& grown1_bx,
                   MFIter* mfi,
                   MultiFab& divtau_aux,
                   MultiFab& ep_g,
                   FArrayBox& delm_fbx,
                   FArrayBox& optmp_fbx,
                   FArrayBox& mask_fbx,
                   const MultiFab* volfrac,
                   const EBCellFlagFab& flags_fab)
{
  Array4<Real> const& divc = divtau_aux.array(*mfi);
  Array4<Real> const& epsilon_g = ep_g.array(*mfi);

  Array4<Real> const& delm = delm_fbx.array();
  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& mask = mask_fbx.array();

  Array4<const Real> const& vfrac = volfrac->array(*mfi);

  Array4<const EBCellFlag> const& flags = flags_fab.array();

  AMREX_FOR_4D(grown1_bx, 3, i, j, k, n,
  {
    if(flags(i,j,k).isSingleValued())
    {
      Real divnc = 0;
      Real vtot = 0;

      Real epvfrac = 0;

      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            if((ii != 0 or jj != 0 or kk != 0) and 
                (flags(i,j,k).isConnected(ii,jj,kk) == 1))
            {
              epvfrac = vfrac(i+ii,j+jj,k+kk) * epsilon_g(i+ii,j+jj,k+kk) * 
                        mask(i+ii,j+jj,k+kk);
              vtot += epvfrac;
              divnc += epvfrac * divc(i+ii,j+jj,k+kk,n);
            }

      divnc /= vtot;
      epvfrac = vfrac(i,j,k) * epsilon_g(i,j,k);

      Real optmp_value = (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k,n));
      optmp(i,j,k,n) = optmp_value;
      delm(i,j,k,n) = (-1 * epvfrac * optmp_value) * mask(i,j,k);
    }
    else
      delm(i,j,k,n) = 0;
  });

  Gpu::synchronize();
}

void
redistribute_mass(const Box& grown1_bx,
                  MFIter* mfi,
                  MultiFab& ep_g,
                  FArrayBox& delm_fbx,
                  FArrayBox& optmp_fbx,
                  FArrayBox& mask_fbx,
                  const MultiFab* volfrac,
                  const EBCellFlagFab& flags_fab)
{
  Array4<Real> const& epsilon_g = ep_g.array(*mfi);

  Array4<Real> const& delm = delm_fbx.array();
  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& mask = mask_fbx.array();

  Array4<const Real> const& vfrac = volfrac->array(*mfi);

  Array4<const EBCellFlag> const& flags = flags_fab.array();

  AMREX_FOR_3D(grown1_bx, i, j, k,
  {
    if(flags(i,j,k).isSingleValued())
    {
      Real wtot = 0;
      
      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            if((ii != 0 or jj != 0 or kk != 0) and
                (flags(i,j,k).isConnected(ii,jj,kk) == 1))
            {
              wtot += epsilon_g(i+ii,j+jj,k+kk) * vfrac(i+ii,j+jj,k+kk) * 
                      mask(i+ii,j+jj,k+kk);
            }

      wtot = 1/wtot;
      
      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            if((ii != 0 or jj != 0 or kk != 0) and
                (flags(i,j,k).isConnected(ii,jj,kk) == 1))
            {
              for(int n(0); n < 3; ++n)
              {
                
                // Note: delm has already been multiplied by mask
#ifdef AMREX_USE_CUDA
                Gpu::Atomic::Add(&optmp(i+ii,j+jj,k+kk,n), delm(i,j,k,n) * wtot);
#else
                optmp(i+ii,j+jj,k+kk,n) += delm(i,j,k,n) * wtot;
#endif
              }
            }
    }
  });

  Gpu::synchronize();
}

} // end namespace redist_diff_aux

using namespace redist_diff_aux;

void
compute_redist_diff(Box& bx,
                    MultiFab& divtau,
                    MultiFab& ep_g,
                    MultiFab& divtau_aux,
                    MFIter* mfi,
                    const EBCellFlagFab& flags_fab,
                    const MultiFab* volfrac,
                    const int cyclic_x,
                    const int cyclic_y,
                    const int cyclic_z,
                    Box& domain,
                    const Real* dx)
{
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& divergence = divtau.array(*mfi);
  Array4<Real> const& divc = divtau_aux.array(*mfi);

  const Box& grown1_bx = amrex::grow(bx,1);
  const Box& grown2_bx = amrex::grow(bx,2);

  FArrayBox delm_fbx(grown1_bx, 3);
  FArrayBox optmp_fbx(grown2_bx, 3);
  FArrayBox mask_fbx(grown2_bx);

  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& mask = mask_fbx.array();

  //
  // Array "mask" is used to sever the link to ghost cells when the BCs are not
  // periodic
  // It is set to 1 when a cell can be used in computations, 0 otherwise
  //
  AMREX_FOR_3D(grown2_bx, i, j, k,
  {
    if(((not cyclic_x) and (i < dom_low.x or i > dom_high.x)) or
       ((not cyclic_y) and (j < dom_low.y or j > dom_high.y)) or
       ((not cyclic_z) and (k < dom_low.z or k > dom_high.z)))
      mask(i,j,k) = 0;
    else
      mask(i,j,k) = 1;
  });

  //
  // Here we do the redistribution steps
  //
    
  // Set this to zero here
  setFabVal(optmp_fbx, 0.0, grown2_bx, 0, 3);

  //
  // Step 2: compute delta M (mass gain or loss) on (lo-1,lo+1)
  //
  compute_delta_mass(grown1_bx, mfi, divtau_aux, ep_g, delm_fbx, optmp_fbx,
      mask_fbx, volfrac, flags_fab);

  //
  // Step 3: redistribute excess/loss of mass
  //
  redistribute_mass(grown1_bx, mfi, ep_g, delm_fbx, optmp_fbx, mask_fbx,
      volfrac, flags_fab);

  //
  // Resume the correct sign, AKA return the negative
  //
  AMREX_FOR_4D(bx, 3, i, j, k, n,
  {
    divergence(i,j,k,n) = divc(i,j,k,n) + optmp(i,j,k,n);
  });
  
  Gpu::synchronize();
}
