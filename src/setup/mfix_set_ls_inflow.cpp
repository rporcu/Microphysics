#include <mfix.H>

namespace set_ls_inflow_aux {

AMREX_GPU_HOST_DEVICE
bool
is_equal_to_any(const int bc,
                const int* bc_types,
                const int size)
{
  for(int i(0); i < size; ++i)
  {
    if(bc == bc_types[i])
      return true;
  }
  return false;
}

} // end namespace set_ls_inflow_aux

using namespace set_ls_inflow_aux;

void 
mfix::set_ls_inflow(const int lev,
                    FArrayBox& ls_phi_fab,
                    const Box& domain,
                    const int* ng,
                    const int& nref,
                    const Real* dx)
{
  const Real offset(1.e-8);

  Array4<double> const& ls_phi = ls_phi_fab.array();

  const IntVect sbx_lo(ls_phi_fab.loVect());
  const IntVect sbx_hi(ls_phi_fab.hiVect());
  const Box sbx(sbx_lo, sbx_hi);

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();
  
  // Here if the level set (slo,shi) is at a finer resolution (by nref) than
  //  the boundary condition routines,
  //  the domain boundaries (domlo,domhi), and dx,
  //     then we make sure to adjust by nref
  const amrex::GpuArray<const Real, AMREX_SPACEDIM> dx_fine = 
    {dx[0]/Real(nref), dx[1]/Real(nref), dx[2]/Real(nref)};

  const int nlft = amrex::max(0, nref*dom_lo[0]-sbx_lo[0]);
  const int nbot = amrex::max(0, nref*dom_lo[1]-sbx_lo[1]);
  const int ndwn = amrex::max(0, nref*dom_lo[2]-sbx_lo[2]);

  const int nrgt = amrex::max(0, sbx_hi[0]-(nref*dom_hi[0]+1));
  const int ntop = amrex::max(0, sbx_hi[1]-(nref*dom_hi[1]+1));
  const int nup  = amrex::max(0, sbx_hi[2]-(nref*dom_hi[2]+1));

  const int minf = bc_list.get_minf();

  if (nlft > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(sbx, i, j, k,
    {
      int bct[4];
      bct[0] = bct_ilo(dom_lo[0]-1,j/nref,k/nref,0);
      bct[1] = bct_ilo(dom_lo[0]-1,j/nref,k/nref+1,0);
      bct[2] = bct_ilo(dom_lo[0]-1,j/nref+1,k/nref,0);
      bct[3] = bct_ilo(dom_lo[0]-1,j/nref+1,k/nref+1,0);

      if(is_equal_to_any(minf, &bct[0], 4))
      {
        if(i < 0)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = Real(i)*(dx_fine[0]) - offset;
          else
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(i)*(dx_fine[0]) - offset);
        }
        if(i == 0)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = - offset;
        }
        if(i > 0)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(i)*(dx_fine[0]) - offset);
        }
      }
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (nrgt > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(sbx, i, j, k,
    {
      int bct[4];
      bct[0] = bct_ihi(dom_hi[0]+1,j/nref,k/nref,0);
      bct[1] = bct_ihi(dom_hi[0]+1,j/nref,k/nref+1,0);
      bct[2] = bct_ihi(dom_hi[0]+1,j/nref+1,k/nref,0);
      bct[3] = bct_ihi(dom_hi[0]+1,j/nref+1,k/nref+1,0);

      if(is_equal_to_any(minf, &bct[0], 4))
      {
        if(i < (dom_hi[0]+1)*nref)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(dom_hi[0]-i+1)*(dx_fine[0]) - offset);
        }
        if(i == (dom_hi[0]+1)*nref)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = - offset;
        }
        if(i > (dom_hi[0]+1)*nref)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = Real(dom_hi[0]-i+1)*(dx_fine[0]) - offset;
          else
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(i)*(dx_fine[0]) - offset);
        }
      }
    });
  }
  
#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (nbot > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(sbx, i, j, k,
    {
      int bct[4];
      bct[0] = bct_jlo(i/nref,dom_lo[1]-1,k/nref,0);
      bct[1] = bct_jlo(i/nref,dom_lo[1]-1,k/nref+1,0);
      bct[2] = bct_jlo(i/nref+1,dom_lo[1]-1,k/nref,0);
      bct[3] = bct_jlo(i/nref+1,dom_lo[1]-1,k/nref+1,0);

      if(is_equal_to_any(minf, &bct[0], 4))
      {
        if(j < 0)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = Real(j)*(dx_fine[1]) - offset;
          else
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(j)*(dx_fine[1]) - offset);
        }
        if(j == 0)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = - offset;
        }
        if(j > 0)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(j)*(dx_fine[1]) - offset);
        }
      }
    });
  }
  
#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (ntop > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(sbx, i, j, k,
    {
      int bct[4];
      bct[0] = bct_jhi(i/nref,dom_hi[1]+1,k/nref,0);
      bct[1] = bct_jhi(i/nref,dom_hi[1]+1,k/nref+1,0);
      bct[2] = bct_jhi(i/nref+1,dom_hi[1]+1,k/nref,0);
      bct[3] = bct_jhi(i/nref+1,dom_hi[1]+1,k/nref+1,0);

      if(is_equal_to_any(minf, &bct[0], 4))
      {
        if(j < (dom_hi[1]+1)*nref)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(dom_hi[1]-j+1)*(dx_fine[1]) - offset);
        }
        if(j == (dom_hi[1]+1)*nref)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = - offset;
        }
        if(j > (dom_hi[1]+1)*nref)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = Real(dom_hi[1]-j+1)*(dx_fine[1]) - offset;
          else
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(j)*(dx_fine[1]) - offset);
        }
      }
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (ndwn > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(sbx, i, j, k,
    {
      int bct[4];
      bct[0] = bct_klo(i/nref,j/nref,dom_lo[2]-1,0);
      bct[1] = bct_klo(i/nref,j/nref+1,dom_lo[2]-1,0);
      bct[2] = bct_klo(i/nref+1,j/nref,dom_lo[2]-1,0);
      bct[3] = bct_klo(i/nref+1,j/nref+1,dom_lo[2]-1,0);

      if(is_equal_to_any(minf, &bct[0], 4))
      {
        if(k < 0)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = Real(k)*(dx_fine[2]) - offset;
          else
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(k)*(dx_fine[2]) - offset);
        }
        if(k == 0)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = - offset;
        }
        if(k > 0)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(k)*(dx_fine[2]) - offset);
        }
      }
    });
  }
  
#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  if (nup > 0)
  {
    AMREX_HOST_DEVICE_FOR_3D(sbx, i, j, k,
    {
      int bct[4];
      bct[0] = bct_khi(i/nref,j/nref,dom_hi[2]+1,0);
      bct[1] = bct_khi(i/nref,j/nref+1,dom_hi[2]+1,0);
      bct[2] = bct_khi(i/nref+1,j/nref,dom_hi[2]+1,0);
      bct[3] = bct_khi(i/nref+1,j/nref+1,dom_hi[2]+1,0);

      if(is_equal_to_any(minf, &bct[0], 4))
      {
        if(k < (dom_hi[2]+1)*nref)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(dom_hi[2]-k+1)*(dx_fine[2]) - offset);
        }
        if(k == (dom_hi[2]+1)*nref)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = - offset;
        }
        if(k > (dom_hi[2]+1)*nref)
        {
          if(ls_phi(i,j,k) > 0)
            ls_phi(i,j,k) = Real(dom_hi[2]-k+1)*(dx_fine[2]) - offset;
          else
            ls_phi(i,j,k) = amrex::min(ls_phi(i,j,k), Real(k)*(dx_fine[2]) - offset);
        }
      }
    });
  }

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif
}
