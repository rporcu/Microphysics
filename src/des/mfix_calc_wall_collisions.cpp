#include <mfix.H>

using namespace amrex;

void MFIXParticleContainer::ls_has_walls(int& has_wall,
                                         const Box& bx,
                                         const FArrayBox& phi_fab,
                                         const Real tolerance)
{
  has_wall = 0;

  Array4<const Real> const& phi = phi_fab.array();

#ifdef AMREX_USE_CUDA
  Gpu::DeviceScalar<int> has_wall_gpu(has_wall);
  int* p_has_wall = has_wall_gpu.dataPtr();
#endif

  amrex::ParallelFor(bx, [phi,tolerance,
#ifdef AMREX_USE_CUDA
      p_has_wall]
#else
      &has_wall]
#endif
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    if(phi(i,j,k) <= tolerance)
#ifdef AMREX_USE_CUDA
      *p_has_wall = 1;
#else
      has_wall = 1;
#endif
  });

#ifdef AMREX_USE_CUDA
  Gpu::synchronize();
  has_wall = has_wall_gpu.dataValue();
#endif
}
