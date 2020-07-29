#include <AMReX.H>
#include "AMReX_Particles.H"
#include "AMReX_RealVect.H"
#include <iostream>
#include <mfix_pc.H>
#include <mfix_dem_parms.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB_F.H>

#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBSupport.H>
#include <AMReX_EBMultiFabUtil.H>

#include <math.h>

using namespace amrex;
using namespace std;

void MFIXParticleContainer::RemoveOutOfRange (int lev,
                                              const EBFArrayBoxFactory * ebfactory,
                                              const MultiFab * ls_phi,
                                              int ls_refinement)
{
    // Only call the routine for wall collisions if we actually have walls
    if (ebfactory != NULL) {

        const Real* cell_size  = Geom(lev).CellSize();

        const RealVect dx(cell_size[0], cell_size[1], cell_size[2]);
        const GpuArray<Real,3> plo = Geom(lev).ProbLoArray();

        // This holds the mesh spacing of the level set, which may be finer than
        // the local mesh spacing
        const GpuArray<Real,3> dx_ls{dx[0]/ls_refinement, dx[1]/ls_refinement, dx[2]/ls_refinement};

        const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());

        for (MFIXParIter pti(* this, lev); pti.isValid(); ++pti)
        {
            // Real particles

            const Box & bx = pti.tilebox();

            // Remove particles outside of or touching the walls
            if ((*flags)[pti].getType(bx) != FabType::regular)
            {
                auto& aos = pti.GetArrayOfStructs();
                ParticleType* pstruct = aos().dataPtr();
                const int np = pti.numParticles();

                if ((*flags)[pti].getType(bx) == FabType::covered)
                {
                    amrex::ParallelFor(np, [pstruct] AMREX_GPU_DEVICE (int ip) noexcept
                    {
                      ParticleType& p = pstruct[ip];
                      p.id() = -1;
                    });
                }
                else
                {
                    const auto& flag_fab =  flags->array(pti);
                    const auto&  phi_fab = ls_phi->array(pti);

                    amrex::ParallelFor(np, [pstruct,plo,dx,flag_fab,dx_ls,phi_fab,cg_dem=DEM::cg_dem]
                      AMREX_GPU_DEVICE (int ip) noexcept
                    {
                        ParticleType& p = pstruct[ip];

                        const Real* plo_ptr = plo.data();
                        int icell = static_cast<int>(amrex::Math::floor( ( p.pos(0) - plo_ptr[0] ) / dx[0] ));
                        int jcell = static_cast<int>(amrex::Math::floor( ( p.pos(1) - plo_ptr[1] ) / dx[1] ));
                        int kcell = static_cast<int>(amrex::Math::floor( ( p.pos(2) - plo_ptr[2] ) / dx[2] ));

                        if (flag_fab(icell,jcell,kcell).isCovered())
                        {
                            p.id() = -1;
                        }
                        else
                        { // Interpolates level-set from nodal phi to position pos

                            Real x = ( p.pos(0) - plo_ptr[0] ) / dx_ls[0];
                            Real y = ( p.pos(1) - plo_ptr[1] ) / dx_ls[1];
                            Real z = ( p.pos(2) - plo_ptr[2] ) / dx_ls[2];

                            int i = static_cast<int>(amrex::Math::floor(x));
                            int j = static_cast<int>(amrex::Math::floor(y));
                            int k = static_cast<int>(amrex::Math::floor(z));

                            Real wx_hi = x - i;
                            Real wy_hi = y - j;
                            Real wz_hi = z - k;

                            Real wx_lo = 1.0 - wx_hi;
                            Real wy_lo = 1.0 - wy_hi;
                            Real wz_lo = 1.0 - wz_hi;

                            Real phi_interp = phi_fab(i,   j,   k  ) * wx_lo * wy_lo * wz_lo
                                            + phi_fab(i+1, j,   k  ) * wx_hi * wy_lo * wz_lo
                                            + phi_fab(i,   j+1, k  ) * wx_lo * wy_hi * wz_lo
                                            + phi_fab(i,   j,   k+1) * wx_lo * wy_lo * wz_hi
                                            + phi_fab(i+1, j+1, k  ) * wx_hi * wy_hi * wz_lo
                                            + phi_fab(i,   j+1, k+1) * wx_lo * wy_hi * wz_hi
                                            + phi_fab(i+1, j,   k+1) * wx_hi * wy_lo * wz_hi
                                            + phi_fab(i+1, j+1, k+1) * wx_hi * wy_hi * wz_hi;

                            amrex::Real radius = p.rdata(realData::radius) *
                                std::cbrt(p.rdata(realData::statwt));

                            if (cg_dem)
                            {
                               radius = radius/std::cbrt(p.rdata(realData::statwt));
                            }

                            if (phi_interp < radius)
                            {
                                 p.id() = -1;
                            }
#if 0
                            else {
                                 std::cout << " 1 "
                                           << p.pos(0) << " "
                                           << p.pos(1) << " "
                                           << p.pos(2) << " "
                                           << p.rdata(realData::radius)  << " "
                                           << p.rdata(realData::density) << " "
                                           << p.rdata(realData::velx)    << " "
                                           << p.rdata(realData::vely)    << " "
                                           << p.rdata(realData::velz) << std::endl;
                            }
#endif
                        }
                    });
                }
            }
        }

        Redistribute();

        long fin_np = 0;
        for (MFIXParIter pti(* this, lev); pti.isValid(); ++pti) {
            long np = pti.numParticles();
            fin_np += np;
        }

        ParallelDescriptor::ReduceLongSum(fin_np,ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "Final number of particles on level "
                       << lev << ": " << fin_np << std::endl;
    }
}
