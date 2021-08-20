#include <AMReX.H>
#include <mfix_des_K.H>

#include <mfix_pc.H>
#include <mfix_dem_parms.H>

using namespace amrex;

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

        const Geometry& gm  = Geom(0);
        const auto      p_lo = gm.ProbLoArray();
        const auto      dxi = gm.InvCellSizeArray();

        // This holds the mesh spacing of the level set, which may be finer than
        // the local mesh spacing
        const GpuArray<Real,3> dx_ls{dx[0]/ls_refinement, dx[1]/ls_refinement, dx[2]/ls_refinement};

        const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());

        const Real inv_ep_cp = (PIC::solve) ? 1.0/PIC::ep_cp : 1.0;

        for (MFIXParIter pti(* this, lev); pti.isValid(); ++pti)
        {
            // Real particles

            const Box & bx = pti.tilebox();

            // Remove particles outside of or touching the walls
            if ((*flags)[pti].getType(bx) != FabType::regular)
            {
                auto& aos = pti.GetArrayOfStructs();
                ParticleType* pstruct = aos().dataPtr();

                auto& soa = pti.GetStructOfArrays();
                auto p_realarray = soa.realarray();

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
                    const auto&  phi_arr = ls_phi->array(pti);

                    amrex::ParallelFor(np, [pstruct,p_realarray,plo,dx,flag_fab,inv_ep_cp,
                    dx_ls,phi_arr, ls_refinement, p_lo, dxi, cg_dem=DEM::cg_dem]
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


                          RealVect pos(p.pos());
                          Real ls_value = interp_level_set(pos, ls_refinement, phi_arr, p_lo, dxi);

                          Real radius = p_realarray[SoArealData::radius][ip] *
                            std::cbrt(p_realarray[SoArealData::statwt][ip] * inv_ep_cp);

                          if (cg_dem) {
                            radius = radius/std::cbrt(p_realarray[SoArealData::statwt][ip]);
                          }

                          const Real overlap = radius - ls_value;

                          if (overlap > 0.0) {
                            p.id() = -1;
                          }
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

        ParallelDescriptor::ReduceLongSum(fin_np);
        amrex::Print() << "Final number of particles on level "
                       << lev << ": " << fin_np << std::endl;
        m_total_numparticle = fin_np;
    }
}
