#include <AMReX.H>
#include "AMReX_Particles.H"
#include "AMReX_RealVect.H"
#include <iostream>
#include <MFIXParticleContainer.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB_F.H>

#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBSupport.H>
#include <AMReX_EBMultiFabUtil.H>

#include <math.h>

#include "mfix_F.H"
#include "mfix_des_F.H"
#include "mfix_eb_F.H"
#include "mfix_util_F.H"

using namespace amrex;
using namespace std;

void MFIXParticleContainer::RemoveOutOfRange(int lev, const EBFArrayBoxFactory * ebfactory,
                                             const MultiFab * ls_phi, int ls_refinement)
{
    // Only call the routine for wall collisions if we actually have walls
    if (ebfactory != NULL) {

        const Real * dx  = Geom(lev).CellSize();
        const Real * plo = Geom(lev).ProbLo();

        // This holds the mesh spacing of the level set, which may be finer than the local mesh spacing
        Real dx_ls[3];
        for (int i = 0; i < 3; i++)
           dx_ls[i] = dx[i] / ls_refinement;

        const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());

        for (MFIXParIter pti(* this, lev); pti.isValid(); ++pti) 
        {
            // Real particles
            const int nrp = NumberOfParticles(pti);

            void * particles  = pti.GetArrayOfStructs().data();

            const Box & bx = pti.tilebox();

            // Remove particles outside of or touching the walls
            if ((*flags)[pti].getType(bx) != FabType::regular)
            {
                if ((*flags)[pti].getType(bx) == FabType::covered)
                {
                    for (int ip = 0; ip < pti.numParticles(); ++ip)
                    {
                        ParticleType& p = pti.GetArrayOfStructs()[ip];
                        p.id() = -1;
                    }
                }
                else
                {
                    const auto& flag_fab =  flags->array(pti);
                    const auto&  phi_fab = ls_phi->array(pti);
                    for (int ip = 0; ip < pti.numParticles(); ++ip)
                    {
                        ParticleType& p = pti.GetArrayOfStructs()[ip];

                        int ic = floor( ( p.pos(0) - plo[0] ) / dx[0] );
                        int jc = floor( ( p.pos(1) - plo[1] ) / dx[1] );
                        int kc = floor( ( p.pos(2) - plo[2] ) / dx[2] );
                        
                        if (flag_fab(ic,jc,kc).isCovered())
                        {
                            p.id() = -1;
                        } 
                        else 
                        { // Interpolates level-set from nodal phi to position pos
                            
                            Real x = ( p.pos(0) - plo[0] ) / dx_ls[0];
                            Real y = ( p.pos(1) - plo[1] ) / dx_ls[1];
                            Real z = ( p.pos(2) - plo[2] ) / dx_ls[2];
                            
                            int i = floor(x);
                            int j = floor(y);
                            int k = floor(z);
                            
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

                            if (phi_interp < p.rdata(realData::radius)) p.id() = -1;
                        }
                    }
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
