#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_RealVect.H>
#include <MFIXParticleContainer.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB_F.H>

#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBSupport.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

#include <cmath>
#include <iostream>

#include "mfix_F.H"
#include "mfix_des_K.H"
#include "MFIX_DEM_Parms.H"
#include <MFIX_PIC_Parms.H>

using namespace amrex;
using namespace std;

void MFIXParticleContainer::Replicate (IntVect& Nrep,
                                       Geometry& geom,
                                       DistributionMapping& dmap,
                                       BoxArray& ba)
{
    int lev = 0;

    RealVect orig_domain_size;
    for (int d = 0; d < BL_SPACEDIM; d++)
      orig_domain_size[d] = (geom.ProbHi(d) - geom.ProbLo(d)) / Nrep[d];

    const int myProc = ParallelDescriptor::MyProc();

    for (int idim = 0; idim < 3; ++idim)
    {
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& particles = pti.GetArrayOfStructs();
            int np = pti.numParticles();

            int np_replicated = np * (Nrep[idim] - 1);

            int new_np = np + np_replicated;
            particles.resize(new_np);

            ParticleType* pstruct = particles().dataPtr();

            for (int i = 0; i < Nrep[idim]; i++)
            {
                if (i == 0) continue; // skip the first copy in each direction

                RealVect shift = {0., 0., 0.};
                shift[idim] = i * orig_domain_size[idim];

                const int nextID = ParticleType::NextID();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int n) noexcept
                {
                    int index = n;
                    int index_repl = i*np + n;

                    ParticleType p = pstruct[index];
                    ParticleType& p_rep = pstruct[index_repl];

                    p_rep.pos(0) = p.pos(0) + shift[0];
                    p_rep.pos(1) = p.pos(1) + shift[1];
                    p_rep.pos(2) = p.pos(2) + shift[2];
                
                    p_rep.rdata(realData::velx) = p.rdata(realData::velx);
                    p_rep.rdata(realData::vely) = p.rdata(realData::vely);
                    p_rep.rdata(realData::velz) = p.rdata(realData::velz);
                
                    // Set other particle properties
                    p_rep.idata(intData::phase)       = p.idata(intData::phase);
                    p_rep.idata(intData::state)       = p.idata(intData::state);
                    p_rep.rdata(realData::volume)     = p.rdata(realData::volume);
                    p_rep.rdata(realData::density)    = p.rdata(realData::density);
                    p_rep.rdata(realData::mass)       = p.rdata(realData::mass);
                    p_rep.rdata(realData::oneOverI)   = p.rdata(realData::oneOverI);
                    p_rep.rdata(realData::radius)     = p.rdata(realData::radius);
                    p_rep.rdata(realData::omegax)     = p.rdata(realData::omegax);
                    p_rep.rdata(realData::omegay)     = p.rdata(realData::omegay);
                    p_rep.rdata(realData::omegaz)     = p.rdata(realData::omegaz);
                    p_rep.rdata(realData::statwt)     = p.rdata(realData::statwt);
                    p_rep.rdata(realData::dragcoeff)  = p.rdata(realData::dragcoeff);
                    p_rep.rdata(realData::dragx)      = p.rdata(realData::dragx);
                    p_rep.rdata(realData::dragy)      = p.rdata(realData::dragy);
                    p_rep.rdata(realData::dragz)      = p.rdata(realData::dragz);
                
                    // Set id and cpu for this particle
                    p_rep.id()  = nextID + n;
                    p_rep.cpu() = myProc;
                 
                }); // p

                ParticleType::NextID(nextID + np);
            } // i

            amrex::Gpu::synchronize();

        } // pti

        Redistribute();
    } // idim

}
