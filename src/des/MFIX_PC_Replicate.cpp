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

    Vector<Real> orig_domain_size;
    orig_domain_size.resize(BL_SPACEDIM);
    for (int d = 0; d < BL_SPACEDIM; d++)
        orig_domain_size[d] = (geom.ProbHi(d) - geom.ProbLo(d)) / Nrep[d];

    for (int idim = 0; idim < 3; ++idim)
    {
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& particles = pti.GetArrayOfStructs();
            int np = pti.numParticles();
            Gpu::HostVector<ParticleType> host_particles(np);
            Gpu::copy(Gpu::deviceToHost, particles.begin(), particles.end(), host_particles.begin());

            Gpu::HostVector<ParticleType> replicated_particles;

            for (int i = 0; i < Nrep[idim]; i++)
            {
                if (i == 0) continue; // skip the first copy in each direction

                RealVect shift = {0., 0., 0.};
                shift[idim] = i * orig_domain_size[idim];

                ParticleType p_rep;
                for (const auto& p: host_particles)
                {
                    p_rep.m_rdata.pos[0] = p.m_rdata.pos[0] + shift[0];
                    p_rep.m_rdata.pos[1] = p.m_rdata.pos[1] + shift[1];
                    p_rep.m_rdata.pos[2] = p.m_rdata.pos[2] + shift[2];

                    p_rep.rdata(realData::velx)   = p.rdata(realData::velx);
                    p_rep.rdata(realData::vely)   = p.rdata(realData::vely);
                    p_rep.rdata(realData::velz)   = p.rdata(realData::velz);

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
                    p_rep.id()  = ParticleType::NextID();
                    p_rep.cpu() = ParallelDescriptor::MyProc();

                    // Add everything to the data structure
                    replicated_particles.push_back(p_rep);
                } // p
            } // i

            if (replicated_particles.size() == 0) continue;

            auto new_np = np + replicated_particles.size();
            particles.resize(new_np);
            Gpu::copy(Gpu::hostToDevice, replicated_particles.begin(), replicated_particles.end(),
                      particles.begin() + np);
        } // pti

        Redistribute();
    } // idim
}
