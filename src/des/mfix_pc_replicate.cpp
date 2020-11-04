#include <mfix_des_K.H>
#include <mfix_solids_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_pic_parms.H>

using namespace amrex;

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

    const int nspecies_s = SOLIDS::nspecies;
    const int nreactions = REACTIONS::nreactions;

    const int idx_X = speciesData::X_sn * nspecies_s;

    const int idx_G = speciesData::count * nspecies_s +
                      reactionsData::G_sn_pg_q * nreactions;

    for (int idim = 0; idim < 3; ++idim)
    {
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& ptile = DefineAndReturnParticleTile(lev,pti);
            int np = pti.numParticles();

            int np_replicated = np * (Nrep[idim] - 1);

            int new_np = np + np_replicated;
            ptile.resize(new_np);

            // Add runtime-added components for species mass fractions
            if (SOLIDS::solve_species) {
              for (int n_s(0); n_s < SOLIDS::nspecies; ++n_s)
                ptile.push_back_real(n_s, np_replicated, 0.);
            }

            // Add runtime-added components for species mass rate transfer for
            // each reaction
            if (SOLIDS::solve_species and REACTIONS::solve) {
              const int gap = SOLIDS::nspecies;

              for (int n_s(0); n_s < SOLIDS::nspecies; ++n_s)
                for (int q(gap); q < (gap+REACTIONS::nreactions); ++q) {
                  const int comp = q + n_s * REACTIONS::nreactions;
                  ptile.push_back_real(comp, np_replicated, 0.);
                }
            }

            auto& particles = ptile.GetArrayOfStructs();
            ParticleType* pstruct = particles().dataPtr();

            //Access to added variables
            auto ptile_data = ptile.getParticleTileData();

            for (int i = 0; i < Nrep[idim]; i++)
            {
                if (i == 0) continue; // skip the first copy in each direction

                RealVect shift = {0., 0., 0.};
                shift[idim] = i * orig_domain_size[idim];

                const int nextID = ParticleType::NextID();

                amrex::ParallelFor(np, [pstruct,ptile_data,np,nextID,myProc,
                    nspecies_s,nreactions,idx_X,idx_G,shift,i]
                  AMREX_GPU_DEVICE (int n) noexcept
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
                    p_rep.idata(intData::phase)        = p.idata(intData::phase);
                    p_rep.idata(intData::state)        = p.idata(intData::state);
                    p_rep.rdata(realData::volume)      = p.rdata(realData::volume);
                    p_rep.rdata(realData::density)     = p.rdata(realData::density);
                    p_rep.rdata(realData::mass)        = p.rdata(realData::mass);
                    p_rep.rdata(realData::oneOverI)    = p.rdata(realData::oneOverI);
                    p_rep.rdata(realData::radius)      = p.rdata(realData::radius);
                    p_rep.rdata(realData::omegax)      = p.rdata(realData::omegax);
                    p_rep.rdata(realData::omegay)      = p.rdata(realData::omegay);
                    p_rep.rdata(realData::omegaz)      = p.rdata(realData::omegaz);
                    p_rep.rdata(realData::statwt)      = p.rdata(realData::statwt);
                    p_rep.rdata(realData::dragcoeff)   = p.rdata(realData::dragcoeff);
                    p_rep.rdata(realData::dragx)       = p.rdata(realData::dragx);
                    p_rep.rdata(realData::dragy)       = p.rdata(realData::dragy);
                    p_rep.rdata(realData::dragz)       = p.rdata(realData::dragz);
                    p_rep.rdata(realData::c_ps)        = p.rdata(realData::c_ps);
                    p_rep.rdata(realData::temperature) = p.rdata(realData::temperature);
                    p_rep.rdata(realData::convection)  = p.rdata(realData::convection);

                    // Set id and cpu for this particle
                    p_rep.id()  = nextID + n;
                    p_rep.cpu() = myProc;

                    // Runtime added variables -- species mass fractions
                    for (int n_s(idx_X); n_s < (idx_X+nspecies_s); n_s++) {
                      // Copy data from particle to replicated one
                      ptile_data.m_runtime_rdata[n_s][index_repl] =
                        ptile_data.m_runtime_rdata[n_s][index];
                    }

                    // Runtime added variables -- species mass transfer rate per
                    // reaction
                    for (int n_s(0); n_s < nspecies_s; n_s++) {
                      for (int q(idx_G); q < (idx_G+nreactions); q++) {
                        // Set the index for the q-th mass txfr rate in the SoA
                        const int txfr_idx = q + n_s * nreactions;

                        // Copy data from particle to replicated one
                        ptile_data.m_runtime_rdata[txfr_idx][index_repl] =
                          ptile_data.m_runtime_rdata[txfr_idx][index];
                      }
                    }

                }); // p

                ParticleType::NextID(nextID + np);
            } // i

            amrex::Gpu::synchronize();

        } // pti

        Redistribute();
    } // idim

}
