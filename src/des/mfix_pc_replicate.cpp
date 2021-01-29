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

    const int idx_X = SoAspeciesData::X_sn * nspecies_s;

    const int idx_G = SoAspeciesData::count * nspecies_s +
                      SoAreactionsData::G_sn_pg_q * nreactions;

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

            auto& soa = ptile.GetStructOfArrays();
            auto p_realarray = soa.realarray();
            auto p_intarray = soa.intarray();

            //Access to added variables
            auto ptile_data = ptile.getParticleTileData();

            for (int i = 0; i < Nrep[idim]; i++)
            {
                if (i == 0) continue; // skip the first copy in each direction

                RealVect shift = {0., 0., 0.};
                shift[idim] = i * orig_domain_size[idim];

                const int nextID = ParticleType::NextID();

                amrex::ParallelFor(np, [pstruct,p_realarray,p_intarray,ptile_data,np,
                    nextID,myProc,nspecies_s,nreactions,idx_X,idx_G,shift,i]
                  AMREX_GPU_DEVICE (int n) noexcept
                {
                    int index = n;
                    int index_repl = i*np + n;

                    ParticleType p = pstruct[index];
                    ParticleType& p_rep = pstruct[index_repl];

                    p_rep.pos(0) = p.pos(0) + shift[0];
                    p_rep.pos(1) = p.pos(1) + shift[1];
                    p_rep.pos(2) = p.pos(2) + shift[2];

                    p_realarray[SoArealData::velx][index_repl] = p_realarray[SoArealData::velx][index];
                    p_realarray[SoArealData::vely][index_repl] = p_realarray[SoArealData::vely][index];
                    p_realarray[SoArealData::velz][index_repl] = p_realarray[SoArealData::velz][index];

                    // Set other particle properties
                    p_intarray[SoAintData::phase][index_repl] = p_intarray[SoAintData::phase][index];
                    p_realarray[SoArealData::volume][index_repl] = p_realarray[SoArealData::volume][index];
                    p_realarray[SoArealData::density][index_repl] = p_realarray[SoArealData::density][index];
                    p_realarray[SoArealData::mass][index_repl] = p_realarray[SoArealData::mass][index];
                    p_realarray[SoArealData::oneOverI][index_repl] = p_realarray[SoArealData::oneOverI][index];
                    p_realarray[SoArealData::radius][index_repl] = p_realarray[SoArealData::radius][index];
                    p_realarray[SoArealData::omegax][index_repl] = p_realarray[SoArealData::omegax][index];
                    p_realarray[SoArealData::omegay][index_repl] = p_realarray[SoArealData::omegay][index];
                    p_realarray[SoArealData::omegaz][index_repl] = p_realarray[SoArealData::omegaz][index];
                    p_realarray[SoArealData::statwt][index_repl] = p_realarray[SoArealData::statwt][index];
                    p_realarray[SoArealData::dragcoeff][index_repl] = p_realarray[SoArealData::dragcoeff][index];
                    p_realarray[SoArealData::dragx][index_repl] = p_realarray[SoArealData::dragx][index];
                    p_realarray[SoArealData::dragy][index_repl] = p_realarray[SoArealData::dragy][index];
                    p_realarray[SoArealData::dragz][index_repl] = p_realarray[SoArealData::dragz][index];
                    p_realarray[SoArealData::c_ps][index_repl] = p_realarray[SoArealData::c_ps][index];
                    p_realarray[SoArealData::temperature][index_repl] = p_realarray[SoArealData::temperature][index];
                    p_realarray[SoArealData::convection][index_repl] = p_realarray[SoArealData::convection][index];

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
