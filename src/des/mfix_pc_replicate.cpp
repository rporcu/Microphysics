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

    const int idx_X_sn = m_runtimeRealData.X_sn;
    const int idx_ro_sn_txfr = m_runtimeRealData.ro_sn_txfr;
    const int idx_vel_s_txfr = m_runtimeRealData.vel_s_txfr;
    const int idx_h_s_txfr = m_runtimeRealData.h_s_txfr;
    const int idx_count = m_runtimeRealData.count;

    for (int idim = 0; idim < 3; ++idim)
    {
        for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& ptile = DefineAndReturnParticleTile(lev,pti);
            int np = pti.numParticles();

            int np_replicated = np * (Nrep[idim] - 1);

            int new_np = np + np_replicated;
            ptile.resize(new_np);

            // Add runtime-added components
            const int start = AoSrealData::count + SoArealData::count;
            for (int comp(0); comp < m_runtimeRealData.count; ++comp)
              ptile.push_back_real(start+comp, np_replicated, 0.);

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

                amrex::ParallelFor(np, [np,pstruct,p_realarray,p_intarray,
                    ptile_data,nextID,myProc,nspecies_s,nreactions,idx_X_sn,
                    idx_ro_sn_txfr,idx_vel_s_txfr,idx_h_s_txfr,
                    idx_count,shift,i]
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
                    p_intarray[SoAintData::state][index_repl] = p_intarray[SoAintData::state][index];
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
                    p_realarray[SoArealData::cp_s][index_repl] = p_realarray[SoArealData::cp_s][index];
                    p_realarray[SoArealData::temperature][index_repl] = p_realarray[SoArealData::temperature][index];
                    p_realarray[SoArealData::convection][index_repl] = p_realarray[SoArealData::convection][index];

                    // Set id and cpu for this particle
                    p_rep.id()  = nextID + n;
                    p_rep.cpu() = myProc;

                    int start_idx = idx_X_sn;
                    int end_idx   = idx_ro_sn_txfr;

                    // Runtime added variables -- species mass fractions
                    for (int idx(start_idx); idx < end_idx; ++idx) {
                      // Copy data from particle to replicated one
                      ptile_data.m_runtime_rdata[idx][index_repl] = ptile_data.m_runtime_rdata[idx][index];
                    }

                    start_idx = end_idx;
                    end_idx   = idx_vel_s_txfr;

                    // Runtime added variables -- species mass txfr rates
                    for (int idx(start_idx); idx < end_idx; ++idx) {
                      // Copy data from particle to replicated one
                      ptile_data.m_runtime_rdata[idx][index_repl] = ptile_data.m_runtime_rdata[idx][index];
                    }

                    start_idx = end_idx;
                    end_idx   = idx_h_s_txfr;

                    // Runtime added variables -- species momentum txfr rate
                    for (int idx(start_idx); idx < end_idx; ++idx) {
                      // Copy data from particle to replicated one
                      ptile_data.m_runtime_rdata[idx][index_repl] = ptile_data.m_runtime_rdata[idx][index];
                    }

                    start_idx = end_idx;
                    end_idx   = idx_count;

                    // Runtime added variables -- species energy txfr rate
                    for (int idx(start_idx); idx < end_idx; ++idx) {
                      // Copy data from particle to replicated one
                      ptile_data.m_runtime_rdata[idx][index_repl] = ptile_data.m_runtime_rdata[idx][index];
                    }
                }); // p

                ParticleType::NextID(nextID + np);
            } // i

            amrex::Gpu::synchronize();

        } // pti

        Redistribute();
    } // idim

}
