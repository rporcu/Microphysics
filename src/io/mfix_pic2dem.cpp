#include <AMReX_AmrCore.H>
#include <AMReX_buildInfo.H>
#include <AMReX_Geometry.H>
#include <AMReX_AmrParticles.H>

#include <mfix.H>
#include <mfix_pc.H>


using MFIXParIter = MFIXParticleContainer::MFIXParIter;

void mfix::PIC_to_DEM(const int lev)
{
  long fin_np = 0;
  for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {
    long np = pti.numParticles();
    fin_np += np;
  }

  ParallelDescriptor::ReduceLongSum(fin_np);
  amrex::Print() << "Initial number of read PIC parcels on level "
    << lev << ": " << fin_np << std::endl;

  pc->clearNeighbors();
  pc->Redistribute(0, 0, 0, 0);
  pc->fillNeighbors();

  pc->buildNeighborList(MFIXCheckFullPair(m_dem.neighborhood()), false);

  for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {
    
    MFIXParticleContainer::PairIndex index(pti.index(), pti.LocalTileIndex());

    auto& plev  = pc->GetParticles(lev);
    auto& ptile = plev[index];
    auto& aos   = ptile.GetArrayOfStructs();
    MFIXParticleContainer::ParticleType* pstruct = aos().dataPtr();

    const int nrp = pc->GetParticles(lev)[index].numRealParticles();

    auto& soa = ptile.GetStructOfArrays();
    auto p_realarray = soa.realarray();
    auto p_intarray = soa.intarray();

    auto m_neighbor_list = pc->get_neighbor_list();
    auto nbor_data = m_neighbor_list[lev][index].data();

    amrex::ParallelFor(nrp, [p_intarray,p_realarray,pstruct,nbor_data]
      AMREX_GPU_DEVICE (int i) noexcept
    {
      const Real p1statwt = p_realarray[SoArealData::statwt][i];
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(p1statwt == 1.,
        "A PIC parcel with statwt != 1 cannot be converted to a DEM particle");

      auto& particle = pstruct[i];

      const Real p1radius = p_realarray[SoArealData::radius][i];
      const Real p1mass = p_realarray[SoArealData::mass][i];

      const Real p1I = p1mass * (4*p1radius*p1radius) / 10.;

      p_realarray[SoArealData::oneOverI][i] = 1. / p1I;
      p_realarray[SoArealData::omegax][i] = 0.;
      p_realarray[SoArealData::omegay][i] = 0.;
      p_realarray[SoArealData::omegaz][i] = 0.;

      RealVect pos1(particle.pos());

      Real max_overlap(-1*p1radius);

      const auto neighbs = nbor_data.getNeighbors(i);
      for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit) {
        const int j = mit.index();
        const auto p2 = *mit;

        const RealVect pos2 = p2.pos();
        const Real p2radius = p_realarray[SoArealData::radius][j];

        const Real distance = (pos1 - pos2).vectorLength();

        const Real min_distance = p1radius + p2radius;

        const Real overlap = min_distance - distance;
        max_overlap = amrex::max(max_overlap, overlap);
      }

      if (max_overlap > 0.) {
        p_intarray[SoAintData::state][i] = 10;
      }
    });
  }

  pc->Redistribute();

  Gpu::synchronize();

  fin_np = 0;
  for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {
    long np = pti.numParticles();
    fin_np += np;
  }

  ParallelDescriptor::ReduceLongSum(fin_np);
  amrex::Print() << "Final number of DEM particles on level "
    << lev << ": " << fin_np << std::endl;
  pc->setTotalNumParticles(fin_np);
}
