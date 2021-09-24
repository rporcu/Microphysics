#include <AMReX_AmrCore.H>
#include <AMReX_buildInfo.H>
#include <AMReX_Geometry.H>
#include <AMReX_AmrParticles.H>

#include <mfix.H>
#include <mfix_pc.H>


using MFIXParIter = MFIXParticleContainer::MFIXParIter;

using MFIXAmrParticleContainer = AmrParticleContainer<AoSrealData::count, AoSintData::count,
                                                      SoArealData::count, SoAintData::count,
                                                      PinnedArenaAllocator>;


void mfix::PIC_to_PIC(const int lev,
                      std::string& restart_file,
                      const int pic_factor)
{
  MFIXAmrParticleContainer accumulator(this);

  MFIXAmrParticleContainer auxiliary(this);
  auxiliary.Restart(restart_file, "particles");

//  {
//    long fin_np = 0;
//    for (MFIXParIter pti(auxiliary, lev); pti.isValid(); ++pti) {
//      long np = pti.numParticles();
//      fin_np += np;
//    }
//
//    ParallelDescriptor::ReduceLongSum(fin_np);
//    amrex::Print() << "Initial number of read PIC parcels on level "
//      << lev << ": " << fin_np << std::endl;
//  }

  for (int n(0); n < pic_factor; ++n) {

    // Add particles to the accumulator and calls redistribute (true)
    accumulator.addParticles(auxiliary, true);
  }

  // Copy particles from the accumulator to the MFIXParticleContainer
  pc->copyParticles(accumulator, true);

  for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {
    int np = pti.numParticles();

    auto& soa = pti.GetStructOfArrays();
    auto p_realarray = soa.realarray();

    auto& aos = pti.GetArrayOfStructs();
    auto pstruct = aos().dataPtr();

    ParallelForRNG(np, [p_realarray,pic_factor,pstruct]
      AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
    {
      auto& part = pstruct[i];

      const Real p_radius = p_realarray[SoArealData::radius][i];

      part.pos(0) += (2*amrex::Random(engine)-1)*p_radius;
      part.pos(1) += (2*amrex::Random(engine)-1)*p_radius;
      part.pos(2) += (2*amrex::Random(engine)-1)*p_radius;

      p_realarray[SoArealData::statwt][i] /= pic_factor;
    });
  }

  {
    long fin_np = 0;
    for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {
      long np = pti.numParticles();
      fin_np += np;
    }

    ParallelDescriptor::ReduceLongSum(fin_np);
    amrex::Print() << "Final number of PIC parcels on level "
      << lev << ": " << fin_np << std::endl;
    pc->setTotalNumParticles(fin_np);
  }
}
