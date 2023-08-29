#include <mfix_pc.H>
#include <mfix_des_K.H>
#include <mfix_pc_interactions_K.H>
#include <mfix_pc_updates_K.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_bc.H>
#include <mfix_solvers.H>
#include <mfix_monitors.H>
#include <mfix_calc_cell.H>

using namespace amrex;
using namespace Solvers;

void
MFIXParticleContainer::
ImposeMean (int const a_lev)
{
  BL_PROFILE("mfix::des::ImoseMean()");

  ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
  ReduceData<Real, Real, Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  std::vector<Real> sumVect(4, 0.);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIXParIter pti(*this, a_lev); pti.isValid(); ++pti) {

    const int np = NumberOfParticles(pti);

    sumVect[3] += static_cast<Real>(np);

    SoA& soa = pti.GetStructOfArrays();
    auto p_realarray = soa.realarray();

    reduce_op.eval(np, reduce_data, [p_realarray]
      AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
    {
      return {p_realarray[SoArealData::velx][p_id],
              p_realarray[SoArealData::vely][p_id],
              p_realarray[SoArealData::velz][p_id]};
    });
  }
  ReduceTuple host_tuple = reduce_data.value();

  sumVect[0] = amrex::get<0>(host_tuple);
  sumVect[1] = amrex::get<1>(host_tuple);
  sumVect[2] = amrex::get<2>(host_tuple);

  ParallelDescriptor::ReduceRealSum(sumVect.data(), sumVect.size());

  if (sumVect[3] > 0.) {

    const RealVect meanVel( sumVect[0]/sumVect[3] + m_constraint[0],
                            sumVect[1]/sumVect[3] + m_constraint[1],
                            sumVect[2]/sumVect[3] + m_constraint[2]);

    Print() << "Scaling particle velocity by the mean.\n";

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIXParIter pti(*this, a_lev); pti.isValid(); ++pti) {

      const int np = NumberOfParticles(pti);

      SoA& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      IntVect& apply = m_use_constraint;

      amrex::ParallelFor(np, [p_realarray, meanVel, apply]
      AMREX_GPU_DEVICE (int p_id) noexcept
      {
        if(apply[0]) { p_realarray[SoArealData::velx][p_id] -= meanVel[0]; }
        if(apply[1]) { p_realarray[SoArealData::vely][p_id] -= meanVel[1]; }
        if(apply[2]) { p_realarray[SoArealData::velz][p_id] -= meanVel[2]; }
      });
    }

  }
}
