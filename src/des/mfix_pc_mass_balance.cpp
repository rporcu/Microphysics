#include <mfix_pc.H>
#include <mfix_solids.H>
#include <mfix_reactions.H>

using namespace amrex;


void MFIXParticleContainer::ResetMassBalance (int const a_n)
{
  m_p_mass_accum[a_n] = m_p_mass_accum[a_n + MFIXSpecies::NMAX];

  m_p_mass_accum[a_n + MFIXSpecies::NMAX] = 0.;
  m_p_mass_inflow[a_n]  = 0.;
  m_p_mass_outflow[a_n] = 0.;
  m_p_mass_prod[a_n]    = 0.;
}


void MFIXParticleContainer::ComputeMassProduction (int const a_lev, Real const a_dt)
{

  const int nspecies_s = solids.nspecies();
  std::vector<Real> prod(nspecies_s, 0.); // Mass produced/consumed

  const int idx_mass_txfr = m_runtimeRealData.mass_txfr;

  // Ideally, we would do this in one loop, but since the number of spcies
  // is unknown, we loop over each one.
  for (int n(0); n < nspecies_s; ++n) {

    // Reduce sum operation for mass and production of n-th species
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIXParIter pti(*this, a_lev); pti.isValid(); ++pti) {

      auto& plev  = GetParticles(a_lev);
      PairIndex index(pti.index(), pti.LocalTileIndex());

      // variables added at runtime
      auto ptile_data = plev[index].getParticleTileData();

      // rate of production
      auto p_txfr_mass_n = ptile_data.m_runtime_rdata[idx_mass_txfr+n];

      const int np = NumberOfParticles(pti);

      reduce_op.eval(np, reduce_data, [p_txfr_mass_n]
      AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
      {
        return {p_txfr_mass_n[p_id]};
      });
    } // MFIXParIter

    ReduceTuple host_tuple = reduce_data.value();
    prod[n] = amrex::get<0>(host_tuple);

  } //loop over species

  // Global sum and copy to global variable
  ParallelDescriptor::ReduceRealSum(prod.data(), nspecies_s);
  for (int n=0; n < nspecies_s; ++n) {
    m_p_mass_prod[n] += a_dt*prod[n];
  }

}



void MFIXParticleContainer::ComputeMassAccum ( int const a_offset )
{

  const int nspecies_s = solids.nspecies();
  std::vector<Real> accum(nspecies_s, 0.);
  m_total_numparticle = 0;

  for (int lev(0); lev < nlev; lev++) {

    const int idx_X_sn = m_runtimeRealData.X_sn;

    // Ideally, we would do this in one loop, but since the number of spcies
    // is unknown, we loop over each one.
    for (int n(0); n < nspecies_s; ++n) {

      // Reduce sum operation for mass and production of n-th species
      ReduceOps<ReduceOpSum> reduce_op;
      ReduceData<Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

        auto& plev  = GetParticles(lev);
        PairIndex index(pti.index(), pti.LocalTileIndex());

        // SoA real variables
        auto& soa = pti.GetStructOfArrays();

        // particle mass
        auto p_mass = soa.GetRealData(SoArealData::mass).data();

        // variables added at runtime
        auto ptile_data = plev[index].getParticleTileData();

        // species mass fraction
        auto p_Xn = ptile_data.m_runtime_rdata[idx_X_sn+n];

        const int np = NumberOfParticles(pti);

        reduce_op.eval(np, reduce_data, [p_mass, p_Xn]
        AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
        {
          return {p_mass[p_id]*p_Xn[p_id]};
        });
      } // MFIXParIter

      ReduceTuple host_tuple = reduce_data.value();
      accum[n] = amrex::get<0>(host_tuple);

    } //loop over species

    // Get an accurate count of all the particles in the system
    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {
      m_total_numparticle += static_cast<long>(NumberOfParticles(pti));
    } // MFIXParIter

  }//lev

  // Global sum and copy to global variable
  ParallelDescriptor::ReduceLongSum(m_total_numparticle);
  ParallelDescriptor::ReduceRealSum(accum.data(), nspecies_s);
  for (int n=0; n < nspecies_s; ++n) {
    m_p_mass_accum[n + a_offset*MFIXSpecies::NMAX] += accum[n];
  }

}



void MFIXParticleContainer::ComputeMassOutflow (int const a_lev)
{
  const int nspecies_s = solids.nspecies();
  std::vector<Real> outflow(nspecies_s, 0.);

  // particle tiles and geometry of this level

  const auto p_lo = Geom(a_lev).ProbLoArray();
  const auto p_hi = Geom(a_lev).ProbHiArray();

  const int idx_X_sn = m_runtimeRealData.X_sn;

  // Ideally, we would do this in one loop, but since the number of spcies
  // is unknown, we loop over each one.
  for (int n(0); n < nspecies_s; ++n) {

    // Reduce sum operation for mass and production of n-th species
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIXParIter pti(*this, a_lev); pti.isValid(); ++pti) {

      PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& plev  = GetParticles(a_lev);

      auto& aos   = plev[index].GetArrayOfStructs();
      ParticleType* pstruct = aos().dataPtr();

      // SoA real variables
      auto& soa = pti.GetStructOfArrays();

      // particle mass
      auto p_mass = soa.GetRealData(SoArealData::mass).data();

      // variables added at runtime
      auto ptile_data = plev[index].getParticleTileData();

      // species mass fraction
      auto p_Xn = ptile_data.m_runtime_rdata[idx_X_sn+n];

      const int np = NumberOfParticles(pti);

      reduce_op.eval(np, reduce_data, [p_lo, p_hi, pstruct, p_mass, p_Xn]
      AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
      {
        auto p = pstruct[p_id];
        Real p_mass_Xn(0);
        if ( p.pos(0) < p_lo[0] || p.pos(0) > p_hi[0] ||
             p.pos(1) < p_lo[1] || p.pos(1) > p_hi[1] ||
             p.pos(2) < p_lo[2] || p.pos(2) > p_hi[2] ) {
          p_mass_Xn = p_mass[p_id]*p_Xn[p_id];
        }
        return {p_mass_Xn};
      });
    } // MFIXParIter

  ReduceTuple host_tuple = reduce_data.value();
  outflow[n] = amrex::get<0>(host_tuple);

  } //loop over species

  // Global sum and copy to global variable
  ParallelDescriptor::ReduceRealSum(outflow.data(), nspecies_s);
  for (int n=0; n < nspecies_s; ++n) {
    m_p_mass_outflow[n] += outflow[n];
  }


}
