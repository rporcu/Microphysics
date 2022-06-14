#include <mfix_bc_list.H>
#include <mfix_bc_parms.H>
#include <mfix_pc.H>
#include <mfix_dem_parms.H>

using namespace amrex;


void MFIXParticleContainer::MFIX_PC_InitCollisionParams ()
{
  Vector<Real> max_dp, max_ro;
  Vector<Real> avg_dp, avg_ro;

  max_dp.resize(solids.NTYPES);
  max_ro.resize(solids.NTYPES);
  avg_dp.resize(solids.NTYPES);
  avg_ro.resize(solids.NTYPES);

  // Cycle through the different phases, starting from 1
  for (const int& phase: solids.phases)
  {
    Real h_pnum  = 0;   //number of particle
    Real h_pdiam = 0.0; //particle diameters
    Real h_pdens = 0.0; //particle density
    Real h_maxdiam = -1.e32;
    Real h_maxdens = -1.e32;

    for (int lev = 0; lev < nlev; lev++)
    {
      // Reduce sum operation for np, pdiam, pdens, maxdiam, maxdens
      ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpMax, ReduceOpMax> reduce_op;
      ReduceData<Real, Real, Real, Real, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        const int np = particles.size();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();
        auto p_intarray = soa.intarray();

        reduce_op.eval(np, reduce_data, [p_realarray,p_intarray,phase]
            AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
        {
          Real l_pnum  = 0._rt;
          Real l_pdiam = 0._rt;
          Real l_pdens = 0._rt;
          Real l_maxdiam = -1.e32;
          Real l_maxdens = -1.e32;

          if (phase == p_intarray[SoAintData::phase][p_id])
          {
            const Real density  = p_realarray[SoArealData::density][p_id];
            const Real diameter = 2.0*p_realarray[SoArealData::radius][p_id];

            l_pnum  = 1._rt;
            l_pdiam = diameter;
            l_pdens = density;
            l_maxdiam = diameter;
            l_maxdens = density;
          }

          return {l_pnum, l_pdiam, l_pdens, l_maxdiam, l_maxdens};
        });
      }

      ReduceTuple host_tuple = reduce_data.value();
      h_pnum  += amrex::get<0>(host_tuple);
      h_pdiam += amrex::get<1>(host_tuple);
      h_pdens += amrex::get<2>(host_tuple);
      h_maxdiam = amrex::max(h_maxdiam, amrex::get<3>(host_tuple));
      h_maxdens = amrex::max(h_maxdens, amrex::get<4>(host_tuple));
    }

    // A single MPI call passes all three variables
    ParallelDescriptor::ReduceRealSum({h_pnum,h_pdiam,h_pdens});
    ParallelDescriptor::ReduceRealMax({h_maxdiam,h_maxdens});

    //calculate averages or set = zero if no particles of that phase
    const int phase_idx = SolidsPhase::phase_to_index(phase);

    if (h_pnum==0) {
       avg_dp[phase_idx] = 0.0;
       avg_ro[phase_idx] = 0.0;

       max_dp[phase_idx] = 0.0;
       max_ro[phase_idx] = 0.0;

    } else {
       avg_dp[phase_idx] = h_pdiam/h_pnum;
       avg_ro[phase_idx] = h_pdens/h_pnum;

       max_dp[phase_idx] = h_maxdiam;
       max_ro[phase_idx] = h_maxdens;
    }
  }

  // Loop over BCs
  for (int bcv(0); bcv < BC::bc.size(); ++bcv) {

    // EB flow with at least one solid
    if (BC::bc[bcv].type == BCList::eb && BC::bc[bcv].solids.size() > 0) {

      const Real tolerance = std::numeric_limits<Real>::epsilon();

      for(int lcs(0); lcs < BC::bc[bcv].solids.size(); lcs++) {

        const SOLIDS_t& bc_solid = BC::bc[bcv].solids[lcs];

        if(bc_solid.volfrac > tolerance) {

          const Real mean_dp_bc = bc_solid.diameter.get_mean();
          const Real  max_dp_bc = bc_solid.diameter.is_constant() ?
                                    mean_dp_bc : bc_solid.diameter.get_max();

          const Real mean_rhop_bc = bc_solid.density.get_mean();
          const Real  max_rhop_bc = bc_solid.density.is_constant() ?
                                      mean_rhop_bc : bc_solid.density.get_max();

          if ( avg_dp[lcs] == 0.0 )
            avg_dp[lcs] = mean_dp_bc;

          if ( avg_ro[lcs] == 0.0 )
            avg_ro[lcs] = mean_rhop_bc;

          max_dp[lcs] = amrex::max(max_dp[lcs], max_dp_bc);
          max_ro[lcs] = amrex::max(max_ro[lcs], max_rhop_bc);

        }
      }
    }
  }


  Real max_max_dp(0.);
  for (const int& phase: solids.phases)
  {
    const int phase_idx = SolidsPhase::phase_to_index(phase);

     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(avg_dp[phase_idx] > 0.0,
        "Average particle diameter cannot be zero");
     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_dp[phase_idx] > 0.0,
        "Maximum particle diameter cannot be zero");
     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(avg_ro[phase_idx] > 0.0,
        "Average particle density cannot be zero");
     AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_ro[phase_idx] > 0.0,
        "Maximum particle density cannot be zero");

     max_max_dp = amrex::max(max_max_dp, max_dp[phase_idx]);
  }

  // (3*max_dp/2)^2
  DEM::neighborhood = 2.25*max_max_dp*max_max_dp;

  Real tcoll(1.0);

  DEM::A2D::array_type host_etan, host_etat, host_en;
  DEM::A1D::array_type host_etan_w, host_etat_w, host_en_w;
#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy_async(&host_en, DEM::en.arrayPtr(), sizeof(DEM::A2D::array_type));
  amrex::Gpu::dtoh_memcpy_async(&host_en_w, DEM::en_w.arrayPtr(), sizeof(DEM::A1D::array_type));
  amrex::Gpu::synchronize();
#else
  host_en = *(DEM::en.arrayPtr());
  host_en_w = *(DEM::en_w.arrayPtr());
#endif

  for (const int& phase: solids.phases) {

    const int phase_idx = SolidsPhase::phase_to_index(phase);

    const Real dp_m = avg_dp[phase_idx];
    const Real mass_m = (M_PI/6.0)*(dp_m*dp_m*dp_m) * avg_ro[phase_idx];

    // Collision parameters between phase_idx M and all other phase_idxs
    for (int other_phase_idx(phase_idx); other_phase_idx < solids.NTYPES; ++other_phase_idx) {

      const Real dp_l = avg_dp[other_phase_idx];
      const Real mass_l = (M_PI/6.0)*(dp_l*dp_l*dp_l) * avg_ro[other_phase_idx];

      const Real mass_eff = (mass_m * mass_l) / (mass_m + mass_l);

      //Calculate the M-L normal and tangential damping coefficients
      host_etan(phase_idx,other_phase_idx) = 2.0*std::sqrt(DEM::kn*mass_eff);
      if(amrex::Math::abs(host_en(phase_idx,other_phase_idx)) > 0.0){
        const Real log_en = std::log(host_en(phase_idx,other_phase_idx));
        host_etan(phase_idx,other_phase_idx) *= amrex::Math::abs(log_en)/(std::sqrt(M_PI*M_PI + log_en*log_en));
      }
      host_etat(phase_idx,other_phase_idx) = DEM::eta_fac * host_etan(phase_idx,other_phase_idx);

      // Store the symmetric components
      host_etan(other_phase_idx,phase_idx) = host_etan(phase_idx,other_phase_idx);
      host_etat(other_phase_idx,phase_idx) = host_etat(phase_idx,other_phase_idx);

      // Collision time scale for M-L interactions
      Real tcoll_ml = M_PI/sqrt(DEM::kn/mass_eff -
         0.25*(host_etan(phase_idx,other_phase_idx)/mass_eff)*(host_etan(phase_idx,other_phase_idx)/mass_eff));

      tcoll = amrex::min(tcoll, tcoll_ml);
    }

    // Collision parameters between phase_idx M and wall
    {
      const Real mass_eff = mass_m;

      //Calculate the M-L normal and tangential damping coefficients
      host_etan_w(phase_idx) = 2.0*std::sqrt(DEM::kn_w*mass_eff);
      if(amrex::Math::abs(host_en_w(phase_idx)) > 0.0){
         const Real log_en = std::log(host_en_w(phase_idx));
         host_etan_w(phase_idx) *= amrex::Math::abs(log_en)/(std::sqrt(M_PI*M_PI + log_en*log_en));
      }
      host_etat_w(phase_idx) = DEM::eta_w_fac * host_etan_w(phase_idx);

      // Calculate the collision time scale.
      Real tcoll_w = M_PI/std::sqrt(DEM::kn_w/mass_eff -
         0.25*(host_etan_w(phase_idx)/mass_eff)*(host_etan_w(phase_idx)/mass_eff));

      tcoll = amrex::min(tcoll, tcoll_w);
    }
  }

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy_async(DEM::etan.arrayPtr(), &host_etan, sizeof(DEM::A2D::array_type));
  amrex::Gpu::htod_memcpy_async(DEM::etat.arrayPtr(), &host_etat, sizeof(DEM::A2D::array_type));
  amrex::Gpu::htod_memcpy_async(DEM::etan_w.arrayPtr(), &host_etan_w, sizeof(DEM::A1D::array_type));
  amrex::Gpu::htod_memcpy_async(DEM::etat_w.arrayPtr(), &host_etat_w, sizeof(DEM::A1D::array_type));
  amrex::Gpu::synchronize();
#else
  *(DEM::etan.arrayPtr()) = host_etan;
  *(DEM::etat.arrayPtr()) = host_etat;
  *(DEM::etan_w.arrayPtr()) = host_etan_w;
  *(DEM::etat_w.arrayPtr()) = host_etat_w;
#endif

  ParmParse pp("mfix");
  Real tcoll_ratio = 50.;
  pp.query("tcoll_ratio", tcoll_ratio);

  DEM::dtsolid = tcoll / tcoll_ratio;
}
