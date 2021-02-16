#include <mfix_pc.H>
#include <mfix_dem_parms.H>

using namespace amrex;


void MFIXParticleContainer::MFIX_PC_InitCollisionParams ()
{

   amrex::Real max_dp[10], max_ro[10];
   amrex::Real avg_dp[10], avg_ro[10];

   // The number of phases was previously hard set at 10, however lowering
   //  this number would make this code faster.
   int num_of_phases_in_use = DEM::NPHASE;

   // Cycle through the different phases, starting from 1
   for (int phase = 1; phase <= num_of_phases_in_use; ++phase)
   {
      Real h_pnum  = 0;   //number of particle
      Real h_pdiam = 0.0; //particle diameters
      Real h_pdens = 0.0; //particle density
      Real h_maxdiam = -1.e32;
      Real h_maxdens = -1.e32;

      for (int lev = 0; lev < nlev; lev++)
      {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
          // Reduce sum operation for np, pdiam, pdens, maxdiam, maxdens
          ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpMax, ReduceOpMax> reduce_op;
          ReduceData<Real, Real, Real, Real, Real> reduce_data(reduce_op);
          using ReduceTuple = typename decltype(reduce_data)::Type;

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
        else
#endif
        {
#ifdef _OPENMP
#pragma omp parallel reduction(+:h_pnum,h_pdiam,h_pdens) \
                     reduction(max:h_maxdiam,h_maxdens) \
                     if (Gpu::notInLaunchRegion())
#endif
          for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
          {
            auto& particles = pti.GetArrayOfStructs();
            const int np = particles.size();

            auto& soa = pti.GetStructOfArrays();
            auto p_realarray = soa.realarray();
            auto p_intarray = soa.intarray();

            for (int p_id(0); p_id < np; ++p_id)
            {
              if (phase == p_intarray[SoAintData::phase][p_id])
              {
                const Real density  = p_realarray[SoArealData::density][p_id];
                const Real diameter = 2.0*p_realarray[SoArealData::radius][p_id];

                h_pnum  += 1.0;
                h_pdiam += diameter;
                h_pdens += density;
                h_maxdiam = amrex::max(h_maxdiam, diameter);
                h_maxdens = amrex::max(h_maxdens, density);
              }
            }
          }
        }
      }

      // A single MPI call passes all three variables
      ParallelDescriptor::ReduceRealSum({h_pnum,h_pdiam,h_pdens});
      ParallelDescriptor::ReduceRealMax({h_maxdiam,h_maxdens});

      //calculate averages or set = zero if no particles of that phase
      if (h_pnum==0){
         avg_dp[phase-1] = 0.0;
         avg_ro[phase-1] = 0.0;

         max_dp[phase-1] = 0.0;
         max_ro[phase-1] = 0.0;

      } else {
         avg_dp[phase-1] = h_pdiam/h_pnum;
         avg_ro[phase-1] = h_pdens/h_pnum;

         max_dp[phase-1] = h_maxdiam;
         max_ro[phase-1] = h_maxdens;
      }
   }

   amrex::Real max_max_dp(0.);
   for (int phase = 0; phase < num_of_phases_in_use; ++phase)
   {
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(avg_dp[phase] > 0.0,
         "Average particle diameter cannot be zero");
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_dp[phase] > 0.0,
         "Maximum particle diameter cannot be zero");
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(avg_ro[phase] > 0.0,
         "Average particle density cannot be zero");
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_ro[phase] > 0.0,
         "Maximum particle density cannot be zero");

      max_max_dp = amrex::max(max_max_dp, max_dp[phase]);
   }

   // (3*max_dp/2)^2
   DEM::neighborhood = 2.25*max_max_dp*max_max_dp;

   amrex::Real tcoll(1.0);

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

   for (int m(0); m < num_of_phases_in_use; ++m)
   {
      const amrex::Real dp_m = avg_dp[m];
      const amrex::Real mass_m = (M_PI/6.0)*(dp_m*dp_m*dp_m) * avg_ro[m];

      // Collision parameters between type M and all other types
      for (int l(m); l < num_of_phases_in_use; ++l)
      {

         const amrex::Real dp_l = avg_dp[l];
         const amrex::Real mass_l = (M_PI/6.0)*(dp_l*dp_l*dp_l) * avg_ro[l];

         const amrex::Real mass_eff = (mass_m * mass_l) / (mass_m + mass_l);

         //Calculate the M-L normal and tangential damping coefficients
         host_etan(m,l) = 2.0*std::sqrt(DEM::kn*mass_eff);
         if(amrex::Math::abs(host_en(m,l)) > 0.0){
            const amrex::Real log_en = std::log(host_en(m,l));
            host_etan(m,l) *= amrex::Math::abs(log_en)/(std::sqrt(M_PI*M_PI + log_en*log_en));
         }
         host_etat(m,l) = DEM::eta_fac * host_etan(m,l);

         // Store the symmetric components
         host_etan(l,m) = host_etan(m,l);
         host_etat(l,m) = host_etat(m,l);

         // Collision time scale for M-L interactions
         amrex::Real tcoll_ml = M_PI/sqrt(DEM::kn/mass_eff -
            0.25*(host_etan(m,l)/mass_eff)*(host_etan(m,l)/mass_eff));

         tcoll = amrex::min(tcoll, tcoll_ml);
      }

      // Collision parameters between type M and wall
      {
         const amrex::Real mass_eff = mass_m;

         //Calculate the M-L normal and tangential damping coefficients
         host_etan_w(m) = 2.0*std::sqrt(DEM::kn_w*mass_eff);
         if(amrex::Math::abs(host_en_w(m)) > 0.0){
            const amrex::Real log_en = std::log(host_en_w(m));
            host_etan_w(m) *= amrex::Math::abs(log_en)/(std::sqrt(M_PI*M_PI + log_en*log_en));
         }
         host_etat_w(m) = DEM::eta_w_fac * host_etan_w(m);

        // Calculate the collision time scale.
        amrex::Real tcoll_w = M_PI/std::sqrt(DEM::kn_w/mass_eff -
           0.25*(host_etan_w(m)/mass_eff)*(host_etan_w(m)/mass_eff));

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
   amrex::Real tcoll_ratio = 50.;
   pp.query("tcoll_ratio", tcoll_ratio);

   DEM::dtsolid = tcoll / tcoll_ratio;
}
