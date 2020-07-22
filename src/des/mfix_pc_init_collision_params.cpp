#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <mfix_pc.H>

#include <math.h>

#include <mfix_dem_parms.H>

using namespace amrex;
using namespace std;


void MFIXParticleContainer::MFIX_PC_InitCollisionParams ()
{

   // amrex::Real max_dp[10], max_ro[10];
   amrex::Real avg_dp[10], avg_ro[10];

   // The number of phases was previously hard set at 10, however lowering
   //  this number would make this code faster.
   int num_of_phases_in_use = DEM::NPHASE;

   // Cycle through the different phases, starting from 1
   for (int phse = 1; phse <= num_of_phases_in_use; ++phse)
   {
      amrex::Real p_num  = 0.0; //number of particle
      amrex::Real p_diam = 0.0; //particle diameters
      amrex::Real p_dens = 0.0; //particle density

      amrex::Real max_diam = -1.0e32;
      amrex::Real max_den  = -1.0e32;

      for (int lev = 0; lev < nlev; lev++)
      {
#ifdef _OPENMP
#pragma omp parallel reduction(+:p_num, p_diam, p_dens) if (Gpu::notInLaunchRegion())
#endif
         for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
         {
            auto& particles = pti.GetArrayOfStructs();

            Gpu::HostVector<ParticleType> host_particles(pti.numParticles());
            Gpu::copy(Gpu::deviceToHost, particles.begin(), particles.end(), host_particles.begin());

            for (const auto& p: host_particles){
               if ( phse==p.idata(intData::phase) )
               {
                  p_num  += 1.0;
                  p_diam += p.rdata(realData::radius) * 2.0;
                  p_dens += p.rdata(realData::density);

                  max_diam = amrex::max(max_diam, p.rdata(realData::radius) * 2.0 );
                  max_den  = amrex::max(max_den,  p.rdata(realData::density) );
               }
            }
        }
     }

      // A single MPI call passes all three variables
      ParallelDescriptor::ReduceRealSum({p_num,p_diam,p_dens});
      ParallelDescriptor::ReduceRealMax({max_diam, max_den});

      //calculate averages or set = zero if no particles of that phase
      if (p_num==0){
         avg_dp[phse-1] = 0.0;
         avg_ro[phse-1] = 0.0;

         // max_dp[phse-1] = 0.0;
         // max_ro[phse-1] = 0.0;

      } else {
         avg_dp[phse-1] = p_diam/p_num;
         avg_ro[phse-1] = p_dens/p_num;

         // max_dp[phse-1] = max_diam;
         // max_ro[phse-1] = max_den;
      }
   }

   amrex::Real max_max_dp(0.);
   for (int phse = 0; phse < num_of_phases_in_use; ++phse)
   {
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(avg_dp[phse] > 0.0,
         "Average particle diameter cannot be zero");
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(avg_ro[phse] > 0.0,
         "Average particle density cannot be zero");

      max_max_dp = amrex::max(max_max_dp, avg_dp[phse]);
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
