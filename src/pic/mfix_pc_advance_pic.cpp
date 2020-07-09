#include <AMReX_EBFArrayBox.H>
#include <AMReX_Box.H>

#include <mfix.H>
#include <MFIX_BC_Parms.H>
#include <MFIX_PIC_Parms.H>

#include <math.h>

void MFIXParticleContainer::MFIX_PC_AdvanceParcels (amrex::Real dt, amrex::RealVect& gravity,
                                                    amrex::Vector< amrex::EBFArrayBoxFactory* > particle_ebfactory,
                                                    amrex::Vector< amrex::MultiFab* >& ep_s_in,
                                                    amrex::Vector< amrex::MultiFab* >& avg_prop_in,
                                                    amrex::Vector< amrex::MultiFab* >& cost,
                                                    std::string& knapsack_weight_type)
{

  BL_PROFILE("MFIXParticleContainer::MFIX_PC_AdvanceParcels()");

  for (int lev = 0; lev < nlev; lev ++ )
  {

    const auto dxi_array = Geom(lev).InvCellSizeArray();

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {

      // Timer used for load-balancing
      amrex::Real wt = ParallelDescriptor::second();

      PairIndex index(pti.index(), pti.LocalTileIndex());

      const int nrp = GetParticles(lev)[index].numRealParticles();

      auto& plev = GetParticles(lev);
      auto& ptile = plev[index];
      auto& aos   = ptile.GetArrayOfStructs();
      ParticleType* pstruct = aos().dataPtr();

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR("pic_time_march()", pic_time_march);
#endif
      /********************************************************************
       * Move particles based on collision forces and torques             *
       *******************************************************************/

      const auto p_lo = Geom(lev).ProbLoArray();
      const auto p_hi = Geom(lev).ProbHiArray();

      int x_lo_bc = BC::domain_bc[0];
      int x_hi_bc = BC::domain_bc[1];
      int y_lo_bc = BC::domain_bc[2];
      int y_hi_bc = BC::domain_bc[3];
      int z_lo_bc = BC::domain_bc[4];
      int z_hi_bc = BC::domain_bc[5];

      const auto& avg_prop_array = avg_prop_in[lev]->array(pti);

      const amrex::Real small_number = std::numeric_limits<Real>::epsilon();

      const amrex::Real velfac = PIC::velfac;
      const amrex::Real en = (PIC::damping_factor + 1.0);

      amrex::ParallelFor(nrp,
         [pstruct,dt,gravity,small_number,p_hi,p_lo, dxi_array, avg_prop_array,
          en,x_lo_bc,x_hi_bc,y_lo_bc,y_hi_bc,z_lo_bc,z_hi_bc,velfac]
        AMREX_GPU_DEVICE (int ip) noexcept
        {
          ParticleType& p = pstruct[ip];

          // position
          const amrex::RealVect pos = p.pos();

          // cell containing particle
          const int i = static_cast<int>(amrex::Math::floor((pos[0] - p_lo[0])*dxi_array[0]));
          const int j = static_cast<int>(amrex::Math::floor((pos[1] - p_lo[1])*dxi_array[1]));
          const int k = static_cast<int>(amrex::Math::floor((pos[2] - p_lo[2])*dxi_array[2]));

          // solids stress gradient:
          // grad_tau_p = (volume * grad_tau_p) / (mass * ep_s))
          amrex::RealVect grad_tau_p;
          grad_tau_p[0] = p.rdata(realData::omegax);
          grad_tau_p[1] = p.rdata(realData::omegay);
          grad_tau_p[2] = p.rdata(realData::omegaz);

          // inverse of particle mass and macroscopic density
          const amrex::Real inv_mass = 1.0 / p.rdata(realData::mass);

          // scale wrt drag coefficient and step size
          const amrex::Real scale = 1.0 / (1.0 + dt*p.rdata(realData::dragcoeff)*inv_mass);

          // updated parcel velocity devoid of the particle normal stress
          amrex::RealVect vel;
          vel[0] = scale*(p.rdata(realData::velx) + dt*(p.rdata(realData::dragx)*inv_mass + gravity[0]));
          vel[1] = scale*(p.rdata(realData::vely) + dt*(p.rdata(realData::dragy)*inv_mass + gravity[1]));
          vel[2] = scale*(p.rdata(realData::velz) + dt*(p.rdata(realData::dragz)*inv_mass + gravity[2]));

          // Slip velocity between the parcel and the bulk. This needs to be the
          // tentative new velocity, not the current velocity.
          amrex::RealVect slip_vel;
          slip_vel[0] = velfac*avg_prop_array(i,j,k,0) - vel[0];
          slip_vel[1] = velfac*avg_prop_array(i,j,k,1) - vel[1];
          slip_vel[2] = velfac*avg_prop_array(i,j,k,2) - vel[2];

          for( int dir(0); dir<3; dir++)
          {
            // solids stress velocity contribution:
            // -dt*( (volume * grad_tau_p) / (mass * ep_s)) / (1 + dt*beta/mass)
            amrex::Real up_prime_tau = -scale*dt*grad_tau_p[dir];

            // Add in contribution from stress gradient.
            // 1) Negative gradient. Push parcel in positive direction.
            // 2) Positive gradient. Push parcel in negative direction.
            // Both cases are limited limited by the slip velocity with the bulk.
            if(grad_tau_p[dir] <= 0.0 ){
              vel[dir] += amrex::max(0.0, amrex::min(up_prime_tau, en*slip_vel[dir]));
            } else {
              vel[dir] += amrex::min(0.0, amrex::max(up_prime_tau, en*slip_vel[dir]));
            }
          }

          // update parcel velocity
          p.rdata(realData::velx) = vel[0];
          p.rdata(realData::vely) = vel[1];
          p.rdata(realData::velz) = vel[2];

          // move the parcels
          p.pos(0) += dt * p.rdata(realData::velx);
          p.pos(1) += dt * p.rdata(realData::vely);
          p.pos(2) += dt * p.rdata(realData::velz);

          // Impose walls at domain extents
          amrex::Real eff_rad = p.rdata(realData::radius)*
            pow(p.rdata(realData::statwt),0.333);

          if (x_lo_bc && p.pos(0) < p_lo[0]+eff_rad)
          {
              p.pos(0) = p_lo[0] + eff_rad + small_number;
              p.rdata(realData::velx) = -p.rdata(realData::velx);
          }
          if (x_hi_bc && p.pos(0) + eff_rad > p_hi[0])
          {
             p.pos(0) = p_hi[0] - eff_rad - small_number;
             p.rdata(realData::velx) = -p.rdata(realData::velx);
          }
          if (y_lo_bc && p.pos(1) < p_lo[1])
          {
              p.pos(1) = p_lo[1] + small_number;
              p.rdata(realData::vely) = -p.rdata(realData::vely);
          }
          if (y_hi_bc && p.pos(1) > p_hi[1])
          {
             p.pos(1) = p_hi[1] - small_number;
             p.rdata(realData::vely) = -p.rdata(realData::vely);
          }
          if (z_lo_bc && p.pos(2) < p_lo[2])
          {
             p.pos(2) = p_lo[2] + small_number;
             p.rdata(realData::velz) = -p.rdata(realData::velz);
          }
          if (z_hi_bc && p.pos(2) > p_hi[2])
          {
             p.pos(2) = p_hi[2] - small_number;
             p.rdata(realData::velz) = -p.rdata(realData::velz);
          }
        });

      Gpu::synchronize();

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR_STOP(pic_time_march);
#endif

      /********************************************************************
       * Update runtime cost (used in load-balancing)                     *
       *******************************************************************/
      if (cost[lev])
      {
        // Runtime cost is either (weighted by tile box size):
        //   * time spent
        //   * number of particles
        const Box& tbx = pti.tilebox();
        if (knapsack_weight_type == "RunTimeCosts")
        {
          wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
        }
        else if (knapsack_weight_type == "NumParticles")
        {
          wt = nrp / tbx.d_numPts();
        }
        (*cost[lev])[pti].plus<RunOn::Device>(wt, tbx);
      }


    } // particle-tile iterator
  } // loop over levels
}
