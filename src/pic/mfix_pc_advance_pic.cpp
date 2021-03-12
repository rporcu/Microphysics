#include <mfix.H>
#include <mfix_bc_parms.H>


void MFIXParticleContainer::MFIX_PC_AdvanceParcels (Real dt, RealVect& gravity,
                                                    Vector< MultiFab* >& avg_prop_in,
                                                    Vector< MultiFab* >& cost,
                                                    std::string& knapsack_weight_type,
                                                    const int advect_enthalpy)
{

  BL_PROFILE("MFIXParticleContainer::MFIX_PC_AdvanceParcels()");

  const Real en = (PIC::damping_factor + 1.0);
  const Real velfac = PIC::velfac;

  const Real three_sqrt_two(3.0*std::sqrt(2.0));

  const int x_lo_bc = BC::domain_bc[0];
  const int x_hi_bc = BC::domain_bc[1];
  const int y_lo_bc = BC::domain_bc[2];
  const int y_hi_bc = BC::domain_bc[3];
  const int z_lo_bc = BC::domain_bc[4];
  const int z_hi_bc = BC::domain_bc[5];

  for (int lev = 0; lev < nlev; lev ++ )
  {

    Box domain(Geom(lev).Domain());

    const auto dxi = Geom(lev).InvCellSizeArray();

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {

      // Timer used for load-balancing
      Real wt = ParallelDescriptor::second();

      PairIndex index(pti.index(), pti.LocalTileIndex());

      const int nrp = GetParticles(lev)[index].numRealParticles();

      auto& plev = GetParticles(lev);
      auto& ptile = plev[index];
      auto& aos   = ptile.GetArrayOfStructs();
      ParticleType* pstruct = aos().dataPtr();

      auto& soa = ptile.GetStructOfArrays();
      auto p_realarray = soa.realarray();

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR("pic_time_march()", pic_time_march);
#endif
      /********************************************************************
       * Move particles based on collision forces and torques             *
       *******************************************************************/

      const auto p_lo = Geom(lev).ProbLoArray();
      const auto p_hi = Geom(lev).ProbHiArray();

      const auto& avg_prop_array = avg_prop_in[lev]->array(pti);

      const Real tolerance = std::numeric_limits<Real>::epsilon();

      amrex::ParallelFor(nrp,
        [pstruct,p_realarray,dt,gravity, p_hi,p_lo, dxi, tolerance,
         avg_prop_array, velfac, en, three_sqrt_two,
         x_lo_bc,x_hi_bc,y_lo_bc,y_hi_bc,z_lo_bc,z_hi_bc]
        AMREX_GPU_DEVICE (int lp) noexcept
        {
          ParticleType& p = pstruct[lp];

          // position
          const RealVect pos = p.pos();

          // solids volume fraction
          const Real eps_p = amrex::max(1.0e-8, p_realarray[SoArealData::oneOverI][lp]);

          // inverse of particle mass
          const Real inv_mass = 1.0 / p_realarray[SoArealData::mass][lp];

          // scale wrt drag coefficient and step size
          const Real scale = 1.0 / (1.0 + dt*p_realarray[SoArealData::dragcoeff][lp]*inv_mass);

          // solids stress gradient:
          // grad_tau_p = (volume * grad_tau_p) / (mass * ep_s))
          RealVect grad_tau_p;
          grad_tau_p[0] = p_realarray[SoArealData::omegax][lp];
          grad_tau_p[1] = p_realarray[SoArealData::omegay][lp];
          grad_tau_p[2] = p_realarray[SoArealData::omegaz][lp];

          // updated parcel velocity devoid of the particle normal stress
          RealVect vel;
          vel[0] = scale*(p_realarray[SoArealData::velx][lp] + dt*(p_realarray[SoArealData::dragx][lp]*inv_mass + gravity[0]));
          vel[1] = scale*(p_realarray[SoArealData::vely][lp] + dt*(p_realarray[SoArealData::dragy][lp]*inv_mass + gravity[1]));
          vel[2] = scale*(p_realarray[SoArealData::velz][lp] + dt*(p_realarray[SoArealData::dragz][lp]*inv_mass + gravity[2]));

          const Real lx = (pos[0] - p_lo[0])*dxi[0] + 0.5;
          const Real ly = (pos[1] - p_lo[1])*dxi[1] + 0.5;
          const Real lz = (pos[2] - p_lo[2])*dxi[2] + 0.5;

          const int i = static_cast<int>(amrex::Math::floor(lx));
          const int j = static_cast<int>(amrex::Math::floor(ly));
          const int k = static_cast<int>(amrex::Math::floor(lz));

          const Real wx_hi(lx - static_cast<Real>(i));
          const Real wy_hi(ly - static_cast<Real>(j));
          const Real wz_hi(lz - static_cast<Real>(k));

          const Real wx_lo(1.0 - wx_hi);
          const Real wy_lo(1.0 - wy_hi);
          const Real wz_lo(1.0 - wz_hi);

          int ii = i;
          if(vel[0]*grad_tau_p[0] > 0.) {
            ii = (vel[0]>0.) ? i+1 : i-1;
          }

          const Real u_face =
            wy_lo*wz_lo*avg_prop_array(ii,j-1,k-1,0) +
            wy_hi*wz_lo*avg_prop_array(ii,j  ,k-1,0) +
            wy_lo*wz_hi*avg_prop_array(ii,j-1,k  ,0) +
            wy_hi*wz_hi*avg_prop_array(ii,j  ,k  ,0);


          int jj = j;
          if(vel[1]*grad_tau_p[1] > 0.0) {
            jj = ( vel[1] > 0.0 ) ? j+1 : j-1;
          }

          const Real v_face =
            wx_lo*wz_lo*avg_prop_array(i-1,jj,k-1,1) +
            wx_hi*wz_lo*avg_prop_array(i  ,jj,k-1,1) +
            wx_lo*wz_hi*avg_prop_array(i-1,jj,k  ,1) +
            wx_hi*wz_hi*avg_prop_array(i  ,jj,k  ,1);


          int kk = k;
          if(vel[2]*grad_tau_p[2] > 0.0) {
            kk = ( vel[2] > 0.0 ) ? k+1 : k-1;
          }

          const Real w_face =
            wx_lo*wy_lo*avg_prop_array(i-1,j-1,kk,2) +
            wx_hi*wy_lo*avg_prop_array(i  ,j-1,kk,2) +
            wx_lo*wy_hi*avg_prop_array(i-1,j  ,kk,2) +
            wx_hi*wy_hi*avg_prop_array(i  ,j  ,kk,2);

          const RealVect avg_vel  = {u_face, v_face, w_face};

          const Real inv_rops = 1.0 / (eps_p*p_realarray[SoArealData::density][lp]);

          const Real vel_limit = p_realarray[SoArealData::radius][lp] / (dt * three_sqrt_two * eps_p);

          for( int dir(0); dir<3; dir++)
          {

            Real bulk_vel;

            if( amrex::Math::abs(avg_vel[dir]) > tolerance and
                amrex::Math::abs(    vel[dir]) > tolerance ){
              bulk_vel = 2./(1./avg_vel[dir] + 1./vel[dir]);
            } else {
              bulk_vel = 0.0;
            }

            // Slip velocity between the parcel and the bulk. This needs to be the
            // tentative new velocity, not the current velocity.
            const Real slip_vel = velfac*bulk_vel - vel[dir];

            // solids stress velocity contribution:
            // -dt*( grad_tau_p / (density * ep_s)) / (1 + dt*beta/mass)
            const Real del_up = -dt*scale*grad_tau_p[dir] * inv_rops;

            // Add in contribution from stress gradient.
            // 1) Negative gradient. Push parcel in positive direction.
            // 2) Positive gradient. Push parcel in negative direction.
            // Both cases are limited limited by the slip velocity with the bulk.
            if(grad_tau_p[dir] <= 0.0 ){
              vel[dir] += amrex::max(0.0, amrex::min(del_up, en*slip_vel));

            } else {
              vel[dir] += amrex::min(0.0, amrex::max(del_up, en*slip_vel));
            }

            // Limit parcel velocity based on the mean free path
            vel[dir] = ( vel[dir] > 0.0 ) ?
              amrex::min(vel[dir], vel_limit):
              amrex::max(vel[dir],-vel_limit);

          }


          // update parcel velocity
          p_realarray[SoArealData::velx][lp] = vel[0];
          p_realarray[SoArealData::vely][lp] = vel[1];
          p_realarray[SoArealData::velz][lp] = vel[2];

          // move the parcels
          p.pos(0) += dt * p_realarray[SoArealData::velx][lp];
          p.pos(1) += dt * p_realarray[SoArealData::vely][lp];
          p.pos(2) += dt * p_realarray[SoArealData::velz][lp];


          // Impose domain constraints. Note that we only make sure that a parcel's
          // centroid is inside. A later routine will take into account the wall
          // collision and redirection of velocity.

          // Take the particle radius as the offset. This should be (much) smaller
          // than the effective parcel radius. We only want to move the parcel a
          // little but it needs to be enough to keep it inside the domain.
          const Real offset = p_realarray[SoArealData::radius][lp];

          if (x_lo_bc && p.pos(0) < p_lo[0]+offset)
              p.pos(0) = p_lo[0] + offset;

          if (x_hi_bc && p.pos(0) + offset > p_hi[0])
            p.pos(0) = p_hi[0] - offset;

          if (y_lo_bc && p.pos(1) < p_lo[1]+offset)
            p.pos(1) = p_lo[1] + offset;

          if (y_hi_bc && p.pos(1) + offset > p_hi[1])
            p.pos(1) = p_hi[1] - offset;

          if (z_lo_bc && p.pos(2) < p_lo[2]+offset)
            p.pos(2) = p_lo[2] + offset;

          if (z_hi_bc && p.pos(2) + offset > p_hi[2])
            p.pos(2) = p_hi[2] - offset;

        });



      if(advect_enthalpy){

        amrex::ParallelFor(nrp, [p_realarray,dt]
            AMREX_GPU_DEVICE (int lp) noexcept
        {
          AMREX_ASSERT(p_realarray[SoArealData::c_ps][lp] > 0.);

          p_realarray[SoArealData::temperature][lp] += dt * p_realarray[SoArealData::convection][lp] /
            (p_realarray[SoArealData::mass][lp] * p_realarray[SoArealData::c_ps][lp]);
        });
      }


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


      Gpu::synchronize();

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR_STOP(pic_time_march);
#endif


    } // particle-tile iterator


  } // loop over levels


}
