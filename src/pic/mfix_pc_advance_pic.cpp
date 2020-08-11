#include <mfix.H>
#include <mfix_bc_parms.H>


void MFIXParticleContainer::MFIX_PC_AdvanceParcels (amrex::Real dt, amrex::RealVect& gravity,
                                                    amrex::Vector< amrex::MultiFab* >& avg_prop_in,
                                                    amrex::Vector< amrex::MultiFab* >& cost,
                                                    std::string& knapsack_weight_type,
                                                    const int advect_enthalpy)
{

  BL_PROFILE("MFIXParticleContainer::MFIX_PC_AdvanceParcels()");

  const amrex::Real small_number = std::numeric_limits<Real>::epsilon();

  const amrex::Real en = (PIC::damping_factor + 1.0);
  const amrex::Real velfac = PIC::velfac;

  const amrex::Real three_sqrt_two(3.0*std::sqrt(2.0));

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
      //amrex::Real wt = ParallelDescriptor::second(); UNUSED

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

      const auto& avg_prop_array = avg_prop_in[lev]->array(pti);

      const Real tolerance = std::numeric_limits<Real>::epsilon();

      amrex::ParallelFor(nrp,
        [pstruct,dt,gravity, small_number,p_hi,p_lo, dxi, tolerance,
         avg_prop_array, velfac, en, three_sqrt_two,
         x_lo_bc,x_hi_bc,y_lo_bc,y_hi_bc,z_lo_bc,z_hi_bc]
        AMREX_GPU_DEVICE (int lp) noexcept
        {
          ParticleType& p = pstruct[lp];

          // position
          const amrex::RealVect pos = p.pos();

          // solids volume fraction
          const amrex::Real eps_p = amrex::max(1.0e-8, p.rdata(realData::oneOverI));

          // inverse of particle mass
          const amrex::Real inv_mass = 1.0 / p.rdata(realData::mass);

          // scale wrt drag coefficient and step size
          const amrex::Real scale = 1.0 / (1.0 + dt*p.rdata(realData::dragcoeff)*inv_mass);

          // solids stress gradient:
          // grad_tau_p = (volume * grad_tau_p) / (mass * ep_s))
          amrex::RealVect grad_tau_p;
          grad_tau_p[0] = p.rdata(realData::omegax);
          grad_tau_p[1] = p.rdata(realData::omegay);
          grad_tau_p[2] = p.rdata(realData::omegaz);

          // updated parcel velocity devoid of the particle normal stress
          amrex::RealVect vel;
          vel[0] = scale*(p.rdata(realData::velx) + dt*(p.rdata(realData::dragx)*inv_mass + gravity[0]));
          vel[1] = scale*(p.rdata(realData::vely) + dt*(p.rdata(realData::dragy)*inv_mass + gravity[1]));
          vel[2] = scale*(p.rdata(realData::velz) + dt*(p.rdata(realData::dragz)*inv_mass + gravity[2]));

          const amrex::Real lx = (pos[0] - p_lo[0])*dxi[0] + 0.5;
          const amrex::Real ly = (pos[1] - p_lo[1])*dxi[1] + 0.5;
          const amrex::Real lz = (pos[2] - p_lo[2])*dxi[2] + 0.5;

          const int i = static_cast<int>(amrex::Math::floor(lx));
          const int j = static_cast<int>(amrex::Math::floor(ly));
          const int k = static_cast<int>(amrex::Math::floor(lz));

          const amrex::Real wx_hi(lx - static_cast<amrex::Real>(i));
          const amrex::Real wy_hi(ly - static_cast<amrex::Real>(j));
          const amrex::Real wz_hi(lz - static_cast<amrex::Real>(k));

          const amrex::Real wx_lo(1.0 - wx_hi);
          const amrex::Real wy_lo(1.0 - wy_hi);
          const amrex::Real wz_lo(1.0 - wz_hi);

          int ii = i;
          if(vel[0]*grad_tau_p[0] > 0.) {
            ii = (vel[0]>0.) ? i+1 : i-1;
          }

          const amrex::Real u_face =
            wy_lo*wz_lo*avg_prop_array(ii,j-1,k-1,0) +
            wy_hi*wz_lo*avg_prop_array(ii,j  ,k-1,0) +
            wy_lo*wz_hi*avg_prop_array(ii,j-1,k  ,0) +
            wy_hi*wz_hi*avg_prop_array(ii,j  ,k  ,0);


          int jj = j;
          if(vel[1]*grad_tau_p[1] > 0.0) {
            jj = ( vel[1] > 0.0 ) ? j+1 : j-1;
          }

          const amrex::Real v_face =
            wx_lo*wz_lo*avg_prop_array(i-1,jj,k-1,1) +
            wx_hi*wz_lo*avg_prop_array(i  ,jj,k-1,1) +
            wx_lo*wz_hi*avg_prop_array(i-1,jj,k  ,1) +
            wx_hi*wz_hi*avg_prop_array(i  ,jj,k  ,1);


          int kk = k;
          if(vel[2]*grad_tau_p[2] > 0.0) {
            kk = ( vel[2] > 0.0 ) ? k+1 : k-1;
          }

          const amrex::Real w_face =
            wx_lo*wy_lo*avg_prop_array(i-1,j-1,kk,2) +
            wx_hi*wy_lo*avg_prop_array(i  ,j-1,kk,2) +
            wx_lo*wy_hi*avg_prop_array(i-1,j  ,kk,2) +
            wx_hi*wy_hi*avg_prop_array(i  ,j  ,kk,2);

          const amrex::RealVect avg_vel  = {u_face, v_face, w_face};

          const amrex::Real inv_rops = 1.0 / (eps_p*p.rdata(realData::density));

          const amrex::Real vel_limit = p.rdata(realData::radius) / (dt * three_sqrt_two * eps_p);

          for( int dir(0); dir<3; dir++)
          {

            amrex::Real bulk_vel;

            if( amrex::Math::abs(avg_vel[dir]) > tolerance and
                amrex::Math::abs(    vel[dir]) > tolerance ){
              bulk_vel = 2./(1./avg_vel[dir] + 1./vel[dir]);
            } else {
              bulk_vel = 0.0;
            }

            // Slip velocity between the parcel and the bulk. This needs to be the
            // tentative new velocity, not the current velocity.
            const amrex::Real slip_vel = velfac*bulk_vel - vel[dir];

            // solids stress velocity contribution:
            // -dt*( grad_tau_p / (density * ep_s)) / (1 + dt*beta/mass)
            const amrex::Real del_up = -dt*scale*grad_tau_p[dir] * inv_rops;

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
          p.rdata(realData::velx) = vel[0];
          p.rdata(realData::vely) = vel[1];
          p.rdata(realData::velz) = vel[2];

          // move the parcels
          p.pos(0) += dt * p.rdata(realData::velx);
          p.pos(1) += dt * p.rdata(realData::vely);
          p.pos(2) += dt * p.rdata(realData::velz);


          // Impose domain constraints. Note that we only make sure that a parcel's
          // centroid is inside. A later routine will take into account the wall
          // collision and redirection of velocity.

          if (x_lo_bc && p.pos(0) < p_lo[0]+small_number)
              p.pos(0) = p_lo[0] + small_number;

          if (x_hi_bc && p.pos(0) + small_number > p_hi[0])
            p.pos(0) = p_hi[0] - small_number;

          if (y_lo_bc && p.pos(1) < p_lo[1]+small_number)
            p.pos(1) = p_lo[1] + small_number;

          if (y_hi_bc && p.pos(1) + small_number > p_hi[1])
            p.pos(1) = p_hi[1] - small_number;

          if (z_lo_bc && p.pos(2) < p_lo[2]+small_number)
            p.pos(2) = p_lo[2] + small_number;

          if (z_hi_bc && p.pos(2) + small_number > p_hi[2])
            p.pos(2) = p_hi[2] - small_number;

        });



      if(advect_enthalpy){

        amrex::ParallelFor(nrp, [pstruct,dt] AMREX_GPU_DEVICE (int lp) noexcept
        {
          ParticleType& p = pstruct[lp];

          AMREX_ASSERT(p.rdata(realData::c_ps) > 0.);

          p.rdata(realData::temperature) += dt * p.rdata(realData::convection) /
            (p.rdata(realData::mass) * p.rdata(realData::c_ps));

        });
      }

      Gpu::synchronize();

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR_STOP(pic_time_march);
#endif


    } // particle-tile iterator


  } // loop over levels


}
