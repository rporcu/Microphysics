#include <mfix.H>
#include <mfix_bc_parms.H>
#include <mfix_solids_parms.H>
#include <mfix_dem_parms.H>
#include <mfix_reactions_parms.H>
#include <mfix_calc_species_coeffs_K.H>

using namespace amrex;

void MFIXParticleContainer::MFIX_PC_AdvanceParcels (Real dt,
                                                    RealVect& gravity,
                                                    Vector< MultiFab* >& avg_prop_in,
                                                    Vector< MultiFab* >& cost,
                                                    std::string& knapsack_weight_type,
                                                    const int advect_enthalpy,
                                                    const Real enthalpy_source)
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

      auto ptile_data = ptile.getParticleTileData();

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

      const int nspecies_s = SOLIDS::nspecies;
      const int nreactions = REACTIONS::nreactions;

      // Particles SoA starting indexes for mass fractions and rate of
      // formations
      const int idx_X_sn       = m_runtimeRealData.X_sn;
      const int idx_ro_sn_txfr  = m_runtimeRealData.ro_sn_txfr;
      const int idx_vel_s_txfr = m_runtimeRealData.vel_s_txfr;
      const int idx_h_s_txfr   = m_runtimeRealData.h_s_txfr;

      const int update_mass           = SOLIDS::solve_species && REACTIONS::solve;
      const int update_temperature    = advect_enthalpy;
      const int local_advect_enthalpy = advect_enthalpy;
      const int local_solve_reactions = REACTIONS::solve;

      const int solid_is_mixture = SOLIDS::is_a_mixture;

      Gpu::AsyncArray<Real> d_cp_sn0_loc(SOLIDS::cp_sn0.dataPtr(), SOLIDS::cp_sn0.size());
      Real* p_cp_sn0_loc = d_cp_sn0_loc.data();

      const Real T_ref = SOLIDS::T_ref;

      const Real tolerance = std::numeric_limits<Real>::epsilon();

      amrex::ParallelFor(nrp,
        [pstruct,p_realarray,ptile_data,dt,gravity,p_hi,p_lo,dxi,tolerance,
         avg_prop_array,velfac,en,three_sqrt_two,x_lo_bc,x_hi_bc,y_lo_bc,y_hi_bc,
         z_lo_bc,z_hi_bc,nspecies_s,nreactions,idx_X_sn,idx_ro_sn_txfr,idx_vel_s_txfr,
         update_mass,update_temperature,local_solve_reactions,idx_h_s_txfr,T_ref,
         p_cp_sn0_loc,solid_is_mixture,local_advect_enthalpy,enthalpy_source]
        AMREX_GPU_DEVICE (int lp) noexcept
      {
        auto& p = pstruct[lp];

        GpuArray<Real,SPECIES::NMAX> X_sn;

        // Get current particle's mass
        Real p_mass_old = p_realarray[SoArealData::mass][lp];
        Real p_mass_new(p_mass_old);

        // Get current particle's density
        const Real p_density_old = p_realarray[SoArealData::density][lp];
        Real p_density_new(p_density_old);

        // Get current particle's volume
        Real p_vol = p_realarray[SoArealData::volume][lp];

        // Flag to stop computing particle's quantities if mass_new < 0, i.e.
        // the particle disappears because of chemical reactions
        int proceed = 1;

        //*********************************************************************
        // First step: update parcels' mass and density
        //*********************************************************************
        if (update_mass) {
          // Total particle density exchange rate
          Real total_ro_rate(0);

          // Loop over species
          for (int n_s(0); n_s < nspecies_s; n_s++)
          {
            // Current species mass fraction
            X_sn[n_s] = ptile_data.m_runtime_rdata[idx_X_sn+n_s][lp];

            // Get the current reaction rate for species n_s
            const Real ro_sn_rate = ptile_data.m_runtime_rdata[idx_ro_sn_txfr+n_s][lp];

            X_sn[n_s] = X_sn[n_s]*p_density_old + dt*ro_sn_rate;

            // Update the total mass exchange rate
            total_ro_rate += ro_sn_rate;
          }

          // Update the total mass of the particle
          p_density_new = p_density_old + dt * total_ro_rate;

          if (p_density_new > 0) {

            Real total_X(0.);

            // Normalize species mass fractions
            for (int n_s(0); n_s < nspecies_s; n_s++) {
              Real X_sn_new = X_sn[n_s] / p_density_new;

              if (X_sn_new < 0) X_sn_new = 0;
              if (X_sn_new > 1) X_sn_new = 1;

              total_X += X_sn_new;
              X_sn[n_s] = X_sn_new;
            }

            for (int n_s(0); n_s < nspecies_s; n_s++) {
              // Divide updated species mass fractions by total_X
              ptile_data.m_runtime_rdata[idx_X_sn+n_s][lp] = X_sn[n_s] / total_X;
            }

            // Write out to global memory particle's mass and density
            p_realarray[SoArealData::density][lp] = p_density_new;
            p_mass_new = p_density_new * p_vol;
            p_realarray[SoArealData::mass][lp] = p_mass_new;
          } else {
            p.id() = -1;
            proceed = 0;
          }
        }

        if (proceed) {
          //*********************************************************************
          // Second step: update parcels' positions and velocities
          //*********************************************************************

          // position
          RealVect pos = p.pos();

          // solids volume fraction
          const amrex::Real eps_p = amrex::max(1.0e-8, p_realarray[SoArealData::oneOverI][lp]);

          // inverse of particle mass
          const Real inv_mass = 1./p_mass_new;

          // scale wrt drag coefficient and step size
          const amrex::Real scale = 1.0 / (1.0 + dt*p_realarray[SoArealData::dragcoeff][lp]*inv_mass);

          // solids stress gradient:
          // grad_tau_p = (volume * grad_tau_p) / (mass * ep_s))
          amrex::RealVect grad_tau_p(0.);
          grad_tau_p[0] = p_realarray[SoArealData::omegax][lp];
          grad_tau_p[1] = p_realarray[SoArealData::omegay][lp];
          grad_tau_p[2] = p_realarray[SoArealData::omegaz][lp];

          const Real vel_coeff = update_mass ? p_mass_old/p_mass_new : 1.;

          // updated parcel velocity devoid of the particle normal stress
          const Real p_velx_old = p_realarray[SoArealData::velx][lp];
          const Real p_vely_old = p_realarray[SoArealData::vely][lp];
          const Real p_velz_old = p_realarray[SoArealData::velz][lp];

          RealVect vel(0.);
          vel[0] = scale*(vel_coeff*p_velx_old +
              dt*(p_realarray[SoArealData::dragx][lp]*inv_mass + vel_coeff*gravity[0]));
          vel[1] = scale*(vel_coeff*p_vely_old +
              dt*(p_realarray[SoArealData::dragy][lp]*inv_mass + vel_coeff*gravity[1]));
          vel[2] = scale*(vel_coeff*p_velz_old +
              dt*(p_realarray[SoArealData::dragz][lp]*inv_mass + vel_coeff*gravity[2]));

          if (local_solve_reactions) {
            const Real inv_density = 1. / p_density_new;
            vel[0] += dt*inv_density*scale*ptile_data.m_runtime_rdata[idx_vel_s_txfr+0][lp];
            vel[1] += dt*inv_density*scale*ptile_data.m_runtime_rdata[idx_vel_s_txfr+1][lp];
            vel[2] += dt*inv_density*scale*ptile_data.m_runtime_rdata[idx_vel_s_txfr+2][lp];
          }

          const Real lx = (pos[0] - p_lo[0])*dxi[0] + 0.5;
          const Real ly = (pos[1] - p_lo[1])*dxi[1] + 0.5;
          const Real lz = (pos[2] - p_lo[2])*dxi[2] + 0.5;

          const int i = static_cast<int>(Math::floor(lx));
          const int j = static_cast<int>(Math::floor(ly));
          const int k = static_cast<int>(Math::floor(lz));

          const Real wx_hi(lx - static_cast<Real>(i));
          const Real wy_hi(ly - static_cast<Real>(j));
          const Real wz_hi(lz - static_cast<Real>(k));

          const Real wx_lo(1 - wx_hi);
          const Real wy_lo(1 - wy_hi);
          const Real wz_lo(1 - wz_hi);

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

          const RealVect avg_vel = {u_face, v_face, w_face};

          const amrex::Real inv_rops = 1.0 / (eps_p*p_density_new);

          const amrex::Real vel_limit = p_realarray[SoArealData::radius][lp] / (dt * three_sqrt_two * eps_p);

          for (int dir(0); dir < 3; dir++)
          {
            Real bulk_vel;

            if( Math::abs(avg_vel[dir]) > tolerance &&
                Math::abs(    vel[dir]) > tolerance ) {
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

          // move the parcels
          pos[0] += dt * vel[0];
          pos[1] += dt * vel[1];
          pos[2] += dt * vel[2];

          // Impose domain constraints. Note that we only make sure that a parcel's
          // centroid is inside. A later routine will take into account the wall
          // collision and redirection of velocity.

          // Take the particle radius as the offset. This should be (much) smaller
          // than the effective parcel radius. We only want to move the parcel a
          // little but it needs to be enough to keep it inside the domain.
          const amrex::Real offset = p_realarray[SoArealData::radius][lp];

          if (x_lo_bc && pos[0] < p_lo[0]+offset)
              pos[0] = p_lo[0] + offset;

          if (x_hi_bc && pos[0] + offset > p_hi[0])
            pos[0] = p_hi[0] - offset;

          if (y_lo_bc && pos[1] < p_lo[1]+offset)
            pos[1] = p_lo[1] + offset;

          if (y_hi_bc && pos[1] + offset > p_hi[1])
            pos[1] = p_hi[1] - offset;

          if (z_lo_bc && pos[2] < p_lo[2]+offset)
            pos[2] = p_lo[2] + offset;

          if (z_hi_bc && pos[2] + offset > p_hi[2])
            pos[2] = p_hi[2] - offset;

          // Update positions
          p.pos(0) = pos[0];
          p.pos(1) = pos[1];
          p.pos(2) = pos[2];

          // update parcel velocity
          p_realarray[SoArealData::velx][lp] = vel[0];
          p_realarray[SoArealData::vely][lp] = vel[1];
          p_realarray[SoArealData::velz][lp] = vel[2];

          //*********************************************************************
          // Third step: update parcels' temperature
          //*********************************************************************
          if(local_advect_enthalpy) {
            const Real cp_s_old = p_realarray[SoArealData::cp_s][lp];
            Real cp_s_new(0);

            if (solid_is_mixture) {
              for (int n_s(0); n_s < nspecies_s; ++n_s)
                cp_s_new += p_cp_sn0_loc[n_s] * ptile_data.m_runtime_rdata[idx_X_sn+n_s][lp];

              p_realarray[SoArealData::cp_s][lp] = cp_s_new;
            } else {
              cp_s_new = cp_s_old;
            }

            AMREX_ASSERT(cp_s_new > 0.);

            if (! update_mass) {
              p_realarray[SoArealData::temperature][lp] +=
                dt*(p_realarray[SoArealData::convection][lp]+enthalpy_source) / (p_mass_new*cp_s_new);
            } else {
              Real p_enthalpy_new =
                p_mass_old*SOLIDS::calc_h_s(cp_s_old, p_realarray[SoArealData::temperature][i], 0, 0) +
                dt*(p_realarray[SoArealData::convection][i]+enthalpy_source);

              p_enthalpy_new -= dt*p_vol*ptile_data.m_runtime_rdata[idx_h_s_txfr][i];

              p_realarray[SoArealData::temperature][i] =
                SOLIDS::calc_T_s(cp_s_new, p_enthalpy_new, 0, 0) / p_mass_new;
            }
          }
        }
      });

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
