#ifndef _MFIX_PC_UPDATES_K_H_
#define _MFIX_PC_UPDATES_K_H_

#include <mfix_des_K.H>
#include <mfix_pc.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_bc.H>
#include <mfix_solvers.H>
#include <mfix_monitors.H>
#include <mfix_calc_cell.H>


//***************************************************************
// First step: update particles' mass and density
//***************************************************************
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
part_mass_update (MFIXParticleContainer::ParticleType& p,
                  const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& ptile_data,
                  const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& p_realarray,
                  const int& i,
                  const int idx_X_sn,
                  const int& nspecies_s,
                  const amrex::Real& dt,
                  const amrex::Real& p_mass_old,
                  amrex::Real& p_mass_new,
                  const amrex::Real* p_oneOverI_old,
                  amrex::Real* p_oneOverI_new,
                  int& proceed,
                  amrex::Real& coeff,
                  const int& idx_mass_txfr,
                  const int& solve_dem)
{
  // Get current particle's density
  const amrex::Real p_density_old = p_realarray[SoArealData::density][i];
  amrex::Real p_density_new(p_density_old);

  // Get current particle's volume
  const amrex::Real p_vol = p_realarray[SoArealData::volume][i];

  // Total particle density exchange rate
  amrex::Real total_mass_rate(0);

  for (int n_s(0); n_s < nspecies_s; ++n_s) {

    // Get the current reaction rate for species n_s
    const amrex::Real mass_sn_rate = ptile_data.m_runtime_rdata[idx_mass_txfr+n_s][i];

    ptile_data.m_runtime_rdata[idx_X_sn+n_s][i] =
      ptile_data.m_runtime_rdata[idx_X_sn+n_s][i]*p_mass_old + dt*mass_sn_rate;

    // Update the total mass exchange rate
    total_mass_rate += mass_sn_rate;
  }

  // Update the total mass of the particle
  p_mass_new = p_mass_old + dt * total_mass_rate;

  if (p_mass_new > 0) {

    amrex::Real total_X(0.);

    // Normalize species mass fractions
    for (int n_s(0); n_s < nspecies_s; n_s++) {
      amrex::Real X_sn_new = ptile_data.m_runtime_rdata[idx_X_sn+n_s][i] / p_mass_new;

      if (X_sn_new < 0) X_sn_new = 0;
      if (X_sn_new > 1) X_sn_new = 1;

      total_X += X_sn_new;
      ptile_data.m_runtime_rdata[idx_X_sn+n_s][i] = X_sn_new;
    }

    for (int n_s(0); n_s < nspecies_s; n_s++) {
      // Divide updated species mass fractions by total_X
      ptile_data.m_runtime_rdata[idx_X_sn+n_s][i] /= total_X;
    }

    // Write out to global memory particle's mass and density
    p_realarray[SoArealData::mass][i] = p_mass_new;
    p_density_new = p_mass_new / p_vol;
    p_realarray[SoArealData::density][i] = p_density_new;

    coeff = p_mass_old / p_mass_new;

    // Write out to global memory particle's moment of inertia
    if (solve_dem) {
      *p_oneOverI_new = coeff * (*p_oneOverI_old);
      p_realarray[SoArealData::oneOverI][i] = *p_oneOverI_new;
    }

  } else {
    p.id() = -1;
    proceed = 0;
  }
}


//***************************************************************
// Second step: update particles' positions and velocities
//***************************************************************
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
part_momentum_update (MFIXParticleContainer::ParticleType& p,
                      const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& ptile_data,
                      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& p_realarray,
                      const amrex::GpuArray<int*, SoAintData::count>& p_intarray,
                      const int& i,
                      const amrex::Real& dt,
                      const amrex::Real& coeff,
                      const amrex::Real& p_mass_new,
                      const amrex::Real& p_oneOverI_new,
                      const int& ntot,
                      amrex::Real* fc_ptr,
                      amrex::Real* tow_ptr,
                      const amrex::RealVect& gravity,
                      const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& p_lo,
                      const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& p_hi,
                      const int& solve_reactions,
                      const int& idx_vel_txfr,
                      const amrex::GpuArray<int,6>& lo_hi_bc)
{
  const amrex::Real p_velx_old = p_realarray[SoArealData::velx][i];
  const amrex::Real p_vely_old = p_realarray[SoArealData::vely][i];
  const amrex::Real p_velz_old = p_realarray[SoArealData::velz][i];

  amrex::Real p_velx_new = coeff*p_velx_old +
    dt*((p_realarray[SoArealData::dragx][i]+fc_ptr[i]) / p_mass_new + coeff*gravity[0]);
  amrex::Real p_vely_new = coeff*p_vely_old +
    dt*((p_realarray[SoArealData::dragy][i]+fc_ptr[i+ntot]) / p_mass_new + coeff*gravity[1]);
  amrex::Real p_velz_new = coeff*p_velz_old +
    dt*((p_realarray[SoArealData::dragz][i]+fc_ptr[i+2*ntot]) / p_mass_new + coeff*gravity[2]);

  if (solve_reactions) {
    p_velx_new += dt*(ptile_data.m_runtime_rdata[idx_vel_txfr+0][i] / p_mass_new);
    p_vely_new += dt*(ptile_data.m_runtime_rdata[idx_vel_txfr+1][i] / p_mass_new);
    p_velz_new += dt*(ptile_data.m_runtime_rdata[idx_vel_txfr+2][i] / p_mass_new);
  }

  const amrex::Real p_omegax_old = p_realarray[SoArealData::omegax][i];
  const amrex::Real p_omegay_old = p_realarray[SoArealData::omegay][i];
  const amrex::Real p_omegaz_old = p_realarray[SoArealData::omegaz][i];

  amrex::Real p_omegax_new = coeff*p_omegax_old + dt * p_oneOverI_new * tow_ptr[i];
  amrex::Real p_omegay_new = coeff*p_omegay_old + dt * p_oneOverI_new * tow_ptr[i+ntot];
  amrex::Real p_omegaz_new = coeff*p_omegaz_old + dt * p_oneOverI_new * tow_ptr[i+2*ntot];

  const amrex::Real p_posx_old = p.pos(0);
  const amrex::Real p_posy_old = p.pos(1);
  const amrex::Real p_posz_old = p.pos(2);

  amrex::Real p_posx_new = p_posx_old + dt * p_velx_new;
  amrex::Real p_posy_new = p_posy_old + dt * p_vely_new;
  amrex::Real p_posz_new = p_posz_old + dt * p_velz_new;

  if (lo_hi_bc[0] && p_posx_new < p_lo[0])
  {
      p_posx_new = p_lo[0] + std::numeric_limits<amrex::Real>::epsilon();
      p_velx_new = -p_velx_new;
  }
  if (lo_hi_bc[1] && p_posx_new > p_hi[0])
  {
      p_posx_new = p_hi[0] - std::numeric_limits<amrex::Real>::epsilon();
      p_velx_new = -p_velx_new;
  }
  if (lo_hi_bc[2] && p_posy_new < p_lo[1])
  {
      p_posy_new = p_lo[1] + std::numeric_limits<amrex::Real>::epsilon();
      p_vely_new = -p_vely_new;
  }
  if (lo_hi_bc[3] && p_posy_new > p_hi[1])
  {
      p_posy_new = p_hi[1] - std::numeric_limits<amrex::Real>::epsilon();
      p_vely_new = -p_vely_new;
  }
  if (lo_hi_bc[4] && p_posz_new < p_lo[2])
  {
      p_posz_new = p_lo[2] + std::numeric_limits<amrex::Real>::epsilon();
      p_velz_new = -p_velz_new;
  }
  if (lo_hi_bc[5] && p_posz_new > p_hi[2])
  {
      p_posz_new = p_hi[2] - std::numeric_limits<amrex::Real>::epsilon();
      p_velz_new = -p_velz_new;
  }

  if (p_intarray[SoAintData::state][i] != 0) {

    // Update positions
    p.pos(0) = p_posx_new;
    p.pos(1) = p_posy_new;
    p.pos(2) = p_posz_new;

    // Update velocities
    p_realarray[SoArealData::velx][i] = p_velx_new;
    p_realarray[SoArealData::vely][i] = p_vely_new;
    p_realarray[SoArealData::velz][i] = p_velz_new;

    // Update angular velocities
    p_realarray[SoArealData::omegax][i] = p_omegax_new;
    p_realarray[SoArealData::omegay][i] = p_omegay_new;
    p_realarray[SoArealData::omegaz][i] = p_omegaz_new;

  } else {

    p.pos(0) += dt * p_velx_old + fc_ptr[i         ];
    p.pos(1) += dt * p_vely_old + fc_ptr[i +   ntot];
    p.pos(2) += dt * p_velz_old + fc_ptr[i + 2*ntot];
  }
}


//***************************************************************
// Third step: update particles' temperature
//***************************************************************
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
part_enthalpy_update (const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& ptile_data,
                      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& p_realarray,
                      const int& i,
                      const int idx_X_sn,
                      const int& nspecies_s,
                      const int& solid_is_a_mixture,
                      const MFIXSolidsParms& solids_parms,
                      const amrex::Real& dt,
                      const amrex::Real& coeff,
                      const amrex::Real& p_enthalpy_old,
                      const amrex::Real& p_mass_new,
                      amrex::Real* cond_ptr,
                      const amrex::Real& enthalpy_source,
                      const int& solve_reactions,
                      const int& idx_h_txfr,
                      const amrex::Real& abstol,
                      const amrex::Real& reltol,
                      const int& maxiter,
                      const int& solve_dem)
{
  amrex::Real p_enthalpy_new(0.);

  if (solve_dem) {
    p_enthalpy_new = coeff*p_enthalpy_old +
      dt*((p_realarray[SoArealData::convection][i]+cond_ptr[i]+enthalpy_source) / p_mass_new);
  } else {
    p_enthalpy_new = coeff*p_enthalpy_old +
      dt*((p_realarray[SoArealData::convection][i]+enthalpy_source) / p_mass_new);
  }

  if (solve_reactions) {
    p_enthalpy_new += dt*(ptile_data.m_runtime_rdata[idx_h_txfr][i] / p_mass_new);
  }

  // ************************************************************
  // Solvers::Newton-Raphson solver for solving implicit equation for
  // temperature
  // ************************************************************
  Solvers::Newton::SolidsEnthalpy::Residue residue(solid_is_a_mixture,
      i, nspecies_s, solids_parms, ptile_data, idx_X_sn, p_enthalpy_new);

  Solvers::Newton::SolidsEnthalpy::Gradient gradient(solid_is_a_mixture,
      i, nspecies_s, solids_parms, ptile_data, idx_X_sn);

  amrex::Real Tp_new(p_realarray[SoArealData::temperature][i]);

  auto output = Solvers::Newton::solve(Tp_new, residue, gradient, abstol, reltol, maxiter);

  if (output.iterations == -1)
    amrex::Abort("Newton solver did not converge");

  p_realarray[SoArealData::temperature][i] = Tp_new;

  // Update cp_s
  amrex::Real cp_s_new(0);

  if (solid_is_a_mixture) {
    for (int n_s(0); n_s < nspecies_s; ++n_s)
      cp_s_new += ptile_data.m_runtime_rdata[idx_X_sn+n_s][i] *
        solids_parms.calc_cp_sn<run_on>(Tp_new,n_s);

  } else {
    cp_s_new = solids_parms.calc_cp_s<run_on>(Tp_new);
  }

  AMREX_ASSERT(cp_s_new > 0.);
  p_realarray[SoArealData::cp_s][i] = cp_s_new;
}

#endif
