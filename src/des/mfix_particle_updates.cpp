#include <mfix_des_K.H>

#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_bc.H>
#include <mfix_solvers.H>
#include <mfix_monitors.H>
#include <mfix_calc_cell.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE
void
MFIXParticleContainer::
ParticleUpdates::update_mass (const int i,
                              const int nspecies_s,
                              GpuArray<Real, MFIXSpecies::NMAX>& X_sn,
                              const Real dt,
                              const SoAruntimerealData& runtimedata_idxs,
                              int& proceed,
                              const int solve_dem) const
{
  ParticleType& particle = m_pstruct[i];

  // Get current particle's mass
  const Real p_mass_old = m_p_realarray[SoArealData::mass][i];
  Real p_mass_new(p_mass_old);

  // Get current particle's density
  const Real p_density_old = m_p_realarray[SoArealData::density][i];
  Real p_density_new(p_density_old);

  // Get current particle's oneOverI
  const Real p_oneOverI_old = m_p_realarray[SoArealData::oneOverI][i];
  Real p_oneOverI_new(p_oneOverI_old);

  // Get current particle's volume
  const Real p_vol = m_p_realarray[SoArealData::volume][i];

  // Total particle density exchange rate
  Real total_mass_rate(0);

  for (int n_s(0); n_s < nspecies_s; ++n_s) {

    // Get the current reaction rate for species n_s
    const int idx = runtimedata_idxs.chem_species_mass_txfr;
    const Real mass_sn_rate = m_ptile_data.m_runtime_rdata[idx+n_s][i];

    X_sn[n_s] = X_sn[n_s]*p_mass_old + dt*mass_sn_rate;

    // Update the total mass exchange rate
    total_mass_rate += mass_sn_rate;
  }

  // Update the total mass of the particle
  p_mass_new = p_mass_old + dt * total_mass_rate;

  if (p_mass_new > 0) {

    Real total_X(0.);

    // Normalize species mass fractions
    for (int n_s(0); n_s < nspecies_s; n_s++) {
      Real X_sn_new = X_sn[n_s] / p_mass_new;

      if (X_sn_new < 0) X_sn_new = 0;
      if (X_sn_new > 1) X_sn_new = 1;

      total_X += X_sn_new;
      X_sn[n_s] = X_sn_new;
    }

    for (int n_s(0); n_s < nspecies_s; n_s++) {
      // Divide updated species mass fractions by total_X
      const int idx = runtimedata_idxs.species_mass_fractions;
      X_sn[n_s] /= total_X;
      m_ptile_data.m_runtime_rdata[idx+n_s][i] = X_sn[n_s];
    }

    // Write out to global memory particle's mass and density
    m_p_realarray[SoArealData::mass][i] = p_mass_new;
    p_density_new = p_mass_new / p_vol;
    m_p_realarray[SoArealData::density][i] = p_density_new;

    if (solve_dem) {
      // Write out to global memory particle's moment of inertia
      p_oneOverI_new = (p_density_old/p_density_new)*p_oneOverI_old;
      m_p_realarray[SoArealData::oneOverI][i] = p_oneOverI_new;
    }

  } else {
    particle.id() = -1;
    proceed = 0;
  }
}


AMREX_GPU_HOST_DEVICE
void
MFIXParticleContainer::
ParticleUpdates::update_momentum (const int i,
                                  const int solve_reactions,
                                  const Real mass_coeff,
                                  Real* fc_ptr,
                                  Real* tow_ptr,
                                  const int ntot,
                                  const Real dt,
                                  const SoAruntimerealData& runtimedata_idxs,
                                  const RealVect& gravity,
                                  const GpuArray<Real,AMREX_SPACEDIM>& p_lo,
                                  const GpuArray<Real,AMREX_SPACEDIM>& p_hi) const
{
  ParticleType& particle = m_pstruct[i];

  const Real p_mass_new = m_p_realarray[SoArealData::mass][i];
  const Real p_oneOverI_new = m_p_realarray[SoArealData::oneOverI][i];

  const Real p_velx_old = m_p_realarray[SoArealData::velx][i];
  const Real p_vely_old = m_p_realarray[SoArealData::vely][i];
  const Real p_velz_old = m_p_realarray[SoArealData::velz][i];

  Real p_velx_new = mass_coeff*p_velx_old +
    dt*((m_p_realarray[SoArealData::dragx][i]+fc_ptr[i]) / p_mass_new + mass_coeff*gravity[0]);
  Real p_vely_new = mass_coeff*p_vely_old +
    dt*((m_p_realarray[SoArealData::dragy][i]+fc_ptr[i+ntot]) / p_mass_new + mass_coeff*gravity[1]);
  Real p_velz_new = mass_coeff*p_velz_old +
    dt*((m_p_realarray[SoArealData::dragz][i]+fc_ptr[i+2*ntot]) / p_mass_new + mass_coeff*gravity[2]);

  if (solve_reactions) {
    const int idx = runtimedata_idxs.chem_vel_txfr;
    p_velx_new += dt*(m_ptile_data.m_runtime_rdata[idx+0][i] / p_mass_new);
    p_vely_new += dt*(m_ptile_data.m_runtime_rdata[idx+1][i] / p_mass_new);
    p_velz_new += dt*(m_ptile_data.m_runtime_rdata[idx+2][i] / p_mass_new);
  }

  const Real p_omegax_old = m_p_realarray[SoArealData::omegax][i];
  const Real p_omegay_old = m_p_realarray[SoArealData::omegay][i];
  const Real p_omegaz_old = m_p_realarray[SoArealData::omegaz][i];

  Real p_omegax_new = mass_coeff*p_omegax_old + dt * p_oneOverI_new * tow_ptr[i];
  Real p_omegay_new = mass_coeff*p_omegay_old + dt * p_oneOverI_new * tow_ptr[i+ntot];
  Real p_omegaz_new = mass_coeff*p_omegaz_old + dt * p_oneOverI_new * tow_ptr[i+2*ntot];

  const Real p_posx_old = particle.pos(0);
  const Real p_posy_old = particle.pos(1);
  const Real p_posz_old = particle.pos(2);

  Real p_posx_new = p_posx_old + dt * p_velx_new;
  Real p_posy_new = p_posy_old + dt * p_vely_new;
  Real p_posz_new = p_posz_old + dt * p_velz_new;

  const Real eps = std::numeric_limits<Real>::epsilon();

  if (m_bc_parms.domain_bc[0] && p_posx_new < p_lo[0])
  {
      p_posx_new = p_lo[0] + eps;
      p_velx_new = -p_velx_new;
  }
  if (m_bc_parms.domain_bc[1] && p_posx_new > p_hi[0])
  {
      p_posx_new = p_hi[0] - eps;
      p_velx_new = -p_velx_new;
  }
  if (m_bc_parms.domain_bc[2] && p_posy_new < p_lo[1])
  {
      p_posy_new = p_lo[1] + eps;
      p_vely_new = -p_vely_new;
  }
  if (m_bc_parms.domain_bc[3] && p_posy_new > p_hi[1])
  {
      p_posy_new = p_hi[1] - eps;
      p_vely_new = -p_vely_new;
  }
  if (m_bc_parms.domain_bc[4] && p_posz_new < p_lo[2])
  {
      p_posz_new = p_lo[2] + eps;
      p_velz_new = -p_velz_new;
  }
  if (m_bc_parms.domain_bc[5] && p_posz_new > p_hi[2])
  {
      p_posz_new = p_hi[2] - eps;
      p_velz_new = -p_velz_new;
  }

  if (m_p_intarray[SoAintData::state][i] != 0) {

    // Update positions
    particle.pos(0) = p_posx_new;
    particle.pos(1) = p_posy_new;
    particle.pos(2) = p_posz_new;

    // Update velocities
    m_p_realarray[SoArealData::velx][i] = p_velx_new;
    m_p_realarray[SoArealData::vely][i] = p_vely_new;
    m_p_realarray[SoArealData::velz][i] = p_velz_new;

    // Update angular velocities
    m_p_realarray[SoArealData::omegax][i] = p_omegax_new;
    m_p_realarray[SoArealData::omegay][i] = p_omegay_new;
    m_p_realarray[SoArealData::omegaz][i] = p_omegaz_new;

  } else {

    particle.pos(0) += dt * p_velx_old + fc_ptr[i         ];
    particle.pos(1) += dt * p_vely_old + fc_ptr[i +   ntot];
    particle.pos(2) += dt * p_velz_old + fc_ptr[i + 2*ntot];
  }
}


template<class S>
AMREX_GPU_HOST_DEVICE
void
MFIXParticleContainer::
ParticleUpdates::update_enthalpy (const int i,
                                  const int nspecies_s,
                                  GpuArray<Real, MFIXSpecies::NMAX>& X_sn,
                                  const int solve_reactions,
                                  const int solid_is_a_mixture,
                                  const Real mass_coeff,
                                  Real* cond_ptr,
                                  const Real dt,
                                  const SoAruntimerealData& runtimedata_idxs,
                                  const Real enthalpy_source,
                                  S nonlinear_solver) const
{
  const Real p_mass_new = m_p_realarray[SoArealData::mass][i];

  Real p_enthalpy_old(0);

  const Real Tp = m_p_realarray[SoArealData::temperature][i];

  if (solid_is_a_mixture) {
    for (int n_s(0); n_s < nspecies_s; ++n_s) {
      p_enthalpy_old += X_sn[n_s]*m_solids_parms.calc_h_sn<run_on>(Tp, n_s);
    }
  } else {
    p_enthalpy_old = m_solids_parms.calc_h_s<run_on>(Tp);
  }

  Real p_enthalpy_new(0.);

  if (cond_ptr != nullptr) {
    p_enthalpy_new = mass_coeff*p_enthalpy_old +
      dt*((m_p_realarray[SoArealData::convection][i]+cond_ptr[i]+enthalpy_source) / p_mass_new);
  } else {
    p_enthalpy_new = mass_coeff*p_enthalpy_old +
      dt*((m_p_realarray[SoArealData::convection][i]+enthalpy_source) / p_mass_new);
  }

  if (solve_reactions) {
    const int idx = runtimedata_idxs.chem_enthalpy_txfr;
    p_enthalpy_new += dt*(m_ptile_data.m_runtime_rdata[idx][i] / p_mass_new);
  }

  // ************************************************************
  // Newton-Raphson solver for solving implicit equation for
  // temperature
  // ************************************************************
  MFIXSolidsParms solids_parms(m_solids_parms);
  
  // Residual computation
  auto R = [=] AMREX_GPU_HOST_DEVICE (Real Tp_arg)
  {
    Real hp_loc(0);

    if (!solid_is_a_mixture) {

      hp_loc = solids_parms.calc_h_s<run_on>(Tp_arg);
    } else {

      for (int n_s(0); n_s < nspecies_s; ++n_s)
        hp_loc += X_sn[n_s]*solids_parms.calc_h_sn<run_on>(Tp_arg,n_s);
    }

    return hp_loc - p_enthalpy_new;
  };

  // Partial derivative computation
  auto partial_R = [=] AMREX_GPU_HOST_DEVICE (Real Tp_arg)
  {
    Real gradient(0);

    if (!solid_is_a_mixture) {

      gradient = solids_parms.calc_partial_h_s<run_on>(Tp_arg);
    } else {

      for (int n_s(0); n_s < nspecies_s; ++n_s) {
        gradient += X_sn[n_s]*solids_parms.calc_partial_h_sn<run_on>(Tp_arg,n_s);
      }
    }

    return gradient;
  };

  //const Real Tp_old = m_p_realarray[SoArealData::temperature][i];
  //Real Tp_new(Tp_old);
  Real Tp_new(m_p_realarray[SoArealData::temperature][i]);

  nonlinear_solver.solve(Tp_new, R, partial_R);

  m_p_realarray[SoArealData::temperature][i] = Tp_new;

  // Update cp_s
  Real cp_s_new(0);

  if (solid_is_a_mixture) {
    for (int n_s(0); n_s < nspecies_s; ++n_s)
      cp_s_new += X_sn[n_s]*m_solids_parms.calc_cp_sn<run_on>(Tp_new,n_s);

  } else {
    cp_s_new = m_solids_parms.calc_cp_s<run_on>(Tp_new);
  }

  AMREX_ASSERT(cp_s_new > 0.);
  m_p_realarray[SoArealData::cp_s][i] = cp_s_new;
}


template
AMREX_GPU_HOST_DEVICE
void
MFIXParticleContainer::
ParticleUpdates::update_enthalpy<MFIXSolvers::Newton> (const int i,
                                                       const int nspecies_s,
                                                       GpuArray<Real, MFIXSpecies::NMAX>& X_sn,
                                                       const int solve_reactions,
                                                       const int solid_is_a_mixture,
                                                       const Real mass_coeff,
                                                       Real* cond_ptr,
                                                       const Real dt,
                                                       const SoAruntimerealData& runtimedata_idxs,
                                                       const Real enthalpy_source,
                                                       MFIXSolvers::Newton nonlinear_solver) const;
