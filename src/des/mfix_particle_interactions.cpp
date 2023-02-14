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
ParticleInteractions::particle_wall (const int i,
                                     const int solve_enthalpy,
                                     const int ls_refinement,
                                     const Array4<const Real>& phiarr,
                                     const GpuArray<Real,AMREX_SPACEDIM>& plo,
                                     const GpuArray<Real,AMREX_SPACEDIM>& dxi,
                                     const Real small_number,
                                     const Real dt,
                                     RealVect& total_force,
                                     RealVect& total_tow_force,
                                     Real& Q_dot) const
{
  const ParticleType& particle = m_pstruct[i];

  const int istate = m_p_intarray[SoAintData::state][i];

  Real rp = m_p_realarray[SoArealData::radius][i];

  RealVect pos(particle.pos());

  Real ls_value = interp_level_set(pos, ls_refinement, phiarr, plo, dxi);

  Real overlap_n = rp - ls_value;

  // PFW conduction
  Real Tp1, Tp2;

  if(solve_enthalpy && m_solids_parms.get_do_pfp_cond<run_on>()) {

    const Real FLPC = m_solids_parms.get_flpc<run_on>();
    Real Rlens      = (1.0 + FLPC)*rp;

    if (ls_value < Rlens) {

      const Real Rough = m_solids_parms.get_min_cond<run_on>();
      Tp1              = m_p_realarray[SoArealData::temperature][i];

      // Construct a point inside the wall (to machine precision)
      RealVect normal(0.);
      level_set_normal(pos, ls_refinement, normal, phiarr, plo, dxi);
      normal[0] *= -1;
      normal[1] *= -1;
      normal[2] *= -1;
      RealVect posw = normal*(ls_value + small_number) + pos;

      // Find BC region this point lives in and get Twall
      Tp2 = Tp1;
      for (int bcv(0); bcv < m_bc_parms.bc_tw_count; ++bcv) {
        if (m_bc_parms.p_bc_rbv[bcv].contains(posw))
          Tp2 = m_bc_parms.p_bc_twv[bcv];
      }

      Q_dot += 2.*des_pfp_conduction(ls_value, rp, Rlens, Rough, m_dem_parms.k_g, Tp1, Tp2);
    }
  }

  if (ls_value < rp) {

    // PP conduction (Tw already found from PFW conduction hit)
    if(solve_enthalpy && m_solids_parms.get_do_pfp_cond<run_on>()) {

      const Real kp1 = m_solids_parms.calc_kp_sn<run_on>(Tp1,0);
      const Real kp2 = m_solids_parms.calc_kp_sn<run_on>(Tp2,0);

      Q_dot += des_pp_conduction(ls_value+rp, rp, rp, kp1, kp2, Tp1, Tp2);
    }

    RealVect normal(0.);
    level_set_normal(pos, ls_refinement, normal, phiarr, plo, dxi);

    normal[0] *= -1;
    normal[1] *= -1;
    normal[2] *= -1;

    RealVect v_rot(0.);
    v_rot[0] = ls_value * m_p_realarray[SoArealData::omegax][i];
    v_rot[1] = ls_value * m_p_realarray[SoArealData::omegay][i];
    v_rot[2] = ls_value * m_p_realarray[SoArealData::omegaz][i];

    RealVect vreltrans(0.);
    RealVect cprod(0.);

    cross_product(v_rot, normal, cprod);
    vreltrans[0] = m_p_realarray[SoArealData::velx][i] + cprod[0];
    vreltrans[1] = m_p_realarray[SoArealData::vely][i] + cprod[1];
    vreltrans[2] = m_p_realarray[SoArealData::velz][i] + cprod[2];

    Real vreltrans_norm = dot_product(vreltrans, normal);

    RealVect vrel_t(0.);
    vrel_t[0] = vreltrans[0] - vreltrans_norm*normal[0];
    vrel_t[1] = vreltrans[1] - vreltrans_norm*normal[1];
    vrel_t[2] = vreltrans[2] - vreltrans_norm*normal[2];

    const int phase = m_p_intarray[SoAintData::phase][i];
    const int phase_idx = MFIXSolidsPhase::phase_to_index(phase);

    Real kn_des_w   = m_dem_parms.kn_w;
    Real etan_des_w = m_dem_parms.etan_w(phase_idx);

    // NOTE - we don't use the tangential components right now,
    // but we might in the future
    // Real kt_des_w = m_dem.kt_w;
    // Real etat_des_w = m_dem.etat_w()[phase_idx];

    RealVect local_fn(0.);
    RealVect local_ft(0.);
    RealVect overlap_t(0.);
    Real mag_overlap_t(0.);

    // calculate the normal contact force
    local_fn[0] = -(kn_des_w*overlap_n*normal[0]
                  + etan_des_w*vreltrans_norm*normal[0]);
    local_fn[1] = -(kn_des_w*overlap_n*normal[1]
                  + etan_des_w*vreltrans_norm*normal[1]);
    local_fn[2] = -(kn_des_w*overlap_n*normal[2]
                  + etan_des_w*vreltrans_norm*normal[2]);

    // calculate the tangential displacement
    overlap_t[0] = dt*vrel_t[0];
    overlap_t[1] = dt*vrel_t[1];
    overlap_t[2] = dt*vrel_t[2];

    mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));

    if (mag_overlap_t > 0.0) {

      Real fnmd = m_dem_parms.mew_w * sqrt(dot_product(local_fn, local_fn));
      RealVect tangent(0.);
      tangent[0] = overlap_t[0]/mag_overlap_t;
      tangent[1] = overlap_t[1]/mag_overlap_t;
      tangent[2] = overlap_t[2]/mag_overlap_t;
      local_ft[0] = -fnmd * tangent[0];
      local_ft[1] = -fnmd * tangent[1];
      local_ft[2] = -fnmd * tangent[2];

    } else {

      local_ft[0] = 0.0;
      local_ft[1] = 0.0;
      local_ft[2] = 0.0;
    }

    if ( istate > 0 ) { // normal particles

      total_force[0] += local_fn[0] + local_ft[0];
      total_force[1] += local_fn[1] + local_ft[1];
      total_force[2] += local_fn[2] + local_ft[2];

      RealVect tow_force(0.);

      cross_product(normal, local_ft, tow_force);

      total_tow_force[0] += ls_value*tow_force[0];
      total_tow_force[1] += ls_value*tow_force[1];
      total_tow_force[2] += ls_value*tow_force[2];

    } else { // entering particles

      Real velx = m_p_realarray[SoArealData::velx][i];
      Real vely = m_p_realarray[SoArealData::vely][i];
      Real velz = m_p_realarray[SoArealData::velz][i];

      Real velmag = std::sqrt(velx*velx + vely*vely + velz*velz);

      Real dotprod = (normal[0] * velx +
                      normal[1] * vely +
                      normal[2] * velz)/velmag;

      // This is to catch particles that are not moving normal to
      // the levelset so that we can adjust their velocity and make sure
      // they fully enter the domain.
      if(Math::abs(1.0 + dotprod) > std::numeric_limits<Real>::epsilon()) {

        m_p_realarray[SoArealData::velx][i] = -velmag*normal[0];
        m_p_realarray[SoArealData::vely][i] = -velmag*normal[1];
        m_p_realarray[SoArealData::velz][i] = -velmag*normal[2];

      }
    }

  // An entering particle is no longer overlapping the wall.
  } else if(istate == 0) {
    //amrex::AllPrint() << "setting particle to normal\n";

    // Set the state to normal so it no longer ignores forces.
    m_p_intarray[SoAintData::state][i] = 1;
  }
}


AMREX_GPU_HOST_DEVICE
void
MFIXParticleContainer::
ParticleInteractions::particle_particle (const int i,
                                         Neighbors<ParticleType>::const_iterator mit,
                                         const int solve_enthalpy,
                                         const Real small_number,
                                         const amrex::Real dt,
                                         const int nrp,
                                         int& has_collisions,
                                         RealVect& total_force_i,
                                         RealVect& total_tow_force_i,
                                         Real& Q_dot_i,
                                         RealVect& total_force_j,
                                         RealVect& total_tow_force_j,
                                         Real& Q_dot_j) const
{
  const ParticleType& particle = m_pstruct[i];

  const int istate = m_p_intarray[SoAintData::state][i];

  RealVect pos1(particle.pos());

  const auto p2 = *mit;
  const int j = mit.index();

  Real dist_x = p2.pos(0) - pos1[0];
  Real dist_y = p2.pos(1) - pos1[1];
  Real dist_z = p2.pos(2) - pos1[2];

  Real r2 = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

  const Real p1radius = m_p_realarray[SoArealData::radius][i];
  const Real p2radius = m_p_realarray[SoArealData::radius][j];

  Real r_lm = p1radius + p2radius;

  AMREX_ASSERT_WITH_MESSAGE(
      !(particle.id() == p2.id() &&
        particle.cpu() == p2.cpu()),
    "A particle should not be its own neighbor!");

  // PFP conduction
  if(solve_enthalpy && m_solids_parms.get_do_pfp_cond<run_on>()) {

    const Real FLPC = m_solids_parms.get_flpc<run_on>();
    Real Rp_eff     = 2.0*(p1radius*p2radius)/(p1radius + p2radius);
    Real Rlens_eff  = (1.0 + FLPC)*Rp_eff;
    Real lens_lm    = 2.0*Rlens_eff;

    if ( r2 <= (lens_lm - small_number)*(lens_lm - small_number) ) {

      const Real Rough  = m_solids_parms.get_min_cond<run_on>();
      const Real Tp1    = m_p_realarray[SoArealData::temperature][i];
      const Real Tp2    = m_p_realarray[SoArealData::temperature][j];

      Real dist_mag_eff = sqrt(r2)/2.0; // Two particles with a midpoint wall

      Real local_Q_dot = des_pfp_conduction(dist_mag_eff, Rp_eff, Rlens_eff,
                                            Rough, m_dem_parms.k_g, Tp1, Tp2);

      Q_dot_i += local_Q_dot;
      Q_dot_j -= local_Q_dot;
    }
  }

  if ( r2 <= (r_lm - small_number)*(r_lm - small_number) )
  {
    has_collisions = 1;

    const int jstate = m_p_intarray[SoAintData::state][j];

    Real dist_mag = sqrt(r2);

    AMREX_ASSERT(dist_mag >= std::numeric_limits<Real>::epsilon());

    // PP conduction
    if(solve_enthalpy && m_solids_parms.get_do_pfp_cond<run_on>()) {

      const Real Tp1 = m_p_realarray[SoArealData::temperature][i];
      const Real Tp2 = m_p_realarray[SoArealData::temperature][j];
      const Real kp1 = m_solids_parms.calc_kp_sn<run_on>(Tp1,0);
      const Real kp2 = m_solids_parms.calc_kp_sn<run_on>(Tp2,0);

      Real local_Q_dot = des_pp_conduction(dist_mag, p1radius, p2radius, kp1,
                                           kp2, Tp1, Tp2);

      Q_dot_i += local_Q_dot;
      Q_dot_j -= local_Q_dot;
    }

    Real dist_mag_inv = 1.e0/dist_mag;

    RealVect normal(0.);
    normal[0] = dist_x * dist_mag_inv;
    normal[1] = dist_y * dist_mag_inv;
    normal[2] = dist_z * dist_mag_inv;

    Real overlap_n(0.);

    if (istate == 10 || jstate == 10) {

      // most of overlaps (99.99%) are in the range [0, 2.5e-8] m
      // which means [0, 5.e-4] radiuses
      // we set max overlap to   2.5e-4*radius
      overlap_n = amrex::min(r_lm - dist_mag, 2.5e-4*p1radius);

    } else {
      overlap_n = r_lm - dist_mag;
    }

    Real vrel_trans_norm;
    RealVect vrel_t(0.);

    RealVect p1vel(m_p_realarray[SoArealData::velx][i],
                   m_p_realarray[SoArealData::vely][i],
                   m_p_realarray[SoArealData::velz][i]);

    RealVect p2vel(m_p_realarray[SoArealData::velx][j],
                   m_p_realarray[SoArealData::vely][j],
                   m_p_realarray[SoArealData::velz][j]);

    RealVect p1omega(m_p_realarray[SoArealData::omegax][i],
                     m_p_realarray[SoArealData::omegay][i],
                     m_p_realarray[SoArealData::omegaz][i]);

    RealVect p2omega(m_p_realarray[SoArealData::omegax][j],
                     m_p_realarray[SoArealData::omegay][j],
                     m_p_realarray[SoArealData::omegaz][j]);

    cfrelvel(p1vel, p2vel, p1radius, p2radius, p1omega, p2omega,
             vrel_trans_norm, vrel_t, normal, dist_mag);

    const int phase1 = m_p_intarray[SoAintData::phase][i];
    const int phase2 = m_p_intarray[SoAintData::phase][j];

    const int phase1_idx = MFIXSolidsPhase::phase_to_index(phase1);
    const int phase2_idx = MFIXSolidsPhase::phase_to_index(phase2);

    Real kn_des   = m_dem_parms.kn;
    Real etan_des = m_dem_parms.etan(phase1_idx, phase2_idx);

    // NOTE - we don't use the tangential components right now,
    // but we might in the future
    // Real kt_des = m_dem.kt;
    // Real etat_des = m_dem.etat[phase1_idx][phase2_idx];

    RealVect local_fn(0.);
    RealVect local_ft(0.);
    RealVect overlap_t(0.);
    Real mag_overlap_t(0.);

    // calculate the normal contact force
    local_fn[0] = -(kn_des*overlap_n*normal[0]
                  + etan_des*vrel_trans_norm*normal[0]);
    local_fn[1] = -(kn_des*overlap_n*normal[1]
                  + etan_des*vrel_trans_norm*normal[1]);
    local_fn[2] = -(kn_des*overlap_n*normal[2]
                  + etan_des*vrel_trans_norm*normal[2]);

    // calculate the tangential overlap
    overlap_t[0] = dt*vrel_t[0];
    overlap_t[1] = dt*vrel_t[1];
    overlap_t[2] = dt*vrel_t[2];
    mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));

    if (mag_overlap_t > 0.0) {
      Real fnmd = m_dem_parms.mew * sqrt(dot_product(local_fn, local_fn));
      RealVect tangent(0.);
      tangent[0] = overlap_t[0]/mag_overlap_t;
      tangent[1] = overlap_t[1]/mag_overlap_t;
      tangent[2] = overlap_t[2]/mag_overlap_t;
      local_ft[0] = -fnmd * tangent[0];
      local_ft[1] = -fnmd * tangent[1];
      local_ft[2] = -fnmd * tangent[2];

    } else {
      local_ft[0] = 0.0;
      local_ft[1] = 0.0;
      local_ft[2] = 0.0;
    }

    Real dist_cl1 = 0.5 * (dist_mag + (p1radius*p1radius - p2radius*p2radius) * dist_mag_inv);
    dist_cl1 = dist_mag - dist_cl1;

    Real dist_cl2 = 0.5 * (dist_mag + (p2radius*p2radius - p1radius*p1radius) * dist_mag_inv);
    dist_cl2 = dist_mag - dist_cl2;

    RealVect local_tow_force(0.);
    cross_product(normal, local_ft, local_tow_force);

    if ( istate > 0 ) {
      total_force_i[0] += local_fn[0] + local_ft[0];
      total_force_i[1] += local_fn[1] + local_ft[1];
      total_force_i[2] += local_fn[2] + local_ft[2];

      total_tow_force_i[0] += dist_cl1*local_tow_force[0];
      total_tow_force_i[1] += dist_cl1*local_tow_force[1];
      total_tow_force_i[2] += dist_cl1*local_tow_force[2];
    }

    if (j < nrp && jstate != 0) {
      total_force_j[0] -= local_fn[0] + local_ft[0];
      total_force_j[1] -= local_fn[1] + local_ft[1];
      total_force_j[2] -= local_fn[2] + local_ft[2];

      total_tow_force_j[0] += dist_cl2*local_tow_force[0];
      total_tow_force_j[1] += dist_cl2*local_tow_force[1];
      total_tow_force_j[2] += dist_cl2*local_tow_force[2];
    }

    // Special case of two entering particles having an overlap
    if (istate == 0 && jstate == 0) {

      const Real shift = 1.0001*overlap_n;
      const RealVect sumvel(p1vel + p2vel);
      const int imove = (( sumvel[0]*normal[0]
                         + sumvel[1]*normal[1]
                         + sumvel[2]*normal[2]) > 0.) ? 1 : 0;

      if (imove) {
        total_force_i[0] -= shift * normal[0];
        total_force_i[1] -= shift * normal[1];
        total_force_i[2] -= shift * normal[2];

      } else if (j < nrp) {
        {
          total_force_j[0] += shift * normal[0];
          total_force_j[1] += shift * normal[1];
          total_force_j[2] += shift * normal[2];
        }
      }
    } // end overlap between entering particles

  } // end overlap
}
