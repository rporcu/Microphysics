#include <AMReX.H>
#include "AMReX_Particles.H"
#include <mfix_pc.H>

#include <mfix_deposition_K.H>
#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_pic_K.H>

#include <mfix_bc.H>
#include <mfix_pic.H>

using namespace amrex;

void MFIXParticleContainer::
MFIX_PC_SolidsVelocityDeposition (int lev,
                                  Array<MultiFab*,3>& vel_s_in,
                                  const FabArray<EBCellFlagFab>* flags)
{
  BL_PROFILE("(MFIXParticleContainer::MFIX_PC_SolidsVelocityDeposition)");

#define USE_HAT_FUNCTION 0

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto p_lo = gm.ProbLoArray();
  const auto  dxi = gm.InvCellSizeArray();

  const auto inv_reg_cell_vol = dxi[0]*dxi[1]*dxi[2];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_u_s_fab;
    FArrayBox local_v_s_fab;
    FArrayBox local_w_s_fab;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();

      FArrayBox& u_s_fab = (*vel_s_in[0])[pti];
      FArrayBox& v_s_fab = (*vel_s_in[1])[pti];
      FArrayBox& w_s_fab = (*vel_s_in[2])[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered ) {

        auto u_s = u_s_fab.array();
        auto v_s = v_s_fab.array();
        auto w_s = w_s_fab.array();

#ifdef _OPENMP
        const int ncomp = vel_s_in[0]->nComp();
        AMREX_ALWAYS_ASSERT(ncomp == vel_s_in[1]->nComp());
        AMREX_ALWAYS_ASSERT(ncomp == vel_s_in[2]->nComp());

        const int ngrow = vel_s_in[0]->nGrow();
        AMREX_ASSERT(ngrow == vel_s_in[1]->nGrow());
        AMREX_ASSERT(ngrow == vel_s_in[1]->nGrow());

        Box const xbx = Box(box).grow(ngrow).surroundingNodes(0);
        Box const ybx = Box(box).grow(ngrow).surroundingNodes(1);
        Box const zbx = Box(box).grow(ngrow).surroundingNodes(2);

        if (Gpu::notInLaunchRegion())
        {
          local_u_s_fab.resize(xbx, ncomp);
          local_u_s_fab.setVal<RunOn::Host>(0.0);
          u_s = local_u_s_fab.array();

          local_v_s_fab.resize(ybx, ncomp);
          local_v_s_fab.setVal<RunOn::Host>(0.0);
          v_s = local_v_s_fab.array();

          local_w_s_fab.resize(zbx, ncomp);
          local_w_s_fab.setVal<RunOn::Host>(0.0);
          w_s = local_w_s_fab.array();
        }
#endif

        amrex::ParallelFor(nrp, [pstruct,p_realarray,p_lo,dxi,inv_reg_cell_vol,u_s, v_s, w_s]
        AMREX_GPU_DEVICE (int ip) noexcept
        {
          const ParticleType& p = pstruct[ip];

          const Real lx = (p.pos(0) - p_lo[0]) * dxi[0] + 0.5;
          const Real ly = (p.pos(1) - p_lo[1]) * dxi[1] + 0.5;
          const Real lz = (p.pos(2) - p_lo[2]) * dxi[2] + 0.5;

          const int i = static_cast<int>(amrex::Math::floor(lx));
          const int j = static_cast<int>(amrex::Math::floor(ly));
          const int k = static_cast<int>(amrex::Math::floor(lz));

          const Real wx_hi(lx - static_cast<Real>(i));
          const Real wy_hi(ly - static_cast<Real>(j));
          const Real wz_hi(lz - static_cast<Real>(k));

          const Real wx_lo(1.0 - wx_hi);
          const Real wy_lo(1.0 - wy_hi);
          const Real wz_lo(1.0 - wz_hi);

          const Real pvol = p_realarray[SoArealData::statwt][ip] *
            p_realarray[SoArealData::volume][ip] * inv_reg_cell_vol;

          {// Deposition of x velocity -- x-face deposition

            const Real pvelx = pvol*p_realarray[SoArealData::velx][ip];
#if USE_HAT_FUNCTION
            const Real lxc = (p.pos(0) - p_lo[0]) * dxi[0];

            const int ii = static_cast<int>(amrex::Math::floor(lxc));

            const Real wu_hi(lxc - static_cast<Real>(ii));
            const Real wu_lo(1.0 - wu_hi);

            HostDevice::Atomic::Add(&u_s(ii,   j-1, k-1, 0),wu_lo*wy_lo*wz_lo*pvelx);
            HostDevice::Atomic::Add(&u_s(ii,   j-1, k  , 0),wu_lo*wy_lo*wz_hi*pvelx);
            HostDevice::Atomic::Add(&u_s(ii,   j,   k-1, 0),wu_lo*wy_hi*wz_lo*pvelx);
            HostDevice::Atomic::Add(&u_s(ii,   j,   k  , 0),wu_lo*wy_hi*wz_hi*pvelx);
            HostDevice::Atomic::Add(&u_s(ii+1, j-1, k-1, 0),wu_hi*wy_lo*wz_lo*pvelx);
            HostDevice::Atomic::Add(&u_s(ii+1, j-1, k  , 0),wu_hi*wy_lo*wz_hi*pvelx);
            HostDevice::Atomic::Add(&u_s(ii+1, j,   k-1, 0),wu_hi*wy_hi*wz_lo*pvelx);
            HostDevice::Atomic::Add(&u_s(ii+1, j,   k  , 0),wu_hi*wy_hi*wz_hi*pvelx);

            HostDevice::Atomic::Add(&u_s(ii,   j-1, k-1, 1),wu_lo*wy_lo*wz_lo*pvol);
            HostDevice::Atomic::Add(&u_s(ii,   j-1, k  , 1),wu_lo*wy_lo*wz_hi*pvol);
            HostDevice::Atomic::Add(&u_s(ii,   j,   k-1, 1),wu_lo*wy_hi*wz_lo*pvol);
            HostDevice::Atomic::Add(&u_s(ii,   j,   k  , 1),wu_lo*wy_hi*wz_hi*pvol);
            HostDevice::Atomic::Add(&u_s(ii+1, j-1, k-1, 1),wu_hi*wy_lo*wz_lo*pvol);
            HostDevice::Atomic::Add(&u_s(ii+1, j-1, k  , 1),wu_hi*wy_lo*wz_hi*pvol);
            HostDevice::Atomic::Add(&u_s(ii+1, j,   k-1, 1),wu_hi*wy_hi*wz_lo*pvol);
            HostDevice::Atomic::Add(&u_s(ii+1, j,   k  , 1),wu_hi*wy_hi*wz_hi*pvol);
#else

            HostDevice::Atomic::Add(&u_s(i   , j-1, k-1, 0),wy_lo*wz_lo*pvelx);
            HostDevice::Atomic::Add(&u_s(i   , j-1, k  , 0),wy_lo*wz_hi*pvelx);
            HostDevice::Atomic::Add(&u_s(i   , j,   k-1, 0),wy_hi*wz_lo*pvelx);
            HostDevice::Atomic::Add(&u_s(i   , j,   k  , 0),wy_hi*wz_hi*pvelx);

            HostDevice::Atomic::Add(&u_s(i   , j-1, k-1, 1),wy_lo*wz_lo*pvol);
            HostDevice::Atomic::Add(&u_s(i   , j-1, k  , 1),wy_lo*wz_hi*pvol);
            HostDevice::Atomic::Add(&u_s(i   , j,   k-1, 1),wy_hi*wz_lo*pvol);
            HostDevice::Atomic::Add(&u_s(i   , j,   k  , 1),wy_hi*wz_hi*pvol);
#endif
          }


          {// Deposition of y velocity -- y-face deposition

            const Real pvely = pvol*p_realarray[SoArealData::vely][ip];
#if USE_HAT_FUNCTION
            const Real lyc = (p.pos(1) - p_lo[1]) * dxi[1];

            const int jj = static_cast<int>(amrex::Math::floor(lyc));

            const Real wv_hi(lyc - static_cast<Real>(jj));
            const Real wv_lo(1.0 - wv_hi);

            HostDevice::Atomic::Add(&v_s(i-1, jj,   k-1, 0),wx_lo*wv_lo*wz_lo*pvely);
            HostDevice::Atomic::Add(&v_s(i-1, jj,   k  , 0),wx_lo*wv_lo*wz_hi*pvely);
            HostDevice::Atomic::Add(&v_s(i-1, jj+1, k-1, 0),wx_lo*wv_hi*wz_lo*pvely);
            HostDevice::Atomic::Add(&v_s(i-1, jj+1, k  , 0),wx_lo*wv_hi*wz_hi*pvely);
            HostDevice::Atomic::Add(&v_s(i,   jj,   k-1, 0),wx_hi*wv_lo*wz_lo*pvely);
            HostDevice::Atomic::Add(&v_s(i,   jj,   k  , 0),wx_hi*wv_lo*wz_hi*pvely);
            HostDevice::Atomic::Add(&v_s(i,   jj+1, k-1, 0),wx_hi*wv_hi*wz_lo*pvely);
            HostDevice::Atomic::Add(&v_s(i,   jj+1, k  , 0),wx_hi*wv_hi*wz_hi*pvely);

            HostDevice::Atomic::Add(&v_s(i-1, jj,   k-1, 1),wx_lo*wv_lo*wz_lo*pvol);
            HostDevice::Atomic::Add(&v_s(i-1, jj,   k  , 1),wx_lo*wv_lo*wz_hi*pvol);
            HostDevice::Atomic::Add(&v_s(i-1, jj+1, k-1, 1),wx_lo*wv_hi*wz_lo*pvol);
            HostDevice::Atomic::Add(&v_s(i-1, jj+1, k  , 1),wx_lo*wv_hi*wz_hi*pvol);
            HostDevice::Atomic::Add(&v_s(i,   jj,   k-1, 1),wx_hi*wv_lo*wz_lo*pvol);
            HostDevice::Atomic::Add(&v_s(i,   jj,   k  , 1),wx_hi*wv_lo*wz_hi*pvol);
            HostDevice::Atomic::Add(&v_s(i,   jj+1, k-1, 1),wx_hi*wv_hi*wz_lo*pvol);
            HostDevice::Atomic::Add(&v_s(i,   jj+1, k  , 1),wx_hi*wv_hi*wz_hi*pvol);
#else
            HostDevice::Atomic::Add(&v_s(i-1, j ,   k-1, 0),wx_lo*wz_lo*pvely);
            HostDevice::Atomic::Add(&v_s(i-1, j ,   k  , 0),wx_lo*wz_hi*pvely);
            HostDevice::Atomic::Add(&v_s(i,   j ,   k-1, 0),wx_hi*wz_lo*pvely);
            HostDevice::Atomic::Add(&v_s(i,   j ,   k  , 0),wx_hi*wz_hi*pvely);

            HostDevice::Atomic::Add(&v_s(i-1, j ,   k-1, 1),wx_lo*wz_lo*pvol);
            HostDevice::Atomic::Add(&v_s(i-1, j ,   k  , 1),wx_lo*wz_hi*pvol);
            HostDevice::Atomic::Add(&v_s(i,   j ,   k-1, 1),wx_hi*wz_lo*pvol);
            HostDevice::Atomic::Add(&v_s(i,   j ,   k  , 1),wx_hi*wz_hi*pvol);

#endif
          }


          {// Deposition of z velocity -- z-face deposition

            const Real pvelz = pvol*p_realarray[SoArealData::velz][ip];
#if USE_HAT_FUNCTION
            const Real lzc = (p.pos(2) - p_lo[2]) * dxi[2];

            const int kk = static_cast<int>(amrex::Math::floor(lzc));

            const Real ww_hi(lzc - static_cast<Real>(kk));
            const Real ww_lo(1.0 - ww_hi);

            HostDevice::Atomic::Add(&w_s(i-1, j-1, kk  , 0),wx_lo*wy_lo*ww_lo*pvelz);
            HostDevice::Atomic::Add(&w_s(i-1, j-1, kk+1, 0),wx_lo*wy_lo*ww_hi*pvelz);
            HostDevice::Atomic::Add(&w_s(i-1, j,   kk  , 0),wx_lo*wy_hi*ww_lo*pvelz);
            HostDevice::Atomic::Add(&w_s(i-1, j,   kk+1, 0),wx_lo*wy_hi*ww_hi*pvelz);
            HostDevice::Atomic::Add(&w_s(i,   j-1, kk  , 0),wx_hi*wy_lo*ww_lo*pvelz);
            HostDevice::Atomic::Add(&w_s(i,   j-1, kk+1, 0),wx_hi*wy_lo*ww_hi*pvelz);
            HostDevice::Atomic::Add(&w_s(i,   j,   kk  , 0),wx_hi*wy_hi*ww_lo*pvelz);
            HostDevice::Atomic::Add(&w_s(i,   j,   kk+1, 0),wx_hi*wy_hi*ww_hi*pvelz);

            HostDevice::Atomic::Add(&w_s(i-1, j-1, kk  , 1),wx_lo*wy_lo*ww_lo*pvol);
            HostDevice::Atomic::Add(&w_s(i-1, j-1, kk+1, 1),wx_lo*wy_lo*ww_hi*pvol);
            HostDevice::Atomic::Add(&w_s(i-1, j,   kk  , 1),wx_lo*wy_hi*ww_lo*pvol);
            HostDevice::Atomic::Add(&w_s(i-1, j,   kk+1, 1),wx_lo*wy_hi*ww_hi*pvol);
            HostDevice::Atomic::Add(&w_s(i,   j-1, kk  , 1),wx_hi*wy_lo*ww_lo*pvol);
            HostDevice::Atomic::Add(&w_s(i,   j-1, kk+1, 1),wx_hi*wy_lo*ww_hi*pvol);
            HostDevice::Atomic::Add(&w_s(i,   j,   kk  , 1),wx_hi*wy_hi*ww_lo*pvol);
            HostDevice::Atomic::Add(&w_s(i,   j,   kk+1, 1),wx_hi*wy_hi*ww_hi*pvol);
#else
            HostDevice::Atomic::Add(&w_s(i-1, j-1, k   , 0),wx_lo*wy_lo*pvelz);
            HostDevice::Atomic::Add(&w_s(i-1, j,   k   , 0),wx_lo*wy_hi*pvelz);
            HostDevice::Atomic::Add(&w_s(i,   j-1, k   , 0),wx_hi*wy_lo*pvelz);
            HostDevice::Atomic::Add(&w_s(i,   j,   k   , 0),wx_hi*wy_hi*pvelz);

            HostDevice::Atomic::Add(&w_s(i-1, j-1, k   , 1),wx_lo*wy_lo*pvol);
            HostDevice::Atomic::Add(&w_s(i-1, j,   k   , 1),wx_lo*wy_hi*pvol);
            HostDevice::Atomic::Add(&w_s(i,   j-1, k   , 1),wx_hi*wy_lo*pvol);
            HostDevice::Atomic::Add(&w_s(i,   j,   k   , 1),wx_hi*wy_hi*pvol);
#endif
          }
        });
#ifdef _OPENMP
        if (Gpu::notInLaunchRegion())
        {
          u_s_fab.atomicAdd<RunOn::Host>(local_u_s_fab, xbx, xbx, 0, 0, ncomp);
          v_s_fab.atomicAdd<RunOn::Host>(local_v_s_fab, ybx, ybx, 0, 0, ncomp);
          w_s_fab.atomicAdd<RunOn::Host>(local_w_s_fab, zbx, zbx, 0, 0, ncomp);
        }
#endif

      }
    }
  }
}


void
MFIXParticleContainer::PICHydroStep (int lev,
                                     const bool apply_forces,
                                     const bool update_parcels,
                                     const bool use_taylor_approx,
                                     const Real advance_vel_p,
                                     Real dt,
                                     RealVect& gravity,
                                     Vector< Array<MultiFab*,3> >& vel_s_in,
                                     MultiFab & ep_s_out,
                                     Array<MultiFab*,3>& vel_s_out,
                                     const MultiFab * volfrac,
                                     const amrex::FabArray<EBCellFlagFab>* flags,
                                     EBFArrayBoxFactory* ebfactory,
                                     const int ls_refinement,
                                     const MultiFab* ls_phi)
{
  BL_PROFILE("MFIXParticleContainer::SolidsVolumeDeposition()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      p_lo = gm.ProbLoArray();
  const auto      p_hi = gm.ProbHiArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

  const auto      reg_cell_vol = dx[0]*dx[1]*dx[2];

  const Real en = (m_pic.damping_factor() + 1.0);

  const Real en_w = m_pic.damping_factor_wall_normal();
  const Real et_w = m_pic.damping_factor_wall_tangent();

  const Real vel_ref_frame = m_pic.vel_ref_frame();

  const Real ep_cp = m_pic.ep_cp();
  const Real inv_ep_cp = 1.0/ep_cp;

  const Real three_sqrt_two(3.0*std::sqrt(2.0));

  const int x_lo_bc = m_boundary_conditions.domain_bc(0);
  const int x_hi_bc = m_boundary_conditions.domain_bc(1);
  const int y_lo_bc = m_boundary_conditions.domain_bc(2);
  const int y_hi_bc = m_boundary_conditions.domain_bc(3);
  const int z_lo_bc = m_boundary_conditions.domain_bc(4);
  const int z_hi_bc = m_boundary_conditions.domain_bc(5);

  const int idx_vel_txfr = m_runtimeRealData.vel_txfr;

  const int solve_reactions = reactions.solve();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_fab_to_be_filled;

    FArrayBox local_u_s_fab;
    FArrayBox local_v_s_fab;
    FArrayBox local_w_s_fab;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

      PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& plev  = GetParticles(lev);
      auto& ptile = plev[index];

      auto& particles = pti.GetArrayOfStructs();
      ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();
      FArrayBox& fab_to_be_filled = ep_s_out[pti];

      FArrayBox& u_s_fab = (*vel_s_out[0])[pti];
      FArrayBox& v_s_fab = (*vel_s_out[1])[pti];
      FArrayBox& w_s_fab = (*vel_s_out[2])[pti];

      const Box& bx  = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(bx) != FabType::covered ) {

        auto volarr = fab_to_be_filled.array();
        const auto& flagsarr = (*flags)[pti].array();
        const auto& vfrac = (*volfrac)[pti].array();

        // Determine if this particle tile actually has any walls
        bool has_walls = false;

        if ((ebfactory != NULL) &&
           ((*flags)[pti].getType(amrex::grow(bx,1)) == FabType::singlevalued))  {

          has_walls = true;

        } else {
          // We need this test for the case of an inflow boundary:
          // inflow does not appear in the EBFactory but
          // the particles see it as a wall

          // Create the nodal refined box based on the current particle tile
          Box refined_box(amrex::convert(amrex::refine(bx,ls_refinement), IntVect{1,1,1}));

          // Set tol to 1/2 dx
          Real tol = amrex::min(dx[0], amrex::min(dx[1], dx[2])) / 2;

          Real ls_min_over_box = ((*ls_phi)[pti]).min<RunOn::Gpu>(refined_box,0);

          if (ls_min_over_box < tol) has_walls = true;

        }

        const auto& phiarr = (has_walls) ? ls_phi->const_array(pti) : Array4<Real const>{};

        //const auto& avg_prop_array = avg_prop_in[lev]->array(pti);
        const auto& u_so = (*vel_s_in[lev][0]).const_array(pti);
        const auto& v_so = (*vel_s_in[lev][1]).const_array(pti);
        const auto& w_so = (*vel_s_in[lev][2]).const_array(pti);

        auto u_s = u_s_fab.array();
        auto v_s = v_s_fab.array();
        auto w_s = w_s_fab.array();
#ifdef _OPENMP
        const int mf_ncomp = ep_s_out.nComp();

        Box tile_box = bx;

        const int vel_ncomp = vel_s_out[0]->nComp();
        AMREX_ALWAYS_ASSERT(vel_ncomp == vel_s_out[1]->nComp());
        AMREX_ALWAYS_ASSERT(vel_ncomp == vel_s_out[2]->nComp());

        const int ngrow = vel_s_in[lev][0]->nGrow();
        Box const xbx = Box(bx).grow(ngrow).surroundingNodes(0);
        Box const ybx = Box(bx).grow(ngrow).surroundingNodes(1);
        Box const zbx = Box(bx).grow(ngrow).surroundingNodes(2);

        if(Gpu::notInLaunchRegion())
        {
          tile_box.grow(ep_s_out.nGrow());
          local_fab_to_be_filled.resize(tile_box, mf_ncomp);
          local_fab_to_be_filled.setVal<RunOn::Host>(0.0);
          volarr = local_fab_to_be_filled.array();

          local_u_s_fab.resize(xbx, vel_ncomp);
          local_u_s_fab.setVal<RunOn::Host>(0.0);
          u_s = local_u_s_fab.array();

          local_v_s_fab.resize(ybx, vel_ncomp);
          local_v_s_fab.setVal<RunOn::Host>(0.0);
          v_s = local_v_s_fab.array();

          local_w_s_fab.resize(zbx, vel_ncomp);
          local_w_s_fab.setVal<RunOn::Host>(0.0);
          w_s = local_w_s_fab.array();
        }
#endif

        auto ptile_data = ptile.getParticleTileData();

        amrex::ParallelFor(nrp,
           [pstruct,p_realarray,p_hi,p_lo,dx,dxi,vfrac,volarr, u_so, v_so, w_so, en, ep_cp,
            reg_cell_vol,flagsarr, dt, gravity, has_walls,ls_refinement,phiarr,
            vel_ref_frame, three_sqrt_two, en_w, et_w, u_s, v_s, w_s, inv_ep_cp,
            apply_forces, update_parcels, use_taylor_approx, advance_vel_p,
            x_lo_bc,x_hi_bc, y_lo_bc,y_hi_bc,z_lo_bc,z_hi_bc, idx_vel_txfr,
            solve_reactions, ptile_data]
          AMREX_GPU_DEVICE (int ip) noexcept
          {
            ParticleType& p = pstruct[ip];

            const amrex::RealVect vel_p_old  = {p_realarray[SoArealData::velx][ip],
                                                p_realarray[SoArealData::vely][ip],
                                                p_realarray[SoArealData::velz][ip]};

            // position
            RealVect pos(p.pos());
            RealVect vel_p(vel_p_old);

            if (apply_forces ) {
              // solids volume fraction
              const amrex::Real eps_p = amrex::max(1.0e-8, p_realarray[SoArealData::oneOverI][ip]);
              const amrex::Real density_p = p_realarray[SoArealData::density][ip];

              // inverse of particle mass
              const Real inv_mass = 1.0/p_realarray[SoArealData::mass][ip];

              // beta_p := (drag coeff * volume) / mass
              //        :=  drag coeff / density
              const Real beta_p = p_realarray[SoArealData::dragcoeff][ip] * inv_mass;

              amrex::RealVect mass_txfr_p(0.);

              if (solve_reactions) {
                mass_txfr_p = {ptile_data.m_runtime_rdata[idx_vel_txfr+0][ip] * inv_mass,
                               ptile_data.m_runtime_rdata[idx_vel_txfr+1][ip] * inv_mass,
                               ptile_data.m_runtime_rdata[idx_vel_txfr+2][ip] * inv_mass};
              }

              // beta_vel_g := ((drag coeff * volume) / mass ) * u_gp
              //            := (drag coeff / density) * u_gp
              const amrex::RealVect beta_vel_fp = {p_realarray[SoArealData::dragx][ip] * inv_mass,
                                                   p_realarray[SoArealData::dragy][ip] * inv_mass,
                                                   p_realarray[SoArealData::dragz][ip] * inv_mass};

              // solids stress gradient:
              // grad_tau_p = (volume * grad_tau_p) / (mass * ep_s))
              const amrex::RealVect grad_tau_p = {p_realarray[SoArealData::omegax][ip],
                                                  p_realarray[SoArealData::omegay][ip],
                                                  p_realarray[SoArealData::omegaz][ip]};

              const Real mfp_vel = (p_realarray[SoArealData::radius][ip] /
                                    (dt * three_sqrt_two * eps_p));


              // Compute the updated PIC velocity
              vel_p = updated_pic_velocity(pos, vel_p_old, density_p,
                 grad_tau_p, beta_p, mass_txfr_p, beta_vel_fp, u_so, v_so, w_so,
                 en, eps_p, ep_cp, vel_ref_frame, mfp_vel, dt, gravity, dxi, p_lo);

              // Update parcel positions
              pos[0] += dt * vel_p[0];
              pos[1] += dt * vel_p[1];
              pos[2] += dt * vel_p[2];

              // Effective radius of the parcel
              Real eff_radius = p_realarray[SoArealData::radius][ip] *
                std::cbrt(p_realarray[SoArealData::statwt][ip] * inv_ep_cp);

              // If this FAB has EB, reflect any parcels overlapping the levelset
              // back into the domain.
              if (has_walls) {

                Real ls_value = interp_level_set(pos, ls_refinement, phiarr, p_lo, dxi);

                const Real overlap = eff_radius - ls_value;

                // The particle intersects the levelset.
                if (overlap > 0.) {

                  RealVect normal(0.);
                  level_set_normal(pos, ls_refinement, normal, phiarr, p_lo, dxi);

                  // Reflect the particle.
                  pos[0] += overlap*normal[0];
                  pos[1] += overlap*normal[1];
                  pos[2] += overlap*normal[2];

                  // Plane ref point
                  const Real Nw_Vp = normal[0]*vel_p[0] + normal[1]*vel_p[1] + normal[2]*vel_p[2];

                  // Parcel normal velocity
                  const RealVect Vpn = {Nw_Vp*normal[0], Nw_Vp*normal[1], Nw_Vp*normal[2]};

                  // Parcel tangential velocity
                  const RealVect Vpt = {vel_p[0]-Vpn[0], vel_p[1]-Vpn[1], vel_p[2]-Vpn[2]};

                  // Rebound parcel if moving towards wall.
                  if(Nw_Vp < 0.) {
                    vel_p[0] = -en_w*Vpn[0] + et_w*Vpt[0];
                    vel_p[1] = -en_w*Vpn[1] + et_w*Vpt[1];
                    vel_p[2] = -en_w*Vpn[2] + et_w*Vpt[2];

                  } else {
                    vel_p[0] = Vpn[0] + et_w*Vpt[0];
                    vel_p[1] = Vpn[1] + et_w*Vpt[1];
                    vel_p[2] = Vpn[2] + et_w*Vpt[2];
                  }

                }
              } // end has_walls

              // Impose domain constraints.
              if (x_lo_bc && pos[0] < p_lo[0]+eff_radius)
                pos[0] = p_lo[0] + eff_radius;

              if (x_hi_bc && pos[0] + eff_radius > p_hi[0])
                pos[0] = p_hi[0] - eff_radius;

              if (y_lo_bc && pos[1] < p_lo[1]+eff_radius)
                pos[1] = p_lo[1] + eff_radius;

              if (y_hi_bc && pos[1] + eff_radius > p_hi[1])
                pos[1] = p_hi[1] - eff_radius;

              if (z_lo_bc && pos[2] < p_lo[2]+eff_radius)
                pos[2] = p_lo[2] + eff_radius;

              if (z_hi_bc && pos[2] + eff_radius > p_hi[2])
                pos[2] = p_hi[2] - eff_radius;
            }

            amrex::Real x = (pos[0] - p_lo[0]) * dxi[0] + 0.5;
            amrex::Real y = (pos[1] - p_lo[1]) * dxi[1] + 0.5;
            amrex::Real z = (pos[2] - p_lo[2]) * dxi[2] + 0.5;

            int i = static_cast<int>(amrex::Math::floor(x));
            int j = static_cast<int>(amrex::Math::floor(y));
            int k = static_cast<int>(amrex::Math::floor(z));

            amrex::GpuArray<amrex::Real,2> wx;
            amrex::GpuArray<amrex::Real,2> wy;
            amrex::GpuArray<amrex::Real,2> wz;

            wx[1] = x - static_cast<Real>(i);
            wy[1] = y - static_cast<Real>(j);
            wz[1] = z - static_cast<Real>(k);

            wx[0] = 1.0 - wx[1];
            wy[0] = 1.0 - wy[1];
            wz[0] = 1.0 - wz[1];

            amrex::Real total_weight = 0.0;
            GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

            if (use_taylor_approx) {

              const amrex::GpuArray<amrex::Real,2> ws({-1.0,1.0});

              const amrex::RealVect dt_vel_dxi =
                {dt*p_realarray[SoArealData::velx][ip]*dxi[0],
                 dt*p_realarray[SoArealData::vely][ip]*dxi[1],
                 dt*p_realarray[SoArealData::velz][ip]*dxi[2]};

              for (int ii = 0; ii <= 1; ++ii)
                for (int jj = 0; jj <= 1; ++jj)
                  for (int kk = 0; kk <= 1; ++kk){
                    if( !flagsarr(i-1+ii,j-1+jj,k-1+kk).isCovered() ) {

                      weights[ii][jj][kk] = wx[ii]*wy[jj]*wz[kk] +
                        ws[ii]*wy[jj]*wz[kk]*dt_vel_dxi[0] +
                        wx[ii]*ws[jj]*wz[kk]*dt_vel_dxi[1] +
                        wx[ii]*wy[jj]*ws[kk]*dt_vel_dxi[2];

                      total_weight += weights[ii][jj][kk];

                    } else {
                      weights[ii][jj][kk] = 0.0;
                    }
                  }

            } else {

              for (int ii = 0; ii <= 1; ++ii)
                for (int jj = 0; jj <= 1; ++jj)
                  for (int kk = 0; kk <= 1; ++kk){
                    if( !flagsarr(i-1+ii,j-1+jj,k-1+kk).isCovered() ) {
                      weights[ii][jj][kk] = wx[ii]*wy[jj]*wz[kk];
                      total_weight += weights[ii][jj][kk];

                    } else {
                      weights[ii][jj][kk] = 0.0;
                    }
                  }
            }

            for (int ii = 0; ii <= 1; ++ii)
              for (int jj = 0; jj <= 1; ++jj)
                for (int kk = 0; kk <= 1; ++kk)
                  weights[ii][jj][kk] /= total_weight;

            const Real pvol = p_realarray[SoArealData::statwt][ip] *
              p_realarray[SoArealData::volume][ip] / reg_cell_vol;

            for (int kk = -1; kk <= 0; ++kk) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int ii = -1; ii <= 0; ++ii) {
                  if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                    continue;

                  Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);
                  HostDevice::Atomic::Add(&volarr(i+ii,j+jj,k+kk), weight_vol*pvol);
                }
              }
            }

            const Real new_fac = advance_vel_p;
            const Real old_fac = 1.0 - new_fac;

            x = ((p.pos(0) + new_fac * dt * vel_p[0]) - p_lo[0]) * dxi[0] + 0.5;
            y = ((p.pos(1) + new_fac * dt * vel_p[1]) - p_lo[1]) * dxi[1] + 0.5;
            z = ((p.pos(2) + new_fac * dt * vel_p[2]) - p_lo[2]) * dxi[2] + 0.5;

            i = static_cast<int>(amrex::Math::floor(x));
            j = static_cast<int>(amrex::Math::floor(y));
            k = static_cast<int>(amrex::Math::floor(z));

            const Real wx_hi = x - static_cast<Real>(i);
            const Real wy_hi = y - static_cast<Real>(j);
            const Real wz_hi = z - static_cast<Real>(k);

            const Real wx_lo = 1.0 - wx_hi;
            const Real wy_lo = 1.0 - wy_hi;
            const Real wz_lo = 1.0 - wz_hi;


            {// Deposition of x velocity -- x-face deposition

              const Real pvelx = pvol*(new_fac*vel_p[0] + old_fac*vel_p_old[0]);
#if USE_HAT_FUNCTION
              const Real lxc = (pos[0] - p_lo[0]) * dxi[0];

              const int ii = static_cast<int>(amrex::Math::floor(lxc));

              const Real wu_hi(lxc - static_cast<Real>(ii));
              const Real wu_lo(1.0 - wu_hi);

              HostDevice::Atomic::Add(&u_s(ii,   j-1, k-1, 0),wu_lo*wy_lo*wz_lo*pvelx);
              HostDevice::Atomic::Add(&u_s(ii,   j-1, k  , 0),wu_lo*wy_lo*wz_hi*pvelx);
              HostDevice::Atomic::Add(&u_s(ii,   j,   k-1, 0),wu_lo*wy_hi*wz_lo*pvelx);
              HostDevice::Atomic::Add(&u_s(ii,   j,   k  , 0),wu_lo*wy_hi*wz_hi*pvelx);
              HostDevice::Atomic::Add(&u_s(ii+1, j-1, k-1, 0),wu_hi*wy_lo*wz_lo*pvelx);
              HostDevice::Atomic::Add(&u_s(ii+1, j-1, k  , 0),wu_hi*wy_lo*wz_hi*pvelx);
              HostDevice::Atomic::Add(&u_s(ii+1, j,   k-1, 0),wu_hi*wy_hi*wz_lo*pvelx);
              HostDevice::Atomic::Add(&u_s(ii+1, j,   k  , 0),wu_hi*wy_hi*wz_hi*pvelx);

              HostDevice::Atomic::Add(&u_s(ii,   j-1, k-1, 1),wu_lo*wy_lo*wz_lo*pvol);
              HostDevice::Atomic::Add(&u_s(ii,   j-1, k  , 1),wu_lo*wy_lo*wz_hi*pvol);
              HostDevice::Atomic::Add(&u_s(ii,   j,   k-1, 1),wu_lo*wy_hi*wz_lo*pvol);
              HostDevice::Atomic::Add(&u_s(ii,   j,   k  , 1),wu_lo*wy_hi*wz_hi*pvol);
              HostDevice::Atomic::Add(&u_s(ii+1, j-1, k-1, 1),wu_hi*wy_lo*wz_lo*pvol);
              HostDevice::Atomic::Add(&u_s(ii+1, j-1, k  , 1),wu_hi*wy_lo*wz_hi*pvol);
              HostDevice::Atomic::Add(&u_s(ii+1, j,   k-1, 1),wu_hi*wy_hi*wz_lo*pvol);
              HostDevice::Atomic::Add(&u_s(ii+1, j,   k  , 1),wu_hi*wy_hi*wz_hi*pvol);
#else

              HostDevice::Atomic::Add(&u_s(i   , j-1, k-1, 0),wy_lo*wz_lo*pvelx);
              HostDevice::Atomic::Add(&u_s(i   , j-1, k  , 0),wy_lo*wz_hi*pvelx);
              HostDevice::Atomic::Add(&u_s(i   , j,   k-1, 0),wy_hi*wz_lo*pvelx);
              HostDevice::Atomic::Add(&u_s(i   , j,   k  , 0),wy_hi*wz_hi*pvelx);

              HostDevice::Atomic::Add(&u_s(i   , j-1, k-1, 1),wy_lo*wz_lo*pvol);
              HostDevice::Atomic::Add(&u_s(i   , j-1, k  , 1),wy_lo*wz_hi*pvol);
              HostDevice::Atomic::Add(&u_s(i   , j,   k-1, 1),wy_hi*wz_lo*pvol);
              HostDevice::Atomic::Add(&u_s(i   , j,   k  , 1),wy_hi*wz_hi*pvol);
#endif
            }


            {// Deposition of y velocity -- y-face deposition

              const Real pvely = pvol*(new_fac*vel_p[1] + old_fac*vel_p_old[1]);
#if USE_HAT_FUNCTION
              const Real lyc = (pos[1] - p_lo[1]) * dxi[1];

              const int jj = static_cast<int>(amrex::Math::floor(lyc));

              const Real wv_hi(lyc - static_cast<Real>(jj));
              const Real wv_lo(1.0 - wv_hi);

              HostDevice::Atomic::Add(&v_s(i-1, jj,   k-1, 0),wx_lo*wv_lo*wz_lo*pvely);
              HostDevice::Atomic::Add(&v_s(i-1, jj,   k  , 0),wx_lo*wv_lo*wz_hi*pvely);
              HostDevice::Atomic::Add(&v_s(i-1, jj+1, k-1, 0),wx_lo*wv_hi*wz_lo*pvely);
              HostDevice::Atomic::Add(&v_s(i-1, jj+1, k  , 0),wx_lo*wv_hi*wz_hi*pvely);
              HostDevice::Atomic::Add(&v_s(i,   jj,   k-1, 0),wx_hi*wv_lo*wz_lo*pvely);
              HostDevice::Atomic::Add(&v_s(i,   jj,   k  , 0),wx_hi*wv_lo*wz_hi*pvely);
              HostDevice::Atomic::Add(&v_s(i,   jj+1, k-1, 0),wx_hi*wv_hi*wz_lo*pvely);
              HostDevice::Atomic::Add(&v_s(i,   jj+1, k  , 0),wx_hi*wv_hi*wz_hi*pvely);

              HostDevice::Atomic::Add(&v_s(i-1, jj,   k-1, 1),wx_lo*wv_lo*wz_lo*pvol);
              HostDevice::Atomic::Add(&v_s(i-1, jj,   k  , 1),wx_lo*wv_lo*wz_hi*pvol);
              HostDevice::Atomic::Add(&v_s(i-1, jj+1, k-1, 1),wx_lo*wv_hi*wz_lo*pvol);
              HostDevice::Atomic::Add(&v_s(i-1, jj+1, k  , 1),wx_lo*wv_hi*wz_hi*pvol);
              HostDevice::Atomic::Add(&v_s(i,   jj,   k-1, 1),wx_hi*wv_lo*wz_lo*pvol);
              HostDevice::Atomic::Add(&v_s(i,   jj,   k  , 1),wx_hi*wv_lo*wz_hi*pvol);
              HostDevice::Atomic::Add(&v_s(i,   jj+1, k-1, 1),wx_hi*wv_hi*wz_lo*pvol);
              HostDevice::Atomic::Add(&v_s(i,   jj+1, k  , 1),wx_hi*wv_hi*wz_hi*pvol);
#else
              HostDevice::Atomic::Add(&v_s(i-1, j ,   k-1, 0),wx_lo*wz_lo*pvely);
              HostDevice::Atomic::Add(&v_s(i-1, j ,   k  , 0),wx_lo*wz_hi*pvely);
              HostDevice::Atomic::Add(&v_s(i,   j ,   k-1, 0),wx_hi*wz_lo*pvely);
              HostDevice::Atomic::Add(&v_s(i,   j ,   k  , 0),wx_hi*wz_hi*pvely);

              HostDevice::Atomic::Add(&v_s(i-1, j ,   k-1, 1),wx_lo*wz_lo*pvol);
              HostDevice::Atomic::Add(&v_s(i-1, j ,   k  , 1),wx_lo*wz_hi*pvol);
              HostDevice::Atomic::Add(&v_s(i,   j ,   k-1, 1),wx_hi*wz_lo*pvol);
              HostDevice::Atomic::Add(&v_s(i,   j ,   k  , 1),wx_hi*wz_hi*pvol);
#endif
            }


            {// Deposition of z velocity -- z-face deposition

              const Real pvelz = pvol*(new_fac*vel_p[2] + old_fac*vel_p_old[2]);
#if USE_HAT_FUNCTION
              const Real lzc = (pos[2] - p_lo[2]) * dxi[2];

              const int kk = static_cast<int>(amrex::Math::floor(lzc));

              const Real ww_hi(lzc - static_cast<Real>(kk));
              const Real ww_lo(1.0 - ww_hi);

              HostDevice::Atomic::Add(&w_s(i-1, j-1, kk  , 0),wx_lo*wy_lo*ww_lo*pvelz);
              HostDevice::Atomic::Add(&w_s(i-1, j-1, kk+1, 0),wx_lo*wy_lo*ww_hi*pvelz);
              HostDevice::Atomic::Add(&w_s(i-1, j,   kk  , 0),wx_lo*wy_hi*ww_lo*pvelz);
              HostDevice::Atomic::Add(&w_s(i-1, j,   kk+1, 0),wx_lo*wy_hi*ww_hi*pvelz);
              HostDevice::Atomic::Add(&w_s(i,   j-1, kk  , 0),wx_hi*wy_lo*ww_lo*pvelz);
              HostDevice::Atomic::Add(&w_s(i,   j-1, kk+1, 0),wx_hi*wy_lo*ww_hi*pvelz);
              HostDevice::Atomic::Add(&w_s(i,   j,   kk  , 0),wx_hi*wy_hi*ww_lo*pvelz);
              HostDevice::Atomic::Add(&w_s(i,   j,   kk+1, 0),wx_hi*wy_hi*ww_hi*pvelz);

              HostDevice::Atomic::Add(&w_s(i-1, j-1, kk  , 1),wx_lo*wy_lo*ww_lo*pvol);
              HostDevice::Atomic::Add(&w_s(i-1, j-1, kk+1, 1),wx_lo*wy_lo*ww_hi*pvol);
              HostDevice::Atomic::Add(&w_s(i-1, j,   kk  , 1),wx_lo*wy_hi*ww_lo*pvol);
              HostDevice::Atomic::Add(&w_s(i-1, j,   kk+1, 1),wx_lo*wy_hi*ww_hi*pvol);
              HostDevice::Atomic::Add(&w_s(i,   j-1, kk  , 1),wx_hi*wy_lo*ww_lo*pvol);
              HostDevice::Atomic::Add(&w_s(i,   j-1, kk+1, 1),wx_hi*wy_lo*ww_hi*pvol);
              HostDevice::Atomic::Add(&w_s(i,   j,   kk  , 1),wx_hi*wy_hi*ww_lo*pvol);
              HostDevice::Atomic::Add(&w_s(i,   j,   kk+1, 1),wx_hi*wy_hi*ww_hi*pvol);
#else
              HostDevice::Atomic::Add(&w_s(i-1, j-1, k   , 0),wx_lo*wy_lo*pvelz);
              HostDevice::Atomic::Add(&w_s(i-1, j,   k   , 0),wx_lo*wy_hi*pvelz);
              HostDevice::Atomic::Add(&w_s(i,   j-1, k   , 0),wx_hi*wy_lo*pvelz);
              HostDevice::Atomic::Add(&w_s(i,   j,   k   , 0),wx_hi*wy_hi*pvelz);

              HostDevice::Atomic::Add(&w_s(i-1, j-1, k   , 1),wx_lo*wy_lo*pvol);
              HostDevice::Atomic::Add(&w_s(i-1, j,   k   , 1),wx_lo*wy_hi*pvol);
              HostDevice::Atomic::Add(&w_s(i,   j-1, k   , 1),wx_hi*wy_lo*pvol);
              HostDevice::Atomic::Add(&w_s(i,   j,   k   , 1),wx_hi*wy_hi*pvol);
#endif
            }


            if (update_parcels) {

              // Update positions
              p.pos(0) = pos[0];
              p.pos(1) = pos[1];
              p.pos(2) = pos[2];

              // update parcel velocity
              p_realarray[SoArealData::velx][ip] = vel_p[0];
              p_realarray[SoArealData::vely][ip] = vel_p[1];
              p_realarray[SoArealData::velz][ip] = vel_p[2];
            }

          });

        Gpu::synchronize();

#ifdef _OPENMP
        if(Gpu::notInLaunchRegion())
        {
          fab_to_be_filled.atomicAdd<RunOn::Host>(local_fab_to_be_filled,
              tile_box, tile_box, 0, 0, mf_ncomp);

          u_s_fab.atomicAdd<RunOn::Host>(local_u_s_fab, xbx, xbx, 0, 0, vel_ncomp);
          v_s_fab.atomicAdd<RunOn::Host>(local_v_s_fab, ybx, ybx, 0, 0, vel_ncomp);
          w_s_fab.atomicAdd<RunOn::Host>(local_w_s_fab, zbx, zbx, 0, 0, vel_ncomp);
        }
#endif

      }
    }
  }
}
