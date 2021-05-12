#include <AMReX.H>
#include "AMReX_Particles.H"
#include <mfix_pc.H>

#include <mfix_deposition_K.H>
#include <mfix.H>
#include <mfix_des_K.H>

using namespace amrex;

void MFIXParticleContainer::
MFIX_PC_SolidsVelocityDeposition (int lev,
                                  Array<MultiFab*,3>& vel_s_in,
                                  const FabArray<EBCellFlagFab>* flags)
{
  BL_PROFILE("(MFIXParticleContainer::MFIX_PC_SolidsVelocityDeposition)");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
  const auto      dxi = gm.InvCellSizeArray();

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

        Box const& xbx = Box(box).grow(ngrow).surroundingNodes(0);
        Box const& ybx = Box(box).grow(ngrow).surroundingNodes(1);
        Box const& zbx = Box(box).grow(ngrow).surroundingNodes(2);

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

        amrex::ParallelFor(nrp, [pstruct,p_realarray,plo,dxi,u_s, v_s, w_s]
        AMREX_GPU_DEVICE (int ip) noexcept
        {
          const ParticleType& p = pstruct[ip];

          const Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
          const Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
          const Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

          const int i = static_cast<int>(amrex::Math::floor(lx));
          const int j = static_cast<int>(amrex::Math::floor(ly));
          const int k = static_cast<int>(amrex::Math::floor(lz));

          const Real wx_hi(lx - static_cast<Real>(i));
          const Real wy_hi(ly - static_cast<Real>(j));
          const Real wz_hi(lz - static_cast<Real>(k));

          const Real wx_lo(1.0 - wx_hi);
          const Real wy_lo(1.0 - wy_hi);
          const Real wz_lo(1.0 - wz_hi);

          const Real pmass = p_realarray[SoArealData::statwt][ip] *
            p_realarray[SoArealData::mass][ip];

          {// Deposition of x velocity -- x-face deposition

            const Real pvelx = pmass*p_realarray[SoArealData::velx][ip];

            const Real lxc = (p.pos(0) - plo[0]) * dxi[0];

            const int ii = static_cast<int>(amrex::Math::floor(lxc));

            const Real wu_hi(lxc - static_cast<Real>(ii));
            const Real wu_lo(1.0 - wu_hi);

            amrex::Gpu::Atomic::Add(&u_s(ii,   j-1, k-1, 0),wu_lo*wy_lo*wz_lo*pvelx);
            amrex::Gpu::Atomic::Add(&u_s(ii,   j-1, k  , 0),wu_lo*wy_lo*wz_hi*pvelx);
            amrex::Gpu::Atomic::Add(&u_s(ii,   j,   k-1, 0),wu_lo*wy_hi*wz_lo*pvelx);
            amrex::Gpu::Atomic::Add(&u_s(ii,   j,   k  , 0),wu_lo*wy_hi*wz_hi*pvelx);
            amrex::Gpu::Atomic::Add(&u_s(ii+1, j-1, k-1, 0),wu_hi*wy_lo*wz_lo*pvelx);
            amrex::Gpu::Atomic::Add(&u_s(ii+1, j-1, k  , 0),wu_hi*wy_lo*wz_hi*pvelx);
            amrex::Gpu::Atomic::Add(&u_s(ii+1, j,   k-1, 0),wu_hi*wy_hi*wz_lo*pvelx);
            amrex::Gpu::Atomic::Add(&u_s(ii+1, j,   k  , 0),wu_hi*wy_hi*wz_hi*pvelx);

            amrex::Gpu::Atomic::Add(&u_s(ii,   j-1, k-1, 1),wu_lo*wy_lo*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&u_s(ii,   j-1, k  , 1),wu_lo*wy_lo*wz_hi*pmass);
            amrex::Gpu::Atomic::Add(&u_s(ii,   j,   k-1, 1),wu_lo*wy_hi*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&u_s(ii,   j,   k  , 1),wu_lo*wy_hi*wz_hi*pmass);
            amrex::Gpu::Atomic::Add(&u_s(ii+1, j-1, k-1, 1),wu_hi*wy_lo*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&u_s(ii+1, j-1, k  , 1),wu_hi*wy_lo*wz_hi*pmass);
            amrex::Gpu::Atomic::Add(&u_s(ii+1, j,   k-1, 1),wu_hi*wy_hi*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&u_s(ii+1, j,   k  , 1),wu_hi*wy_hi*wz_hi*pmass);
          }


          {// Deposition of y velocity -- y-face deposition

            const Real pvely = pmass*p_realarray[SoArealData::vely][ip];

            const Real lyc = (p.pos(1) - plo[1]) * dxi[1];

            const int jj = static_cast<int>(amrex::Math::floor(lyc));

            const Real wv_hi(lyc - static_cast<Real>(jj));
            const Real wv_lo(1.0 - wv_hi);

            amrex::Gpu::Atomic::Add(&v_s(i-1, jj,   k-1, 0),wx_lo*wv_lo*wz_lo*pvely);
            amrex::Gpu::Atomic::Add(&v_s(i-1, jj,   k  , 0),wx_lo*wv_lo*wz_hi*pvely);
            amrex::Gpu::Atomic::Add(&v_s(i-1, jj+1, k-1, 0),wx_lo*wv_hi*wz_lo*pvely);
            amrex::Gpu::Atomic::Add(&v_s(i-1, jj+1, k  , 0),wx_lo*wv_hi*wz_hi*pvely);
            amrex::Gpu::Atomic::Add(&v_s(i,   jj,   k-1, 0),wx_hi*wv_lo*wz_lo*pvely);
            amrex::Gpu::Atomic::Add(&v_s(i,   jj,   k  , 0),wx_hi*wv_lo*wz_hi*pvely);
            amrex::Gpu::Atomic::Add(&v_s(i,   jj+1, k-1, 0),wx_hi*wv_hi*wz_lo*pvely);
            amrex::Gpu::Atomic::Add(&v_s(i,   jj+1, k  , 0),wx_hi*wv_hi*wz_hi*pvely);

            amrex::Gpu::Atomic::Add(&v_s(i-1, jj,   k-1, 1),wx_lo*wv_lo*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&v_s(i-1, jj,   k  , 1),wx_lo*wv_lo*wz_hi*pmass);
            amrex::Gpu::Atomic::Add(&v_s(i-1, jj+1, k-1, 1),wx_lo*wv_hi*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&v_s(i-1, jj+1, k  , 1),wx_lo*wv_hi*wz_hi*pmass);
            amrex::Gpu::Atomic::Add(&v_s(i,   jj,   k-1, 1),wx_hi*wv_lo*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&v_s(i,   jj,   k  , 1),wx_hi*wv_lo*wz_hi*pmass);
            amrex::Gpu::Atomic::Add(&v_s(i,   jj+1, k-1, 1),wx_hi*wv_hi*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&v_s(i,   jj+1, k  , 1),wx_hi*wv_hi*wz_hi*pmass);
          }


          {// Deposition of z velocity -- z-face deposition

            const Real pvelz = pmass*p_realarray[SoArealData::velz][ip];
            const Real lzc = (p.pos(2) - plo[2]) * dxi[2];

            const int kk = static_cast<int>(amrex::Math::floor(lzc));

            const Real ww_hi(lzc - static_cast<Real>(kk));
            const Real ww_lo(1.0 - ww_hi);

            amrex::Gpu::Atomic::Add(&w_s(i-1, j-1, kk  , 0),wx_lo*wy_lo*ww_lo*pvelz);
            amrex::Gpu::Atomic::Add(&w_s(i-1, j-1, kk+1, 0),wx_lo*wy_lo*ww_hi*pvelz);
            amrex::Gpu::Atomic::Add(&w_s(i-1, j,   kk  , 0),wx_lo*wy_hi*ww_lo*pvelz);
            amrex::Gpu::Atomic::Add(&w_s(i-1, j,   kk+1, 0),wx_lo*wy_hi*ww_hi*pvelz);
            amrex::Gpu::Atomic::Add(&w_s(i,   j-1, kk  , 0),wx_hi*wy_lo*ww_lo*pvelz);
            amrex::Gpu::Atomic::Add(&w_s(i,   j-1, kk+1, 0),wx_hi*wy_lo*ww_hi*pvelz);
            amrex::Gpu::Atomic::Add(&w_s(i,   j,   kk  , 0),wx_hi*wy_hi*ww_lo*pvelz);
            amrex::Gpu::Atomic::Add(&w_s(i,   j,   kk+1, 0),wx_hi*wy_hi*ww_hi*pvelz);

            amrex::Gpu::Atomic::Add(&w_s(i-1, j-1, kk  , 1),wx_lo*wy_lo*ww_lo*pmass);
            amrex::Gpu::Atomic::Add(&w_s(i-1, j-1, kk+1, 1),wx_lo*wy_lo*ww_hi*pmass);
            amrex::Gpu::Atomic::Add(&w_s(i-1, j,   kk  , 1),wx_lo*wy_hi*ww_lo*pmass);
            amrex::Gpu::Atomic::Add(&w_s(i-1, j,   kk+1, 1),wx_lo*wy_hi*ww_hi*pmass);
            amrex::Gpu::Atomic::Add(&w_s(i,   j-1, kk  , 1),wx_hi*wy_lo*ww_lo*pmass);
            amrex::Gpu::Atomic::Add(&w_s(i,   j-1, kk+1, 1),wx_hi*wy_lo*ww_hi*pmass);
            amrex::Gpu::Atomic::Add(&w_s(i,   j,   kk  , 1),wx_hi*wy_hi*ww_lo*pmass);
            amrex::Gpu::Atomic::Add(&w_s(i,   j,   kk+1, 1),wx_hi*wy_hi*ww_hi*pmass);
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






#if 0
        if ( cell_centered ) {

          amrex::ParallelFor(nrp, [pstruct,p_realarray,plo,dxi,vel_s_arr]
          AMREX_GPU_DEVICE (int ip) noexcept
          {
            const ParticleType& p = pstruct[ip];

            const Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
            const Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
            const Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

            const int i = static_cast<int>(amrex::Math::floor(lx));
            const int j = static_cast<int>(amrex::Math::floor(ly));
            const int k = static_cast<int>(amrex::Math::floor(lz));

            const Real wx_hi(lx - static_cast<Real>(i));
            const Real wy_hi(ly - static_cast<Real>(j));
            const Real wz_hi(lz - static_cast<Real>(k));

            const Real wx_lo(1.0 - wx_hi);
            const Real wy_lo(1.0 - wy_hi);
            const Real wz_lo(1.0 - wz_hi);

            const Real pmass = p_realarray[SoArealData::statwt][ip] *
              p_realarray[SoArealData::mass][ip];

            const Real pvelx = pmass*p_realarray[SoArealData::velx][ip];
            const Real pvely = pmass*p_realarray[SoArealData::vely][ip];
            const Real pvelz = pmass*p_realarray[SoArealData::velz][ip];

            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, k-1, 0),wx_lo*wy_lo*wz_lo*pvelx);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, k  , 0),wx_lo*wy_lo*wz_hi*pvelx);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   k-1, 0),wx_lo*wy_hi*wz_lo*pvelx);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   k  , 0),wx_lo*wy_hi*wz_hi*pvelx);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, k-1, 0),wx_hi*wy_lo*wz_lo*pvelx);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, k  , 0),wx_hi*wy_lo*wz_hi*pvelx);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   k-1, 0),wx_hi*wy_hi*wz_lo*pvelx);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   k  , 0),wx_hi*wy_hi*wz_hi*pvelx);

            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, k-1, 1),wx_lo*wy_lo*wz_lo*pvely);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, k  , 1),wx_lo*wy_lo*wz_hi*pvely);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   k-1, 1),wx_lo*wy_hi*wz_lo*pvely);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   k  , 1),wx_lo*wy_hi*wz_hi*pvely);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, k-1, 1),wx_hi*wy_lo*wz_lo*pvely);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, k  , 1),wx_hi*wy_lo*wz_hi*pvely);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   k-1, 1),wx_hi*wy_hi*wz_lo*pvely);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   k  , 1),wx_hi*wy_hi*wz_hi*pvely);

            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, k-1, 0),wx_lo*wy_lo*wz_lo*pvelz);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, k  , 0),wx_lo*wy_lo*wz_hi*pvelz);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   k-1, 0),wx_lo*wy_hi*wz_lo*pvelz);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   k  , 0),wx_lo*wy_hi*wz_hi*pvelz);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, k-1, 0),wx_hi*wy_lo*wz_lo*pvelz);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, k  , 0),wx_hi*wy_lo*wz_hi*pvelz);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   k-1, 0),wx_hi*wy_hi*wz_lo*pvelz);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   k  , 0),wx_hi*wy_hi*wz_hi*pvelz);

            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, k-1, 3),wx_lo*wy_lo*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j-1, k  , 3),wx_lo*wy_lo*wz_hi*pmass);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   k-1, 3),wx_lo*wy_hi*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i-1, j,   k  , 3),wx_lo*wy_hi*wz_hi*pmass);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, k-1, 3),wx_hi*wy_lo*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j-1, k  , 3),wx_hi*wy_lo*wz_hi*pmass);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   k-1, 3),wx_hi*wy_hi*wz_lo*pmass);
            amrex::Gpu::Atomic::Add(&vel_s_arr(i,   j,   k  , 3),wx_hi*wy_hi*wz_hi*pmass);

          });
#endif



void MFIXParticleContainer::
PredictPICVolumeDeposition (int lev,
                            Real dt,
                            RealVect& gravity,
                            EBFArrayBoxFactory* ebfactory,
                            const int ls_refinement,
                            const MultiFab* ls_phi,
                            MultiFab & mf_to_be_filled,
                            const MultiFab * volfrac,
                            const amrex::FabArray<EBCellFlagFab>* flags)
{
  BL_PROFILE("MFIXParticleContainer::SolidsVolumeDeposition()");

  // We always use the coarse dx
  const Geometry& gm  = Geom(0);
  const auto      plo = gm.ProbLoArray();
  const auto      dx  = gm.CellSizeArray();
  const auto      dxi = gm.InvCellSizeArray();

  const auto      reg_cell_vol = dx[0]*dx[1]*dx[2];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_fab_to_be_filled;

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();
      FArrayBox& fab_to_be_filled = mf_to_be_filled[pti];

      const Box& bx  = pti.tilebox(); // I need a box without ghosts


      if ((*flags)[pti].getType(bx) != FabType::covered ) {

        auto volarr = fab_to_be_filled.array();
        const auto& flagsarr = (*flags)[pti].array();
        const auto& vfrac = (*volfrac)[pti].array();

        // Determine if this particle tile actually has any walls
        bool has_wall = false;

        if ((ebfactory != NULL) &&
           ((*flags)[pti].getType(amrex::grow(bx,1)) == FabType::singlevalued))  {

          has_wall = true;

        } else {
          // We need this test for the case of an inflow boundary:
          // inflow does not appear in the EBFactory but
          // the particles see it as a wall

          // Create the nodal refined box based on the current particle tile
          Box refined_box(amrex::convert(amrex::refine(bx,ls_refinement), IntVect{1,1,1}));

          // Set tol to 1/2 dx
          Real tol = amrex::min(dx[0], amrex::min(dx[1], dx[2])) / 2;

          Real ls_min_over_box = ((*ls_phi)[pti]).min<RunOn::Gpu>(refined_box,0);

          if (ls_min_over_box < tol) has_wall = true;

        }

        const auto& phiarr = (has_wall) ? ls_phi->const_array(pti) : Array4<Real const>{};

#ifdef _OPENMP
        const int ncomp = mf_to_be_filled.nComp();
        Box tile_box = bx;

        if(Gpu::notInLaunchRegion())
        {
          tile_box.grow(mf_to_be_filled.nGrow());
          local_fab_to_be_filled.resize(tile_box, ncomp);
          local_fab_to_be_filled.setVal<RunOn::Host>(0.0);
          volarr = local_fab_to_be_filled.array();
        }
#endif

        amrex::ParallelFor(nrp,
          [pstruct,p_realarray,plo,dx,dxi,vfrac,volarr,
           reg_cell_vol,flagsarr, dt, gravity, has_wall,ls_refinement,phiarr]
          AMREX_GPU_DEVICE (int ip) noexcept
          {
            const ParticleType& p = pstruct[ip];

            GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

            // inverse of particle mass
            const Real inv_mass = 1.0/p_realarray[SoArealData::mass][ip];

            // scale wrt drag coefficient and step size
            const amrex::Real scale = 1.0 / (1.0 + dt*p_realarray[SoArealData::dragcoeff][ip]*inv_mass);

            amrex::RealVect pos;

            // updated parcel velocity devoid of the particle normal stress
            pos[0] = p.pos(0) + dt * scale*(p_realarray[SoArealData::velx][ip] +
                        dt*(p_realarray[SoArealData::dragx][ip]*inv_mass + gravity[0]));

            pos[1] = p.pos(1) + dt * scale*(p_realarray[SoArealData::vely][ip] +
                        dt*(p_realarray[SoArealData::dragy][ip]*inv_mass + gravity[1]));

            pos[2] = p.pos(2) + dt * scale*(p_realarray[SoArealData::velz][ip] +
                        dt*(p_realarray[SoArealData::dragz][ip]*inv_mass + gravity[2]));

            if ( has_wall ) {

              const Real radius = p_realarray[SoArealData::radius][ip] *
                  std::cbrt(p_realarray[SoArealData::statwt][ip]);

              Real ls_value = interp_level_set(pos, ls_refinement, phiarr, plo, dxi);

              const Real overlap = radius - ls_value;

              // The particle is actually touching the wall. Reflect it.
              if (overlap > 0.)
              {

                RealVect normal(0.);
                level_set_normal(pos, ls_refinement, normal, phiarr, plo, dxi);

                // Reflect the particle.
                pos[0] += overlap*normal[0];
                pos[1] += overlap*normal[1];
                pos[2] += overlap*normal[2];
              }

            }
            amrex::Real x = (pos[0] - plo[0]) * dxi[0] + 0.5;
            amrex::Real y = (pos[1] - plo[1]) * dxi[1] + 0.5;
            amrex::Real z = (pos[2] - plo[2]) * dxi[2] + 0.5;

            int i = static_cast<int>(amrex::Math::floor(x));
            int j = static_cast<int>(amrex::Math::floor(y));
            int k = static_cast<int>(amrex::Math::floor(z));

            amrex::GpuArray<amrex::Real,2> wx;
            amrex::GpuArray<amrex::Real,2> wy;
            amrex::GpuArray<amrex::Real,2> wz;

            wx[1] = x - i;
            wy[1] = y - j;
            wz[1] = z - k;

            wx[0] = 1.0 - wx[1];
            wy[0] = 1.0 - wy[1];
            wz[0] = 1.0 - wz[1];

            amrex::Real total_weight = 0.0;

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

            for (int ii = 0; ii <= 1; ++ii)
              for (int jj = 0; jj <= 1; ++jj)
                for (int kk = 0; kk <= 1; ++kk)
                  weights[ii][jj][kk] /= total_weight;

            Real pvol = p_realarray[SoArealData::statwt][ip] * p_realarray[SoArealData::volume][ip] / reg_cell_vol;

            for (int kk = -1; kk <= 0; ++kk) {
              for (int jj = -1; jj <= 0; ++jj) {
                for (int ii = -1; ii <= 0; ++ii) {
                  if (flagsarr(i+ii,j+jj,k+kk).isCovered())
                    continue;

                  Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);

                  amrex::Gpu::Atomic::Add(&volarr(i+ii,j+jj,k+kk), weight_vol*pvol);
                }
              }
            }
          });

#ifdef _OPENMP
        if(Gpu::notInLaunchRegion())
        {
          fab_to_be_filled.atomicAdd<RunOn::Host>(local_fab_to_be_filled,
              tile_box, tile_box, 0, 0, ncomp);
        }
#endif

      }
    }
  }
}

