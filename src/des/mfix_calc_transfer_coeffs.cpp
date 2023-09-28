#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_shepard_K.H>
#include <mfix_des_drag_K.H>
#include <mfix_des_conv_coeff_K.H>
#include <mfix_des_usr_reactions_rates_K.H>
#include <mfix_mf_helpers.H>
#include <mfix_algorithm.H>

namespace txfr_aux {

class Accessor
{
  public:
    AMREX_GPU_HOST_DEVICE
    Accessor ()
      : m_gap(-1)
      , m_flag(-1)
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set (const int gap,
              const int flag)
    {
      m_gap = gap;
      m_flag = flag;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator () (const int row,
                     const int column) const
    {
      return row*m_gap + column*m_flag;
    }

  private:
    int m_gap;
    int m_flag;
};

} // end namespace txfr_aux

using namespace txfr_aux;

void mfix::mfix_calc_transfer_coeffs (const Real time,
                                      Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      Vector< MultiFab* > const& X_gk_in,
                                      Vector< Real* > const& pressure_g_in)
{
  const amrex::Real small_num = m_dem.small_number();
  const amrex::Real large_num = m_dem.large_number();
  const amrex::Real eps = m_dem.eps();

  if (m_drag_type == DragType::WenYu) {
    mfix_calc_transfer_coeffs(time, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              ComputeDragWenYu(small_num, large_num, eps));
  }
  else if (m_drag_type == DragType::Gidaspow) {
    mfix_calc_transfer_coeffs(time, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              ComputeDragGidaspow(small_num, large_num, eps));
  }
  else if (m_drag_type == DragType::BVK2) {
    mfix_calc_transfer_coeffs(time, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              ComputeDragBVK2(small_num, large_num, eps));
  }
  else if (m_drag_type == DragType::SyamOBrien) {

    const amrex::Real c1 = m_SyamOBrien_coeff_c1;
    const amrex::Real d1 = m_SyamOBrien_coeff_d1;

    mfix_calc_transfer_coeffs(time, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              ComputeDragSyamOBrien1988(small_num, large_num, eps, c1, d1));
  }
  else if (m_drag_type == DragType::UserDrag) {
    mfix_calc_transfer_coeffs(time, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              ComputeDragUser(small_num, large_num, eps));
  }
  else {
    amrex::Abort("Invalid Drag Type.");
  }
}


template <typename F1>
void mfix::mfix_calc_transfer_coeffs (const Real time,
                                      Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      Vector< MultiFab* > const& X_gk_in,
                                      Vector< Real* > const& pressure_g_in,
                                      F1 DragFunc)
{
  if (fluid.solve_enthalpy())
  {
    if (m_convection_type == ConvectionType::RanzMarshall) {
      mfix_calc_transfer_coeffs(time, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in, DragFunc,
                                ComputeConvRanzMarshall(m_dem.small_number(),m_dem.large_number(),m_dem.eps()));
    }
    else if (m_convection_type == ConvectionType::Gunn) {
      mfix_calc_transfer_coeffs(time, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in, DragFunc,
                                ComputeConvGunn(m_dem.small_number(),m_dem.large_number(),m_dem.eps()));
    }
    else if (m_convection_type == ConvectionType::NullConvection) {
      mfix_calc_transfer_coeffs(time, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in, DragFunc,
                                NullConvectionCoeff());
    }
    else {
      amrex::Abort("Invalid Convection Type.");
    }
  }
  else {
    mfix_calc_transfer_coeffs(time, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in, DragFunc,
                              NullConvectionCoeff());
  }
}


template <typename F1, typename F2>
void mfix::mfix_calc_transfer_coeffs (const Real time,
                                      Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      Vector< MultiFab* > const& X_gk_in,
                                      Vector< Real* > const& pressure_g_in,
                                      F1 DragFunc,
                                      F2 ConvectionCoeff)
{
  if (m_reaction_rates_type == ReactionRatesType::RRatesUser) {
    mfix_calc_transfer_coeffs(time, ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, pressure_g_in,
                              DragFunc, ConvectionCoeff, HeterogeneousRatesUser());
  } else {
    amrex::Abort("Invalid Reaction Rates Type.");
  }
}


template <typename F1, typename F2, typename F3>
void mfix::mfix_calc_transfer_coeffs (const Real time,
                                      Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      Vector< MultiFab* > const& X_gk_in,
                                      Vector< Real* > const& pressure_g_in,
                                      F1 DragFunc,
                                      F2 ConvectionCoeff,
                                      F3 HeterogeneousRRates)
{
  using PairIndex = MFIXParticleContainer::PairIndex;
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  BL_PROFILE("mfix::mfix_calc_transfer_coeff()");

  //***************************************************************************
  // Data for chemical reactions
  //***************************************************************************
  // Fluid species data
  const int nspecies_g = fluid.nspecies();

  // Fluid species data
  const int nspecies_s = solids.nspecies();

  // MFIXReactions data
  const int nreactions = reactions.nreactions();

  // Particles SoA starting indexes for mass fractions and rate of formations
  const int idx_X_sn   = (pc->m_runtimeRealData).X_sn;
  const int idx_mass_txfr = (pc->m_runtimeRealData).mass_txfr;
  const int idx_vel_txfr = (pc->m_runtimeRealData).vel_txfr;
  const int idx_h_txfr = (pc->m_runtimeRealData).h_txfr;

  const auto& fluid_parms = fluid.parameters();
  const auto& solids_parms = solids.parameters();
  const auto& reactions_parms = reactions.parameters();

  //***************************************************************************
  //
  //***************************************************************************

  // We copy the value inside the domain to the outside to avoid
  // unphysical volume fractions.
  const int dir_bc_in = 2;
  m_boundary_conditions.set_epg_bcs(time, ep_g_in, dir_bc_in);

  auto& aux = m_interphase_txfr_deposition->get_aux();
  aux.clear();
  aux.resize(nlev);

  for (int lev = 0; lev < nlev; lev++) {
    aux[lev].clear();

    for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {

      PairIndex index(pti.index(), pti.LocalTileIndex());
      aux[lev][index] = Gpu::DeviceVector<Real>();
    }
  }

  for (int lev = 0; lev < nlev; lev++) {

    // This is just a sanity check to make sure we're not using covered values
    // We can remove these lines once we're confident in the algorithm
    EB_set_covered(*vel_g_in[lev], 0, 3, 1, covered_val);
    EB_set_covered(*ep_g_in[lev], 0, 1, 1, covered_val);
    EB_set_covered(*ro_g_in[lev], 0, 1, 1, covered_val);

    if (fluid.solve_enthalpy()) {
      EB_set_covered(*T_g_in[lev], 0, 1, 1, covered_val);
    }

    if (fluid.solve_species()) {
      EB_set_covered(*X_gk_in[lev], 0, fluid.nspecies(), 1, covered_val);
    }

    bool OnSameGrids = ((dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                        (grids[lev].CellEqual(pc->ParticleBoxArray(lev))));

    MultiFab* interp_ptr(nullptr);
    MultiFab* X_gk_interp_ptr(nullptr);

    const int interp_ng = 1;    // Only one layer needed for interpolation
    const int interp_comp = 5 + // 3 vel_g + 1 ep_g + 1 ro_g
                            1*int(fluid.solve_enthalpy()) +  // 1 T_g
                            1*int(reactions.solve());  //  1 pressure_g

    MultiFab pressure_cc(ep_g_in[lev]->boxArray(), dmap[lev], 1, interp_ng);
    pressure_cc.setVal(0.);

    if (reactions.solve()) {
      MultiFab pressure_nd(m_leveldata[lev]->p0_g->boxArray(), dmap[lev], 1, interp_ng);
      pressure_nd.setVal(0.);
      EB_set_covered(pressure_nd, 0, 1, 1, covered_val);

      MultiFab::Copy(pressure_nd, *m_leveldata[lev]->p0_g, 0, 0, 1, interp_ng);
      pressure_nd.plus(*pressure_g_in[lev], interp_ng);

      pressure_cc.setVal(0.);
      EB_set_covered(pressure_cc, 0, 1, 1, covered_val);

      amrex::average_node_to_cellcenter(pressure_cc, 0, pressure_nd, 0, 1);
    }

    if (OnSameGrids) {

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(grids[lev], dmap[lev], interp_comp, interp_ng, MFInfo(), *ebfactory[lev]);

      int components_count(0);

      // Copy fluid velocity
      MultiFab::Copy(*interp_ptr, *vel_g_in[lev], 0, components_count, 3, interp_ng);
      components_count += 3;

      // Copy volume fraction
      MultiFab::Copy(*interp_ptr, *ep_g_in[lev],  0, components_count, 1, interp_ng);
      components_count += 1;

      // Copy fluid density
      MultiFab::Copy(*interp_ptr, *ro_g_in[lev],  0, components_count, 1, interp_ng);
      components_count += 1;

      if (fluid.solve_enthalpy()) {
        // Copy fluid temperature
        MultiFab::Copy(*interp_ptr, *T_g_in[lev],  0, components_count, 1, interp_ng);
        components_count += 1;
      }

      if (reactions.solve()) {
        // Copy thermodynamic pressure
        MultiFab::Copy(*interp_ptr, pressure_cc,  0, components_count, 1, interp_ng);
        components_count += 1;
      }

      if (fluid.solve_species()) {
        X_gk_interp_ptr = new MultiFab(grids[lev], dmap[lev], nspecies_g, interp_ng, MFInfo(), *ebfactory[lev]);

        // Copy fluid temperature
        MultiFab::Copy(*X_gk_interp_ptr, *X_gk_in[lev], 0, 0, nspecies_g, interp_ng);
      }
    }
    else
    {
      const BoxArray&            pba = pc->ParticleBoxArray(lev);
      const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

      // Temporary arrays  -- copies with no ghost cells
      //const int ng_to_copy = 0; UNUSED VARIABLE

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_comp, interp_ng, MFInfo(), *particle_ebfactory[lev]);

      int components_count(0);

      // Copy fluid velocity
      interp_ptr->ParallelCopy(*vel_g_in[lev], 0, components_count, 3, interp_ng, interp_ng);
      components_count += 3;

      // Copy volume fraction
      interp_ptr->ParallelCopy(*ep_g_in[lev],  0, components_count, 1, interp_ng, interp_ng);
      components_count += 1;

      // Copy fluid density
      interp_ptr->ParallelCopy(*ro_g_in[lev],  0, components_count, 1, interp_ng, interp_ng);
      components_count += 1;

      if (fluid.solve_enthalpy()) {
        // Copy fluid temperature
        interp_ptr->ParallelCopy(*T_g_in[lev],  0, components_count, 1, interp_ng, interp_ng);
        components_count += 1;
      }

      if (reactions.solve()) {
        // Copy fluid species
        interp_ptr->ParallelCopy(pressure_cc,  0, components_count, 1, interp_ng, interp_ng);
        components_count += 1;
      }

      if (fluid.solve_species()) {
        X_gk_interp_ptr = new MultiFab(pba, pdm, nspecies_g, interp_ng, MFInfo(), *particle_ebfactory[lev]);

        // Copy fluid temperature
        X_gk_interp_ptr->ParallelCopy(*X_gk_in[lev], 0, 0, nspecies_g, interp_ng, interp_ng);
      }
    }

    // FillBoundary on interpolation MultiFab
    interp_ptr->FillBoundary(geom[lev].periodicity());

    if (X_gk_interp_ptr != nullptr) {
      X_gk_interp_ptr->FillBoundary(geom[lev].periodicity());
    }

    {
      const auto dxi = geom[lev].InvCellSizeArray();
      const auto dx  = geom[lev].CellSizeArray();
      const auto plo = geom[lev].ProbLoArray();

      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(interp_ptr->Factory());

      const auto cellcent = &(factory.getCentroid());
      const auto bndrycent = &(factory.getBndryCent());
      const auto areafrac = factory.getAreaFrac();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto& ptile = pc->GetParticles(lev)[index];

        //Access to added variables
        auto ptile_data = ptile.getParticleTileData();

        auto& particles = pti.GetArrayOfStructs();
        MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();

        const int np = particles.size();

        aux[lev][index].clear();
        aux[lev][index].resize(fluid.nspecies()*np, 0.);
        Real* aux_ptr = aux[lev][index].dataPtr();

        Box bx = pti.tilebox();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  interp_fab = static_cast<EBFArrayBox const&>((*interp_ptr)[pti]);
        const EBCellFlagFab&  flags = interp_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered) {

          const Array4<const Real> empty_array;

          const auto& interp_array = interp_ptr->const_array(pti);
          const auto& X_gk_array = (X_gk_interp_ptr != nullptr) ? X_gk_interp_ptr->const_array(pti) : empty_array;

          const auto& flags_array  = flags.array();

          const int grown_bx_is_regular = (flags.getType(amrex::grow(bx,1)) == FabType::regular);

          // Cell centroids
          const auto& ccent_fab = grown_bx_is_regular ? empty_array : cellcent->const_array(pti);
          // Centroid of EB
          const auto& bcent_fab = grown_bx_is_regular ? empty_array : bndrycent->const_array(pti);
          // Area fractions
          const auto& apx_fab = grown_bx_is_regular ? empty_array : areafrac[0]->const_array(pti);
          const auto& apy_fab = grown_bx_is_regular ? empty_array : areafrac[1]->const_array(pti);
          const auto& apz_fab = grown_bx_is_regular ? empty_array : areafrac[2]->const_array(pti);

          const int solve_enthalpy = fluid.solve_enthalpy();
          const int fluid_is_a_mixture = fluid.isMixture();
          const int solve_reactions = reactions.solve();

          auto local_cg_dem = m_dem.cg_dem();

          const Real mu_g0 = fluid.mu_g();

          Real* X_gk_ptr = (nspecies_g > 0) ? aux[lev][index].dataPtr() : nullptr;

          Accessor accessor;

#ifdef AMREX_USE_GPU
          int nthreads_per_block = (nreactions > 0) ? 16 : 256;
          int nblocks = std::ceil(Real(np) / nthreads_per_block);

          std::size_t needed_sm_bytes = nreactions*nthreads_per_block*sizeof(Real);
          std::size_t available_sm_per_block = Gpu::Device::sharedMemPerBlock();

          int use_shared_mem(1);
          accessor.set(nthreads_per_block, 1);

          Gpu::DeviceVector<Real> glob_mem;
          Real* glob_mem_ptr(nullptr);

          if (needed_sm_bytes > available_sm_per_block) {
            nthreads_per_block = 256;
            nblocks = std::ceil(Real(np) / nthreads_per_block);

            use_shared_mem = 0;
            needed_sm_bytes = 0;

            if (nreactions > 0) {
              glob_mem.resize(nreactions*np);
              glob_mem_ptr = glob_mem.dataPtr();
            }

            accessor.set(np, 1);
          }

          amrex::launch(nblocks, nthreads_per_block, needed_sm_bytes, Gpu::gpuStream(),
#else
          int use_shared_mem(0);
          Real* glob_mem_ptr(nullptr);

          accessor.set(1, 0);

          amrex::ParallelFor(np,
#endif
              [pstruct,p_realarray,interp_array,DragFunc,ConvectionCoeff,
               HeterogeneousRRates,plo,dxi,solve_enthalpy,fluid_is_a_mixture,
               nspecies_g,interp_comp,local_cg_dem,ptile_data,nreactions,
               nspecies_s,idx_X_sn,idx_mass_txfr,idx_vel_txfr,idx_h_txfr,
               fluid_parms,solids_parms,reactions_parms,flags_array,mu_g0,
               grown_bx_is_regular,dx,ccent_fab,bcent_fab,apx_fab,apy_fab,
               apz_fab,solve_reactions,aux_ptr,np,X_gk_ptr,X_gk_array,accessor,
               use_shared_mem,glob_mem_ptr]
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
            AMREX_GPU_DEVICE (Gpu::Handler const& handler) noexcept
          {
            int p_id = handler.globalIdx();

            Real* R_q_heterogeneous(nullptr);

            if (nreactions > 0) {
              if (use_shared_mem) {
                R_q_heterogeneous = (Real*)handler.sharedMemory();
              } else {
                R_q_heterogeneous = glob_mem_ptr;
              }
            }

            if (p_id < np) {
#else
            AMREX_GPU_DEVICE () noexcept
          {
            int p_id = blockDim.x*blockIdx.x+threadIdx.x;

            Real* R_q_heterogeneous(nullptr);
            Gpu::SharedMemory<Real> gsm;

            if (nreactions > 0) {
              if (use_shared_mem) {
                R_q_heterogeneous = gsm.dataPtr();
              } else {
                R_q_heterogeneous = glob_mem_ptr;
              }
            }

            if (p_id < np) {
#endif
#else
            (int p_id) noexcept
          {
            Vector<Real> R_q_vec(nreactions, 0.);
            Real* R_q_heterogeneous = R_q_vec.dataPtr();

            {
#endif
            MFIXParticleContainer::ParticleType& particle = pstruct[p_id];

            GpuArray<Real, 7> interp_loc; // vel_g, ep_g, ro_g, T_g, p_g
            interp_loc.fill(0.);

            // Indices of cell where particle is located
            const int iloc = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
            const int jloc = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
            const int kloc = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

            if (grown_bx_is_regular) {

              trilinear_interp(particle.pos(), plo, dxi,
                  interp_array, interp_loc.data(), interp_comp,
                  X_gk_array, X_gk_ptr, nspecies_g, np, p_id);

            } else { // FAB not all regular

              // No drag force for particles in covered cells.
              if (flags_array(iloc,jloc,kloc).isCovered()) {

                // drag variable
                p_realarray[SoArealData::dragcoeff][p_id] = 0.;

                // convection-related enthalpy txfr variable
                if (solve_enthalpy) {
                  p_realarray[SoArealData::convection][p_id] = 0.;
                }

                // chemical reaction txfr variables
                if (solve_reactions) {
                  for (int n_g(0); n_g < nspecies_g; n_g++)
                    aux_ptr[n_g*np + p_id] = 0.;

                  ptile_data.m_runtime_rdata[idx_vel_txfr+0][p_id] = 0.;
                  ptile_data.m_runtime_rdata[idx_vel_txfr+1][p_id] = 0.;
                  ptile_data.m_runtime_rdata[idx_vel_txfr+2][p_id] = 0.;

                  // Write the result in the enthalpy transfer space
                  ptile_data.m_runtime_rdata[idx_h_txfr][p_id] = 0.;
                }

                // Nothing else to do. Return
                return;

              // Cut or regular cell and none of the cells in the stencil is
              // covered (Note we can't assume regular cell has no covered
              // cells in the stencil because of the diagonal case)
              } else {

                // Upper cell in trilinear stencil
                int i = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0] + 0.5));
                int j = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1] + 0.5));
                int k = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2] + 0.5));

                // All cells in the stencil are regular. Use
                // traditional trilinear interpolation
                if (flags_array(i-1,j-1,k-1).isRegular() &&
                    flags_array(i  ,j-1,k-1).isRegular() &&
                    flags_array(i-1,j  ,k-1).isRegular() &&
                    flags_array(i  ,j  ,k-1).isRegular() &&
                    flags_array(i-1,j-1,k  ).isRegular() &&
                    flags_array(i  ,j-1,k  ).isRegular() &&
                    flags_array(i-1,j  ,k  ).isRegular() &&
                    flags_array(i  ,j  ,k  ).isRegular()) {

                  trilinear_interp(particle.pos(), plo, dxi,
                                   interp_array, interp_loc.data(), interp_comp,
                                   X_gk_array, X_gk_ptr, nspecies_g, np, p_id);

                // At least one of the cells in the stencil is cut or covered
                } else {
#if 0
                  // TODO: This was initially split for variables that may have known
                  // EB values (e.g., no-slip velocity). However, the results changed
                  // more than expected so now EB values are not used.
                  {
                    const int srccomp = 0;
                    const int dstcomp = 0;
                    const int numcomp = 3;

                    shepard_interp_eb(particle.pos(), iloc, jloc, kloc, dx, dxi, plo,
                                      flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                      interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);

                  }
                  {
                    const int srccomp = 3;
                    const int dstcomp = 3;
                    const int numcomp = interp_comp-3; // ep_g, ro_g, T_g, X_gk

                    shepard_interp(particle.pos(), iloc, jloc, kloc, dx, dxi, plo,
                                   flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                   interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);
                  }
#else
                  const int srccomp = 0;
                  const int dstcomp = 0;
                  const int numcomp = interp_comp; // vel_g, ep_g, ro_g, T_g, p_g

                  shepard_interp(particle.pos(), iloc, jloc, kloc, dx, dxi, plo,
                                 flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                 interp_array, interp_loc.data(), srccomp, dstcomp, numcomp,
                                 X_gk_array, X_gk_ptr, srccomp, dstcomp, nspecies_g, np, p_id);

#endif
                } // Cut cell
              } // Not covered
            } // type of FAB

            RealVect vel_g(interp_loc[0], interp_loc[1], interp_loc[2]);
            const Real ep_g = interp_loc[3];
            const Real ro_g = interp_loc[4];

            int comp_count(5);

            Real T_g(0);
            if (solve_enthalpy) {
              T_g = interp_loc[comp_count];
              comp_count += 1;
            }

            Real mu_g(0);

            if (solve_enthalpy)
              mu_g = fluid_parms.calc_mu_g(T_g);
            else
              mu_g = mu_g0;

            Real rad = p_realarray[SoArealData::radius][p_id];
            Real vol = p_realarray[SoArealData::volume][p_id];

            int pID = particle.id();

            RealVect pvel(0.);
            pvel[0] = p_realarray[SoArealData::velx][p_id];
            pvel[1] = p_realarray[SoArealData::vely][p_id];
            pvel[2] = p_realarray[SoArealData::velz][p_id];

            Real rop_g = ro_g * ep_g;
            RealVect vslp(0.);
            vslp[0] = vel_g[0] - pvel[0];
            vslp[1] = vel_g[1] - pvel[1];
            vslp[2] = vel_g[2] - pvel[2];

            Real vrel = sqrt(dot_product(vslp, vslp));
            Real dp = 2.0*rad;

            if (local_cg_dem) {
               dp = dp/std::cbrt(p_realarray[SoArealData::statwt][p_id]);
            }

            Real ep_s = 1.0 - ep_g;
            Real beta = vol*DragFunc(ep_g, mu_g, rop_g, vrel, dp, dp, ep_s,
                                     vel_g[0], vel_g[1], vel_g[2],
                                     iloc, jloc, kloc, pID);

            p_realarray[SoArealData::dragcoeff][p_id] = beta;

            if(solve_enthalpy) {
              Real k_g = fluid_parms.calc_k_g(T_g);

              Real cp_g(0.);
              if (!fluid_is_a_mixture)
                cp_g = fluid_parms.calc_cp_g<run_on>(T_g);
              else {
                for (int n_g(0); n_g < nspecies_g; ++n_g) {
                  cp_g += X_gk_ptr[n_g*np+p_id]*fluid_parms.calc_cp_gk<run_on>(T_g,n_g);
                }
              }

              Real gamma = ConvectionCoeff(ep_g, mu_g, k_g, cp_g, rop_g, vrel,
                                           dp, iloc, jloc, kloc, pID);

              p_realarray[SoArealData::convection][p_id] = 4.0*M_PI*rad*rad*gamma;
            }

            if (solve_reactions) {
              // Extract interpolated thermodynamic pressure
              const Real p_g = interp_loc[comp_count];

              const Real ro_p = p_realarray[SoArealData::density][p_id];

              int idx = p_id;

#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
              if (use_shared_mem)
                idx = handler.threadIdx();
#else
              if (use_shared_mem)
                idx = threadIdx.x;
#endif
              for (int q(0); q < nreactions; ++q)
                R_q_heterogeneous[accessor(q,idx)] = 0.;
#else
#endif

              const Real T_p = p_realarray[SoArealData::temperature][p_id];
              const Real DP  = 2. * p_realarray[SoArealData::radius][p_id];

              HeterogeneousRRates.template operator()<run_on>(R_q_heterogeneous,
                                                              accessor, idx,
                                                              reactions_parms, np, p_id,
                                                              solids_parms, ptile_data, idx_X_sn,
                                                              ro_p, ep_s, T_p,
                                                              pvel, fluid_parms,
                                                              X_gk_ptr, ro_g, ep_g,
                                                              T_g, vel_g, DP, p_g);

              //***************************************************************
              // Loop over solids species for computing solids txfr rates
              //***************************************************************
              Real G_m_p_heterogeneous(0.);

              for (int n_s(0); n_s < nspecies_s; n_s++) {
                Real G_m_sn_heterogeneous(0.);

                // Loop over reactions to compute each contribution
                for (int q(0); q < nreactions; q++) {

                  const Real stoich_coeff = solids_parms.get_stoich_coeff<run_on>(n_s, q);

                  // Compute solids species n_s transfer rate for reaction q
                  const Real MW_sn = solids_parms.get_MW_sn<run_on>(n_s);
                  Real G_m_sn_q = stoich_coeff * MW_sn * R_q_heterogeneous[accessor(q,idx)];

                  G_m_sn_heterogeneous += G_m_sn_q;
                }

                ptile_data.m_runtime_rdata[idx_mass_txfr+n_s][p_id] = G_m_sn_heterogeneous;

                G_m_p_heterogeneous += G_m_sn_heterogeneous;
              }

              //***************************************************************
              // Loop over fluid species for computing fluid txfr rates
              //***************************************************************
              Real G_m_g_heterogeneous(0.);
              Real G_H_g_heterogeneous(0.);

              for (int n_g(0); n_g < nspecies_g; n_g++) {
                Real G_m_gk_heterogeneous(0.);

                // Loop over reactions to compute each contribution
                for (int q(0); q < nreactions; q++) {

                  Real stoich_coeff = fluid_parms.get_stoich_coeff<run_on>(n_g, q);

                  // Compute fluid species n_g transfer rate for reaction q
                  const Real MW_gk = fluid_parms.get_MW_gk<run_on>(n_g);
                  Real G_m_gk_q = stoich_coeff * MW_gk * R_q_heterogeneous[accessor(q,idx)];

                  G_m_gk_heterogeneous += G_m_gk_q;

                  // Contribution to the particle
                  const Real h_gk_T_p = fluid_parms.calc_h_gk<run_on>(T_p, n_g);
                  const Real h_gk_T_g = fluid_parms.calc_h_gk<run_on>(T_g, n_g);

                  const Real G_H_pk_q = h_gk_T_g * amrex::min(0., G_m_gk_q);
                  const Real G_H_gk_q = h_gk_T_p * amrex::max(0., G_m_gk_q);

                  G_H_g_heterogeneous += (G_H_pk_q + G_H_gk_q);
                }

                aux_ptr[n_g*np+p_id] = G_m_gk_heterogeneous;

                G_m_g_heterogeneous += G_m_gk_heterogeneous;
              }

              AMREX_ASSERT(std::abs(G_m_p_heterogeneous + G_m_g_heterogeneous) < 1.e-15);

              //***************************************************************
              //
              //***************************************************************
              ptile_data.m_runtime_rdata[idx_vel_txfr][p_id] = G_m_p_heterogeneous;

              // Write the result in the enthalpy transfer space
              ptile_data.m_runtime_rdata[idx_h_txfr][p_id] = -G_H_g_heterogeneous;
            }
          }

          }); // pid

          if ((nreactions > 0) && (!use_shared_mem))
            Gpu::synchronize();

        } // if entire FAB not covered
      } // pti
    } // GPU region

    // Free up memory
    delete interp_ptr;

    if (X_gk_interp_ptr != nullptr) {
      delete X_gk_interp_ptr;
    }

  } // lev

  // Reset the volume fractions back to the correct values at
  // inflow faces.
  const int dir_bc_out = 1;
  m_boundary_conditions.set_epg_bcs(time, ep_g_in, dir_bc_out);
}
