#include <restarter.H>
#include <deposition/mfix_deposition_K.H>

#include <AMReX_ParmParse.H>

using namespace amrex;


namespace IdxsAux {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
int get_g_idx (const int i,
               const int j,
               const int k,
               const amrex::IntVect& bx_size,
               const amrex::IntVect& bx_lo)
{
  return (i-bx_lo[0]) + (j-bx_lo[1])*bx_size[0] + (k-bx_lo[2])*bx_size[0]*bx_size[1];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_l_idxs (const int n,
                 int& i,
                 int& j,
                 int& k,
                 const amrex::IntVect& bx_size,
                 const amrex::IntVect& bx_lo)
{
  k = amrex::Math::floor<int>(n / (bx_size[0]*bx_size[1])) + bx_lo[2];
  j = amrex::Math::floor<int>((n - (k-bx_lo[2])*bx_size[0]*bx_size[1]) / bx_size[0]) + bx_lo[1];
  i = n - (k-bx_lo[2])*bx_size[0]*bx_size[1] - (j-bx_lo[1])*bx_size[0] + bx_lo[0];
}

} // end of namespace IdxsAux


MFIXRestarter::MFIXRestarter (const int nlev_in)
  : m_refinement_ratio(1)
  , m_eps_overflow(1.)
  , nlev(nlev_in)
  , avgdPIC_coarse(nlev_in, nullptr)
  , avgdPIC_fine(nlev_in, nullptr)
{
  ParmParse pp("mfix");

  pp.query("refinement_ratio", m_refinement_ratio);
  pp.query("eps_overflow", m_eps_overflow);

  m_PIC_deposition = new MFIXPICDeposition;
}


MFIXRestarter::~MFIXRestarter ()
{
  for (int lev(0); lev < nlev; ++lev) {
    if (avgdPIC_coarse[lev] != nullptr)
      delete avgdPIC_coarse[lev];

    if (avgdPIC_fine[lev] != nullptr)
      delete avgdPIC_fine[lev];
  }

  delete m_PIC_deposition;
}


void
MFIXRestarter::allocate_coarse_arrays (const mfix* mfix_coarse)
{
  auto& solids = mfix_coarse->solids;
  Transfer txfr_idxs(solids);

  for (int lev(0); lev < nlev; ++lev) {
    const auto* epg_coarse = mfix_coarse->m_leveldata[lev]->ep_g;

    avgdPIC_coarse[lev] = new MultiFab(epg_coarse->boxArray(), epg_coarse->DistributionMap(),
        txfr_idxs.count, epg_coarse->nGrow(), MFInfo(), epg_coarse->Factory());
  }
}


void
MFIXRestarter::allocate_fine_arrays (const mfix* mfix_fine)
{
  auto& solids = mfix_fine->solids;
  Transfer txfr_idxs(solids);

  for (int lev(0); lev < nlev; ++lev) {
    const auto* epg_fine = mfix_fine->m_leveldata[lev]->ep_g;

    avgdPIC_fine[lev] = new MultiFab(epg_fine->boxArray(), epg_fine->DistributionMap(),
        txfr_idxs.count, epg_fine->nGrow(), MFInfo(), epg_fine->Factory());
  }
}


void
MFIXRestarter::set_fine_objects (mfix* mfix_fine,
                                 const mfix* mfix_coarse) const
{
  const int nlev = mfix_fine->nlev;

  Vector<Geometry> geom(nlev, Geometry());
  Vector<BoxArray> ba(nlev, BoxArray());

  for (int lev(0); lev < nlev; ++lev) {
    geom[lev] = mfix_coarse->Geom(lev);
    ba[lev] = mfix_coarse->boxArray(lev);

    geom[lev].refine(IntVect(m_refinement_ratio));
    ba[lev].refine(m_refinement_ratio);

    mfix_fine->SetGeometry(lev, geom[lev]);
    mfix_fine->SetBoxArray(lev, ba[lev]);
    mfix_fine->SetDistributionMap(lev, mfix_coarse->DistributionMap(lev));
  }
}


void
MFIXRestarter::txfr_fluid_data (const mfix* mfix_coarse,
                                mfix* mfix_fine) const
{
  Transfer txfr_idxs(mfix_coarse->solids);
  const int txfr_count = txfr_idxs.count;

  const int solve_enthalpy = mfix_coarse->fluid.solve_enthalpy();
  const int solve_species  = mfix_coarse->fluid.solve_species();

  const int nspecies_g = mfix_coarse->fluid.nspecies();

  TxfrAuxiliary aux;

  for (int lev(0); lev < nlev; ++lev) {

    auto& ld_coarse = *(mfix_coarse->m_leveldata[lev]);
    auto& ld_fine   = *(mfix_fine->m_leveldata[lev]);

    const auto& flags_coarse = mfix_coarse->EBFactory(lev).getMultiEBCellFlagFab();
    const auto& flags_fine   = mfix_fine->EBFactory(lev).getMultiEBCellFlagFab();

    const amrex::MultiFab& volfrac_coarse = mfix_coarse->EBFactory(lev).getVolFrac();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*(ld_fine.ep_g),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& box = mfi.tilebox();

      GeometryData geom_coarse_data = mfix_coarse->Geom(lev).data();
      GeometryData geom_fine_data   = mfix_fine->Geom(lev).data();

      Array4<Real> dummy_arr;
      Array4<const Real> dummy_const_arr;

      const auto& ep_g_coarse_arr  = ld_coarse.ep_g->const_array(mfi);
      const auto& ro_g_coarse_arr  = ld_coarse.ro_g->const_array(mfi);
      const auto& trac_coarse_arr  = ld_coarse.trac->const_array(mfi);
      const auto& vel_g_coarse_arr = ld_coarse.vel_g->const_array(mfi);
      const auto& T_g_coarse_arr   = solve_enthalpy ? ld_coarse.T_g->array(mfi) : dummy_const_arr;
      const auto& X_gk_coarse_arr  = solve_species ? ld_coarse.X_gk->array(mfi) : dummy_const_arr;

      const auto& flags_coarse_arr   = flags_coarse.const_array(mfi);
      const auto& volfrac_coarse_arr = volfrac_coarse.const_array(mfi);

      const auto& ep_g_fine_arr  = ld_fine.ep_g->array(mfi);
      const auto& ro_g_fine_arr  = ld_fine.ro_g->array(mfi);
      const auto& trac_fine_arr  = ld_fine.trac->array(mfi);
      const auto& vel_g_fine_arr = ld_fine.vel_g->array(mfi);
      const auto& T_g_fine_arr   = solve_enthalpy ? ld_fine.T_g->array(mfi) : dummy_arr;
      const auto& X_gk_fine_arr  = solve_species ? ld_fine.X_gk->array(mfi) : dummy_arr;

      const auto& flags_fine_arr = flags_fine.const_array(mfi);

      const auto& avgdPIC_coarse_arr = avgdPIC_coarse[lev]->const_array(mfi);
      const auto& avgdPIC_fine_arr = avgdPIC_fine[lev]->array(mfi);

      amrex::ParallelFor(box, [geom_coarse_data,geom_fine_data,ep_g_coarse_arr,
          ro_g_coarse_arr,trac_coarse_arr,vel_g_coarse_arr,T_g_coarse_arr,ep_g_fine_arr,
          ro_g_fine_arr,trac_fine_arr,vel_g_fine_arr,T_g_fine_arr,avgdPIC_coarse_arr,
          avgdPIC_fine_arr,txfr_count,solve_enthalpy,solve_species,X_gk_fine_arr,
          X_gk_coarse_arr,nspecies_g,flags_coarse_arr,flags_fine_arr,
          volfrac_coarse_arr,aux]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        // Create IntVect from fine indexes
        IntVect ijk_fine(i,j,k);

        if (!flags_fine_arr(i,j,k).isCovered()) {

          // Get coordinates corresponding to fine indexes
          RealVect coords = aux.get_coordinates_cc(geom_fine_data, ijk_fine);

          // Get coarse indexes
          IntVect ijk_coarse = aux.get_indexes_cc(geom_coarse_data, coords);

          GpuArray<GpuArray<GpuArray<Real,3>,3>,3> weights;

          const Box& domain_coarse = geom_coarse_data.domain;

          aux.get_weights_cc(weights, volfrac_coarse_arr, ep_g_coarse_arr,
              domain_coarse, ijk_coarse, flags_coarse_arr, ijk_fine, flags_fine_arr);

          ep_g_fine_arr(ijk_fine) = aux.sum_up(weights, domain_coarse,
              flags_coarse_arr, ep_g_coarse_arr, ijk_coarse);

          ro_g_fine_arr(ijk_fine) = aux.sum_up(weights, domain_coarse,
              flags_coarse_arr, ro_g_coarse_arr, ijk_coarse);

          trac_fine_arr(ijk_fine) = aux.sum_up(weights, domain_coarse,
              flags_coarse_arr, trac_coarse_arr, ijk_coarse);

          for (int n(0); n < AMREX_SPACEDIM; ++n) {
            vel_g_fine_arr(ijk_fine,n) = aux.sum_up(weights, domain_coarse,
                flags_coarse_arr, vel_g_coarse_arr, ijk_coarse, n);
          }

          if (solve_enthalpy)
            T_g_fine_arr(ijk_fine) = aux.sum_up(weights, domain_coarse,
                flags_coarse_arr, T_g_coarse_arr, ijk_coarse);

          if (solve_species) {
            Real X_sum(0.);

            for (int n_g(0); n_g < nspecies_g; ++n_g) {
              X_gk_fine_arr(ijk_fine,n_g) = aux.sum_up(weights, domain_coarse,
                  flags_coarse_arr, X_gk_coarse_arr, ijk_coarse, n_g);
              X_sum += X_gk_fine_arr(ijk_fine,n_g);
            }

            AMREX_ALWAYS_ASSERT(X_sum > 1.e-15);

            for (int n_g(0); n_g < nspecies_g; ++n_g) {
              X_gk_fine_arr(ijk_fine,n_g) /= X_sum;
            }
          }

          for (int n(0); n < txfr_count; ++n)
            avgdPIC_fine_arr(ijk_fine,n) = aux.sum_up(weights, domain_coarse,
                flags_coarse_arr, avgdPIC_coarse_arr, ijk_coarse, n);

        } else {

          ep_g_fine_arr(ijk_fine) = mfix::covered_val;
          ro_g_fine_arr(ijk_fine) = mfix::covered_val;
          trac_fine_arr(ijk_fine) = mfix::covered_val;

          for (int n(0); n < AMREX_SPACEDIM; ++n) {
            vel_g_fine_arr(ijk_fine,n) = mfix::covered_val;
          }

          if (solve_enthalpy)
            T_g_fine_arr(ijk_fine) = mfix::covered_val;

          if (solve_species)
            for (int n_g(0); n_g < nspecies_g; ++n_g)
              X_gk_fine_arr(ijk_fine,n_g) = mfix::covered_val;

          for (int n(0); n < txfr_count; ++n)
            avgdPIC_fine_arr(ijk_fine,n) = mfix::covered_val;
        }
      });
    }

    EB_set_covered(*ld_fine.ep_g,  0, ld_fine.ep_g->nComp(),  ld_fine.ep_g->nGrow(),  mfix::covered_val);
    EB_set_covered(*ld_fine.ro_g,  0, ld_fine.ro_g->nComp(),  ld_fine.ro_g->nGrow(),  mfix::covered_val);
    EB_set_covered(*ld_fine.trac,  0, ld_fine.trac->nComp(),  ld_fine.trac->nGrow(),  mfix::covered_val);
    EB_set_covered(*ld_fine.vel_g, 0, ld_fine.vel_g->nComp(), ld_fine.vel_g->nGrow(), mfix::covered_val);

    if (solve_enthalpy)
      EB_set_covered(*ld_fine.T_g,   0, ld_fine.T_g->nComp(),   ld_fine.T_g->nGrow(),   mfix::covered_val);

    if (solve_species)
      EB_set_covered(*ld_fine.X_gk,  0, ld_fine.X_gk->nComp(),  ld_fine.X_gk->nGrow(),  mfix::covered_val);

    EB_set_covered(*avgdPIC_fine[lev], 0, avgdPIC_fine[lev]->nComp(), avgdPIC_fine[lev]->nGrow(), mfix::covered_val);
  }
}


void
MFIXRestarter::calc_txfr (const mfix* mfix_coarse,
                          const Vector<MultiFab*>& avgdPIC,
                          const Real time)
{
  for (int lev = 0; lev < nlev; lev++) {
    avgdPIC[lev]->setVal(0);
  }

  m_PIC_deposition->set_eps_size(nlev);

  mfix_coarse->mfix_deposit_particles(m_PIC_deposition, avgdPIC, time);

  Vector< MultiFab* > const& tmp_eps = m_PIC_deposition->get_eps();

  // Divide
  for (int lev = 0; lev < nlev; lev++) {

    const auto& solids = mfix_coarse->solids;

    Transfer txfr_idxs(solids);
    const int idx_eps = txfr_idxs.idx_eps;
    const int idx_radius = txfr_idxs.idx_radius;
    const int idx_density = txfr_idxs.idx_density;
    const int idx_vel = txfr_idxs.idx_vel;
    const int idx_temp = txfr_idxs.idx_temp;
    const int idx_species = txfr_idxs.idx_species;

    const int nspecies_s = solids.nspecies();

    EB_set_covered(*avgdPIC[lev], mfix::covered_val);
//    EB_set_covered(*avgdPIC[lev], 0.);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*avgdPIC[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& box = mfi.tilebox();

      const auto& avgdPIC_arr = avgdPIC[lev]->array(mfi);

      const int solve_enthalpy = solids.solve_enthalpy();
      const int solve_species  = solids.solve_species();

      amrex::ParallelFor(box, [avgdPIC_arr,solve_enthalpy,solve_species,idx_eps,
          idx_vel,idx_species,idx_temp,nspecies_s,idx_radius,idx_density]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real eps = avgdPIC_arr(i,j,k,idx_eps);

//        if (eps < 0 || eps >= 1.) {
//          printf("i,j,k = (%d,%d,%d), eps = %e\n", i,j,k,eps);
//          amrex::Abort("Error 11");
//        }

        if ((eps > 1.e-15) && (avgdPIC_arr(i,j,k,idx_radius) < 1.e-15)) {
          printf("i,j,k = (%d,%d,%d), eps = %e, pradius = %e\n", i,j,k,eps,
            avgdPIC_arr(i,j,k,idx_radius));
          amrex::Abort("Error 22");
        }

        if ((eps > 1.e-15) && (avgdPIC_arr(i,j,k,idx_density) < 1.e-15)) {
          printf("i,j,k = (%d,%d,%d), eps = %e, pdensity = %e\n", i,j,k,eps,
            avgdPIC_arr(i,j,k,idx_density));
          amrex::Abort("Error 22");
        }

        if (eps < 1.e-15) {

          avgdPIC_arr(i,j,k,idx_radius) = 0.;
          avgdPIC_arr(i,j,k,idx_density) = 0.;
          avgdPIC_arr(i,j,k,idx_vel+0) = 0.;
          avgdPIC_arr(i,j,k,idx_vel+1) = 0.;
          avgdPIC_arr(i,j,k,idx_vel+2) = 0.;

          if (solve_enthalpy)
            avgdPIC_arr(i,j,k,idx_temp) = 0.;

          if (solve_species)
            for (int n_s(0); n_s < nspecies_s; ++n_s)
              avgdPIC_arr(i,j,k,idx_species+n_s) = 0.;

        } else {

          avgdPIC_arr(i,j,k,idx_radius) /= eps;
          avgdPIC_arr(i,j,k,idx_density) /= eps;
          avgdPIC_arr(i,j,k,idx_vel+0) /= eps;
          avgdPIC_arr(i,j,k,idx_vel+1) /= eps;
          avgdPIC_arr(i,j,k,idx_vel+2) /= eps;

          if (solve_enthalpy)
            avgdPIC_arr(i,j,k,idx_temp) /= eps;

          if (solve_species)
            for (int n_s(0); n_s < nspecies_s; ++n_s)
              avgdPIC_arr(i,j,k,idx_species+n_s) /= eps;

          // TODO TODO TODO
          // make sure X_sum = 1
        }
      });
    }

    EB_set_covered(*avgdPIC[lev], mfix::covered_val);
  }

  m_PIC_deposition->clear_eps();
}


void
MFIXPICDeposition::deposit (int lev,
                            const Geometry& geom,
                            MFIXParticleContainer* pc,
                            const MultiFab* volfrac,
                            const amrex::FabArray<EBCellFlagFab>* flags,
                            MultiFab* txfr_mf,
                            MultiFab* eps_mf)
{
  if (mfix::m_deposition_scheme == DepositionScheme::trilinear) {

    deposit(TrilinearDeposition(), lev, geom, pc, volfrac, flags, txfr_mf, eps_mf);

  } else if (mfix::m_deposition_scheme == DepositionScheme::square_dpvm) {

    deposit(TrilinearDPVMSquareDeposition(), lev, geom, pc, volfrac, flags, txfr_mf, eps_mf);

  } else if (mfix::m_deposition_scheme == DepositionScheme::true_dpvm) {

    deposit(TrueDPVMDeposition(), lev, geom, pc, volfrac, flags, txfr_mf, eps_mf);

  } else if (mfix::m_deposition_scheme == DepositionScheme::centroid) {

    deposit(CentroidDeposition(), lev, geom, pc, volfrac, flags, txfr_mf, eps_mf);

  } else {

    amrex::Abort("Don't know this deposition_scheme!");
  }
}


template <typename F>
void
MFIXPICDeposition::deposit (F WeightFunc,
                            int lev,
                            const Geometry& geom,
                            MFIXParticleContainer* pc,
                            const MultiFab* volfrac,
                            const amrex::FabArray<EBCellFlagFab>* flags,
                            MultiFab* txfr_mf,
                            MultiFab* eps_mf)
{
  // We always use the coarse dx
  const auto plo = geom.ProbLoArray();
  const auto dx  = geom.CellSizeArray();
  const auto dxi = geom.InvCellSizeArray();
  const Real dV = dx[0]*dx[1]*dx[2];

  const auto& solids = pc->get_solids();

  const int nspecies_s = solids.nspecies();

  Transfer txfr_idxs(solids);
  const int idx_eps     = txfr_idxs.idx_eps;
  const int idx_radius  = txfr_idxs.idx_radius;
  const int idx_density = txfr_idxs.idx_density;
  const int idx_vel     = txfr_idxs.idx_vel;
  const int idx_temp    = txfr_idxs.idx_temp;
  const int idx_species = txfr_idxs.idx_species;

  const int idx_X_sn = pc->m_runtimeRealData.X_sn;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox local_eps;
    FArrayBox local_txfr;

    for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {

      PairIndex index(pti.index(), pti.LocalTileIndex());
      auto& ptile = pc->GetParticles(lev)[index];

      //Access to added variables
      auto ptile_data = ptile.getParticleTileData();

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();
      
      const long nrp = pti.numParticles();

      FArrayBox dummy_fab;

      FArrayBox& eps_fab = (*eps_mf)[pti];
      FArrayBox& txfr_fab = (*txfr_mf)[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered) {

        const auto& eps_arr = eps_fab.array();
        const auto& txfr_arr = txfr_fab.array();
        const auto& flags_arr = flags->const_array(pti);
        const auto& vfrac_arr = volfrac->const_array(pti);

        const Real deposition_scale_factor = mfix::m_deposition_scale_factor;

#ifdef _OPENMP
        Box eps_tile_box = box;
        {
          const int ncomp = eps_mf->nComp();

          if (Gpu::notInLaunchRegion()) {
            eps_tile_box.grow(eps_mf->nGrow());
            local_eps.resize(eps_tile_box, ncomp);
            local_eps.setVal<RunOn::Host>(0.0);
            eps_arr = local_eps.array();
          }
        }

        Box txfr_tile_box = box;
        {
          const int ncomp = txfr_mf->nComp();

          if (Gpu::notInLaunchRegion()) {
            txfr_tile_box.grow(txfr_mf->nGrow());
            local_txfr.resize(txfr_tile_box, ncomp);
            local_txfr.setVal<RunOn::Host>(0.0);
            txfr_arr = local_txfr.array();
          }
        }
#endif

        const int solve_enthalpy = solids.solve_enthalpy();
        const int solve_species  = solids.solve_species();

        amrex::ParallelFor(nrp, [pstruct,p_realarray,plo,dx,dxi,ptile_data,dV,idx_eps,
            deposition_scale_factor,WeightFunc,flags_arr,txfr_arr,solve_enthalpy,
            idx_vel,idx_temp,idx_species,nspecies_s,idx_X_sn,solve_species,vfrac_arr,
            idx_radius,idx_density,eps_arr]
          AMREX_GPU_DEVICE (int ip) noexcept
        {
          const ParticleType& p = pstruct[ip];

          int i; int j; int k;

          const Real pradius = p_realarray[SoArealData::radius][ip];
          const Real pdensity = p_realarray[SoArealData::density][ip];

          const Real statwt  = p_realarray[SoArealData::statwt][ip];

          GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

          WeightFunc(plo, dx, dxi, flags_arr, p.pos(), pradius, i, j, k, weights,
            deposition_scale_factor);

          const Real pvol = p_realarray[SoArealData::volume][ip];

          const Real pvel_x = p_realarray[SoArealData::velx][ip];
          const Real pvel_y = p_realarray[SoArealData::vely][ip];
          const Real pvel_z = p_realarray[SoArealData::velz][ip];

          Real ptemp(0.);
          if (solve_enthalpy)
            ptemp = p_realarray[SoArealData::temperature][ip];

          GpuArray<Real, MFIXSpecies::NMAX> pspecies;
          pspecies.fill(0);

          if (solve_species) {
            for (int n_s(0); n_s < nspecies_s; ++n_s) {
              pspecies[n_s] = ptile_data.m_runtime_rdata[idx_X_sn+n_s][ip];
            }
          }

          // Deposition
          for (int ii = -1; ii <= 0; ++ii) {
            for (int jj = -1; jj <= 0; ++jj) {
              for (int kk = -1; kk <= 0; ++kk) {
                if (flags_arr(i+ii,j+jj,k+kk).isCovered())
                  continue;

                Real weight = weights[ii+1][jj+1][kk+1] * (pvol/(vfrac_arr(i+ii,j+jj,k+kk)*dV));

                // TODO improve the following
                HostDevice::Atomic::Add(&eps_arr(i+ii,j+jj,k+kk), statwt*weight);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_eps), statwt*weight);

                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_radius), statwt*weight*pradius);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_density), statwt*weight*pdensity);

                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel+0), statwt*weight*pvel_x);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel+1), statwt*weight*pvel_y);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel+2), statwt*weight*pvel_z);

                if (solve_enthalpy) {
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_temp), statwt*weight*ptemp);
                }

                if (solve_species) {
                  for (int n_s(0); n_s < nspecies_s; ++n_s) {
                    HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_species+n_s), statwt*weight*pspecies[n_s]);
                  }
                }
              }
            }
          }
        });

#ifdef _OPENMP
        if (Gpu::notInLaunchRegion()) {
          const int ncomp = eps_mf->nComp();
          eps_fab.atomicAdd<RunOn::Host>(local_eps, eps_tile_box, eps_tile_box, 0, 0, ncomp);
        }

        if (Gpu::notInLaunchRegion()) {
          const int ncomp = txfr_mf->nComp();
          txfr_fab.atomicAdd<RunOn::Host>(local_txfr, txfr_tile_box, txfr_tile_box, 0, 0, ncomp);
        }
#endif
      }
    }
  }

  if (m_eps[lev] != nullptr)
    delete m_eps[lev];

  m_eps[lev] = new MultiFab(eps_mf->boxArray(), eps_mf->DistributionMap(),
                            eps_mf->nComp(), eps_mf->nGrow(), MFInfo(),
                            eps_mf->Factory());
  
  MultiFab::Copy(*m_eps[lev], *eps_mf, 0, 0, eps_mf->nComp(), eps_mf->nGrow());
}


void
MFIXRestarter::generate_particles (const mfix* mfix_coarse,
                                   mfix* mfix_fine) const
{
  // Allocate the particle data
  if (mfix_fine->m_dem.solve())
  {
    MFIXParticleContainer* pc_fine = mfix_fine->pc;

    Real strt_init_part = ParallelDescriptor::second();

    pc_fine->AllocData();

    amrex::Print() << "Auto generating particles ..." << std::endl;

    init_particles(mfix_coarse, pc_fine);

    pc_fine->Redistribute();

    Real end_init_part = ParallelDescriptor::second() - strt_init_part;
    ParallelDescriptor::ReduceRealMax(end_init_part, ParallelDescriptor::IOProcessorNumber());
    Print() << "Time spent in initializing particles " << end_init_part << std::endl;
  }
}


void
MFIXRestarter::init_particles (const mfix* mfix_coarse,
                               MFIXParticleContainer* pc_fine) const
{
  int lev = 0;

  const auto& solids = mfix_coarse->solids;

  Transfer txfr_idxs(solids);
  const int idx_eps = txfr_idxs.idx_eps;
  const int idx_radius = txfr_idxs.idx_radius;

  MultiFab* eps_ptr;
  MultiFab* pradius_ptr;

  MultiFab* epg_ptr;

  const auto& geom_coarse = mfix_coarse->Geom(lev);

  const auto& f_ba = mfix_coarse->boxArray(lev);
  const auto& f_dmap = mfix_coarse->DistributionMap(lev);

  const BoxArray& p_ba = mfix_coarse->ParticleEBFactory(lev).boxArray();
  const DistributionMapping& p_dmap = mfix_coarse->ParticleEBFactory(lev).DistributionMap();

  bool OnSameGrids = ((f_dmap == p_dmap) && (f_ba.CellEqual(p_ba)));

  if (OnSameGrids) {

    eps_ptr = new MultiFab(f_ba, f_dmap, 1, 0, MFInfo());
    eps_ptr->setVal(0.);
    MultiFab::Copy(*eps_ptr, *avgdPIC_coarse[lev], idx_eps, 0, 1, 0);

    pradius_ptr = new MultiFab(f_ba, f_dmap, 1, 0, MFInfo());
    pradius_ptr->setVal(0.);
    MultiFab::Copy(*pradius_ptr, *avgdPIC_coarse[lev], idx_radius, 0, 1, 0);

    epg_ptr = mfix_coarse->m_leveldata[lev]->ep_g;

  } else {

    eps_ptr = new MultiFab(p_ba, p_dmap, 1, 0);
    eps_ptr->setVal(0.);
    eps_ptr->ParallelCopy(*avgdPIC_coarse[lev], idx_eps, 0, 1, 0, 0);

    pradius_ptr = new MultiFab(p_ba, p_dmap, 1, 0);
    pradius_ptr->setVal(0.);
    pradius_ptr->ParallelCopy(*avgdPIC_coarse[lev], idx_radius, 0, 1, 0, 0);

    epg_ptr = new MultiFab(p_ba, p_dmap, 1, 0);
    epg_ptr->setVal(0.);
    epg_ptr->ParallelCopy(*(mfix_coarse->m_leveldata[lev]->ep_g), 0, 0, 1, 0, 0);

  }

  eps_ptr->FillBoundary(geom_coarse.periodicity());
  EB_set_covered(*eps_ptr, mfix::covered_val);

  pradius_ptr->FillBoundary(geom_coarse.periodicity());
  EB_set_covered(*pradius_ptr, mfix::covered_val);

  epg_ptr->FillBoundary(geom_coarse.periodicity());
  EB_set_covered(*epg_ptr, mfix::covered_val);

//  {
//    std::ofstream file;
//    file.open("eps.fab");
//    (*eps_ptr)[0].writeOn(file, 0, 1);
//    file.close();
//  }

//  const RealVect dxi(geom_coarse.InvCellSize());
  const RealVect dx(geom_coarse.CellSize());
  const RealVect plo(geom_coarse.ProbLo());

  long total_np = 0;

  const auto& flags = mfix_coarse->ParticleEBFactory(lev).getMultiEBCellFlagFab();

  using PairIndex = MFIXParticleContainer::PairIndex;

  std::map<PairIndex, Gpu::DeviceVector<Hex_ClosePack>> hcp_vector;
  std::map<PairIndex, Gpu::DeviceVector<Long>> gen_indexes;
  std::map<PairIndex, Gpu::DeviceVector<Long>> gen_number;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(*epg_ptr,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    const IntVect bx_smallEnd = bx.smallEnd();
    const IntVect bx_size = bx.size();
    const int bx_numPts = bx.numPts();

    AMREX_ALWAYS_ASSERT(bx_numPts == (bx_size[0]*bx_size[1]*bx_size[2]));

    PairIndex index(mfi.index(), mfi.LocalTileIndex());

    hcp_vector[index].resize(bx_numPts, Hex_ClosePack());
    gen_indexes[index].resize(bx_numPts, 0);
    gen_number[index].resize(bx_numPts, 0);

    EBCellFlagFab const& flags_fab = flags[mfi];

    if (flags_fab.getType(bx) != FabType::covered) {

      Hex_ClosePack* hcp_vector_ptr = hcp_vector[index].dataPtr();
      Long* gen_indexes_ptr = gen_indexes[index].dataPtr();
      Long* gen_number_ptr = gen_number[index].dataPtr();

      Array4<const Real> const& eps_arr = eps_ptr->const_array(mfi);
      Array4<const Real> const& pradius_arr = pradius_ptr->const_array(mfi);

      Array4<const Real> const& epg_arr = epg_ptr->const_array(mfi);

      Array4<EBCellFlag const> const& flags_arr = flags.const_array(mfi);

      const amrex::Real eps_overflow = m_eps_overflow;

      // compute nb of particles for each cell
      amrex::ParallelFor(bx, [plo,dx,bx_smallEnd,bx_size,eps_arr,flags_arr,
          hcp_vector_ptr,pradius_arr,epg_arr,gen_number_ptr,eps_overflow]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        amrex::Real eps = eps_arr(i,j,k)*eps_overflow;

        if (!flags_arr(i,j,k).isCovered() && (eps > 1.e-15 && eps < 1.)) {

          amrex::RealVect dlo(plo[0] + dx[0]*i, plo[1] + dx[1]*j, plo[2] + dx[2]*k);
          amrex::RealVect dhi(plo[0] + dx[0]*(i+1), plo[1] + dx[1]*(j+1), plo[2] + dx[2]*(k+1));

          const Long n = IdxsAux::get_g_idx(i, j, k, bx_size, bx_smallEnd);

          hcp_vector_ptr[n].initialize(plo, dx);

          RealBox region(dlo.begin(), dhi.begin());

          Real diameter = 2*pradius_arr(i,j,k);

          hcp_vector_ptr[n].setup({i,j,k}, {i,j,k}, region, diameter, eps);

          gen_number_ptr[n] = hcp_vector_ptr[n].get_particles_number();
        }
      });
    }
  }

  // No OMP Parallel on this loop!
  for (MFIter mfi(*epg_ptr,false); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    PairIndex index(mfi.index(), mfi.LocalTileIndex());

    const int vec_size = gen_indexes[index].size();

    AMREX_ASSERT(vec_size == bx.numPts());

    BL_PROFILE_VAR("Inclusive_scan", inclusive_scan);

    Gpu::inclusive_scan(gen_number[index].begin(), gen_number[index].end(),
        gen_indexes[index].begin());

    Gpu::synchronize();

    BL_PROFILE_VAR_STOP(inclusive_scan);

    Long local_np(0);

#ifdef AMREX_USE_GPU
    Gpu::HostVector<Long> gen_indexes_host(vec_size, 0);
    Gpu::copy(Gpu::deviceToHost, gen_indexes[index].begin(), gen_indexes[index].end(),
        gen_indexes_host.begin());

    local_np = gen_indexes_host[vec_size-1];
#else
    local_np = gen_indexes[index][vec_size-1];
#endif

    Long* gen_indexes_ptr = gen_indexes[index].dataPtr();

    ParallelFor(vec_size, [gen_indexes_ptr,total_np]
      AMREX_GPU_DEVICE (int n) noexcept
    {
      gen_indexes_ptr[n] += total_np;
    });

    total_np += local_np;

    // Now that we know pcount, go ahead and create a particle container for this
    // grid and add the particles to it
    auto& particles = pc_fine->DefineAndReturnParticleTile(lev, mfi);

    const int current_size = particles.numParticles();
    particles.resize(current_size + local_np);

    // Add components for each of the runtime variables
    const int start = SoArealData::count;
    for (int comp(0); comp < pc_fine->m_runtimeRealData.count; ++comp)
      particles.push_back_real(start+comp, local_np, 0.);

  } // MFIter loop

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(*eps_ptr,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    PairIndex index(mfi.index(), mfi.LocalTileIndex());

    Hex_ClosePack* hcp_vector_ptr = hcp_vector[index].dataPtr();
    Long* gen_number_ptr = gen_number[index].dataPtr();
    Long* gen_indexes_ptr = gen_indexes[index].dataPtr();

    // Now that we know pcount, go ahead and create a particle container for this
    // grid and add the particles to it
    auto& particles = pc_fine->ParticlesAt(lev, mfi);

    const int id = MFIXParticleContainer::ParticleType::NextID();
    const int cpu = ParallelDescriptor::MyProc();

    const int idx_species = pc_fine->m_runtimeRealData.X_sn;

    // generate positions
    generate_particles(total_np, bx, particles, hcp_vector_ptr, gen_number_ptr,
                       gen_indexes_ptr, id, cpu, mfix_coarse->solids, idx_species);

    // TODO TODO TODO fix this
    // Update the particles NextID
    MFIXParticleContainer::ParticleType::NextID(id+total_np);
  }

  delete eps_ptr;

  ParallelDescriptor::ReduceLongSum(total_np);

  amrex::Print() << "Total number of generated particles: " << total_np << std::endl;
  pc_fine->setTotalNumParticles(total_np);

  // We shouldn't need this if the particles are tiled with one tile per grid, but otherwise
  // we do need this to move particles from tile 0 to the correct tile.
  pc_fine->Redistribute();
}


void
MFIXRestarter::generate_particles (const Long particles_count,
                                   const Box& bx,
                                   MFIXParticleContainer::ParticleTileType& particles,
                                   const Hex_ClosePack* hcp_vector_ptr,
                                   const Long* gen_number_ptr,
                                   const Long* gen_indexes_ptr,
                                   const int id,
                                   const int cpu,
                                   const MFIXSolidsPhase& solids,
                                   const int idx_species) const
{
  auto ptile_data = particles.getParticleTileData();

  auto& aos = particles.GetArrayOfStructs();
  MFIXParticleContainer::ParticleType* pstruct = aos().dataPtr();

  auto& soa = particles.GetStructOfArrays();
  auto p_realarray = soa.realarray();
  auto p_intarray = soa.intarray();

  const int nspecies_s = solids.nspecies();

  const IntVect bx_size = bx.size();
  const IntVect bx_lo = bx.smallEnd();

  amrex::ParallelForRNG(bx, [pstruct,p_realarray,p_intarray,id,cpu,nspecies_s,
      idx_species,ptile_data,hcp_vector_ptr,gen_number_ptr,gen_indexes_ptr,bx_lo,
      bx_size]
    AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
  {
    const Long n = IdxsAux::get_g_idx(i, j, k, bx_size, bx_lo);

    const int np = gen_number_ptr[n];

    const int glob_id = gen_indexes_ptr[n] - np;

    for (int nn = 0; nn < np; ++nn) {

      const Long p_tot = glob_id + nn;

      MFIXParticleContainer::ParticleType& part = pstruct[p_tot];

      part.id() = id + p_tot;
      part.cpu() = cpu;

      RealVect position = hcp_vector_ptr[n].template operator()<run_on>(nn, engine);

      part.pos(0) = position[0];
      part.pos(1) = position[1];
      part.pos(2) = position[2];

      p_realarray[SoArealData::velx][p_tot] = 9.87654321e32;
      p_realarray[SoArealData::vely][p_tot] = 9.87654321e32;
      p_realarray[SoArealData::velz][p_tot] = 9.87654321e32;

      p_realarray[SoArealData::statwt][p_tot] = 1;

      p_realarray[SoArealData::radius][p_tot] = 0.0;
      p_realarray[SoArealData::density][p_tot] = 0.0;

      p_realarray[SoArealData::volume][p_tot] = 0.0;
      p_realarray[SoArealData::mass][p_tot] = 0.0;
      p_realarray[SoArealData::oneOverI][p_tot] = 0.0;

      p_realarray[SoArealData::omegax][p_tot] = 0.0;
      p_realarray[SoArealData::omegay][p_tot] = 0.0;
      p_realarray[SoArealData::omegaz][p_tot] = 0.0;

      p_realarray[SoArealData::dragcoeff][p_tot] = 0.0;

      p_realarray[SoArealData::dragx][p_tot] = 0.0;
      p_realarray[SoArealData::dragy][p_tot] = 0.0;
      p_realarray[SoArealData::dragz][p_tot] = 0.0;

      for (int n_s(0); n_s < nspecies_s; ++n_s) {
        ptile_data.m_runtime_rdata[idx_species+n_s][p_tot] = 9.87654321e32;
      }

      p_intarray[SoAintData::phase][p_tot] = 1;
      p_intarray[SoAintData::state][p_tot] = 1;
    }
  });

  return;
}


void
MFIXRestarter::init_particles_data (mfix* mfix_fine) const
{
  MFIXParticleContainer* pc_fine = mfix_fine->pc;
  const auto& solids = mfix_fine->solids;

  const int solve_species = solids.solve_species();
  const int solve_enthalpy = solids.solve_enthalpy();

  const int nspecies_s = solids.nspecies();
  const int idx_species = pc_fine->m_runtimeRealData.X_sn;

  Transfer fld_transfer(solids);

  for (int lev = 0; lev < nlev; lev++) {

    const auto& geom = mfix_fine->Geom(lev);

    const auto& f_ba = mfix_fine->boxArray(lev);
    const auto& f_dmap = mfix_fine->DistributionMap(lev);

    const BoxArray& p_ba = mfix_fine->ParticleEBFactory(lev).boxArray();
    const DistributionMapping& p_dmap = mfix_fine->ParticleEBFactory(lev).DistributionMap();

    MultiFab* avgdPIC_ptr(nullptr);

    bool OnSameGrids = ((f_dmap == p_dmap) && (f_ba.CellEqual(p_ba)));

    if (OnSameGrids) {

      avgdPIC_ptr = avgdPIC_fine[lev];

    } else {

      avgdPIC_ptr = new MultiFab(p_ba, p_dmap, avgdPIC_fine[lev]->nComp(), 0);
      avgdPIC_ptr->setVal(0.);

      avgdPIC_ptr->ParallelCopy(*avgdPIC_fine[lev], 0, 0, avgdPIC_fine[lev]->nComp(), 0, 0);
    }

    avgdPIC_ptr->FillBoundary(geom.periodicity());

    const auto dxi = geom.InvCellSizeArray();
    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIXParticleContainer::MFIXParIter pti(*pc_fine,lev); pti.isValid(); ++pti) {

      MFIXParticleContainer::PairIndex index(pti.index(), pti.LocalTileIndex());
      auto& ptile = pc_fine->GetParticles(lev)[index];
      auto ptile_data = ptile.getParticleTileData();

      auto& particles = pti.GetArrayOfStructs();
      MFIXParticleContainer::ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const int np = particles.size();

      const auto& avgdPIC_array = avgdPIC_ptr->const_array(pti);

      ParallelFor(np, [pstruct,p_realarray,ptile_data,solve_species,plo,
          avgdPIC_array,nspecies_s,idx_species,fld_transfer,dxi,solve_enthalpy]
        AMREX_GPU_DEVICE (int p) noexcept
      {
        MFIXParticleContainer::ParticleType& part = pstruct[p];

        const RealVect& pos = part.pos();

        const RealVect lx((pos[0] - plo[0])*dxi[0],
                          (pos[1] - plo[1])*dxi[1],
                          (pos[2] - plo[2])*dxi[2]);

        const IntVect ijk = lx.floor();

        int i = ijk[0]; int j = ijk[1]; int k = ijk[2];

        p_realarray[SoArealData::radius][p] = avgdPIC_array(i,j,k,fld_transfer.idx_radius);
        p_realarray[SoArealData::density][p] = avgdPIC_array(i,j,k,fld_transfer.idx_density);

        const Real radius = p_realarray[SoArealData::radius][p];
        const Real volume = (4./3.)*M_PI*(radius*radius*radius);
        const Real density = p_realarray[SoArealData::density][p];
        const Real mass = volume*density;
        const Real oneOverI = 10. / (mass*4.*radius*radius);

        p_realarray[SoArealData::volume][p] = volume;
        p_realarray[SoArealData::mass][p] = mass;
        p_realarray[SoArealData::oneOverI][p] = oneOverI;

        p_realarray[SoArealData::velx][p] = avgdPIC_array(i,j,k,fld_transfer.idx_vel+0);
        p_realarray[SoArealData::vely][p] = avgdPIC_array(i,j,k,fld_transfer.idx_vel+1);
        p_realarray[SoArealData::velz][p] = avgdPIC_array(i,j,k,fld_transfer.idx_vel+2);

        // Add temperature and species initialization
        if (solve_enthalpy)
          p_realarray[SoArealData::temperature][p] = avgdPIC_array(i,j,k,fld_transfer.idx_temp);

        if (solve_species) {
          Real X_sum(0);

          for (int n_s(0); n_s < nspecies_s; ++n_s) {
            ptile_data.m_runtime_rdata[idx_species+n_s][p] = avgdPIC_array(i,j,k,fld_transfer.idx_species+n_s);
            X_sum += ptile_data.m_runtime_rdata[idx_species+n_s][p];
          }

          AMREX_ALWAYS_ASSERT(X_sum > 1.e-15);

          for (int n_s(0); n_s < nspecies_s; ++n_s) {
            ptile_data.m_runtime_rdata[idx_species+n_s][p] /= X_sum;
          }
        }
      });
    } // ParIter loop
  } // lev loop

  return;
}
