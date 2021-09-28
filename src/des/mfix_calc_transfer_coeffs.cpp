#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_shepard_K.H>
#include <mfix_des_drag_K.H>
#include <mfix_des_conv_coeff_K.H>
#include <mfix_mf_helpers.H>

void mfix::mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      Vector< MultiFab* > const& X_gk_in)
{
  if (m_drag_type == DragType::WenYu) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in,
                              ComputeDragWenYu(DEM::small_number, DEM::large_number, DEM::eps));
  }
  else if (m_drag_type == DragType::Gidaspow) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in,
                              ComputeDragGidaspow(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else if (m_drag_type == DragType::BVK2) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in,
                              ComputeDragBVK2(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else if (m_drag_type == DragType::UserDrag) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in,
                              ComputeDragUser(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else {
    amrex::Abort("Invalid Drag Type.");
  }
}

template <typename F1>
void mfix::mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      Vector< MultiFab* > const& X_gk_in,
                                      F1 DragFunc)
{
  if (advect_enthalpy)
  {
    if (m_convection_type == ConvectionType::RanzMarshall) {
        mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, DragFunc,
                                  ComputeConvRanzMarshall(DEM::small_number,DEM::large_number,DEM::eps));
    }
    else if (m_convection_type == ConvectionType::Gunn) {
      mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, DragFunc,
                                ComputeConvGunn(DEM::small_number,DEM::large_number,DEM::eps));
    }
    else {
      amrex::Abort("Invalid Convection Type.");
    }
  }
  else {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, X_gk_in, DragFunc,
                              NullConvectionCoeff());
  }
}

template <typename F1, typename F2>
void mfix::mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      Vector< MultiFab* > const& X_gk_in,
                                      F1 DragFunc,
                                      F2 ConvectionCoeff)
{
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  BL_PROFILE("mfix::mfix_calc_transfer_coeff()");

  const int run_on_device = Gpu::inLaunchRegion() ? 1 : 0;

  // We copy the value inside the domain to the outside to avoid
  // unphysical volume fractions.
  const int dir_bc_in = 2;
  mfix_set_epg_bcs(ep_g_in, dir_bc_in);

  // This is just a sanity check to make sure we're not using covered values
  // We can remove these lines once we're confident in the algorithm
  EB_set_covered(*vel_g_in[0], 0, 3, 1, covered_val);
  EB_set_covered(*ep_g_in[0], 0, 1, 1, covered_val);
  EB_set_covered(*ro_g_in[0], 0, 1, 1, covered_val);

  if (advect_enthalpy) {
    EB_set_covered(*T_g_in[0], 0, 1, 1, covered_val);
  }

  if (advect_fluid_species) {
    EB_set_covered(*X_gk_in[0], 0, fluid.nspecies, 1, covered_val);
  }

  for (int lev = 0; lev < nlev; lev++) {
    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    MultiFab* interp_ptr;

    const int interp_ng = 1;    // Only one layer needed for interpolation
    const int interp_comp = 5 + // (3 vel_g + 1 ep_g + ro_g + T_g + X_gk)
                            1*int(advect_enthalpy) +
                            fluid.nspecies*int(advect_fluid_species);

    if (OnSameGrids)
    {
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

      if (advect_enthalpy) {
        // Copy fluid temperature
        MultiFab::Copy(*interp_ptr, *T_g_in[lev],  0, components_count, 1, interp_ng);
        components_count += 1;
      }

      if (advect_fluid_species) {
        // Copy volume fraction
        MultiFab::Copy(*interp_ptr, *X_gk_in[lev],  0, components_count, fluid.nspecies, interp_ng);
        components_count += fluid.nspecies;
      }

      interp_ptr->FillBoundary(geom[lev].periodicity());
    }
    else
    {
      const BoxArray&            pba = pc->ParticleBoxArray(lev);
      const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

      EBFArrayBoxFactory ebfactory_loc(*eb_levels[lev], geom[lev], pba, pdm,
                                      {nghost_eb_basic(), nghost_eb_volume(), nghost_eb_full()}, 
                                       EBSupport::full);

      // Temporary arrays  -- copies with no ghost cells
      //const int ng_to_copy = 0; UNUSED VARIABLE

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_comp, interp_ng, MFInfo(), ebfactory_loc);

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

      if (advect_enthalpy) {
        // Copy fluid temperature
        interp_ptr->ParallelCopy(*T_g_in[lev],  0, components_count, 1, interp_ng, interp_ng);
        components_count += 1;
      }

      if (advect_fluid_species) {
        // Copy fluid species
        interp_ptr->ParallelCopy(*X_gk_in[lev],  0, components_count, fluid.nspecies, interp_ng, interp_ng);
        components_count += fluid.nspecies;
      }

      interp_ptr->FillBoundary(geom[lev].periodicity());
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      const auto dxi_array = geom[lev].InvCellSizeArray();
      const auto dx_array  = geom[lev].CellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
      const RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
      const RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(interp_ptr->Factory());

      const auto cellcent = &(factory.getCentroid());
      const auto bndrycent = &(factory.getBndryCent());
      const auto areafrac = factory.getAreaFrac();

      auto& fluid_parms = *fluid.parameters;

      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();

        const int np = particles.size();

        Box bx = pti.tilebox();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  interp_fab = static_cast<EBFArrayBox const&>((*interp_ptr)[pti]);
        const EBCellFlagFab&  flags = interp_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered)
        {
          const auto& interp_array = interp_ptr->const_array(pti);

          const auto& flags_array  = flags.array();

          auto particles_ptr = particles().dataPtr();

          const int adv_enthalpy = advect_enthalpy;
          const int fluid_is_a_mixture = fluid.is_a_mixture;

          const Real mu_g0 = fluid.mu_g0;

          const int nspecies_g = fluid.nspecies;

          auto local_cg_dem = DEM::cg_dem;

          if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
          {
            amrex::ParallelFor(np,
              [particles_ptr,p_realarray,interp_array,DragFunc,ConvectionCoeff,
               plo,dxi,adv_enthalpy,fluid_is_a_mixture,mu_g0,fluid_parms,nspecies_g,
               interp_comp,local_cg_dem,run_on_device]
              AMREX_GPU_DEVICE (int ip) noexcept
            {
              MFIXParticleContainer::ParticleType& particle = particles_ptr[ip];

              GpuArray<Real,6+SPECIES::NMAX> interp_loc; // vel_g, ep_g, ro_g, T_g, X_gk
              interp_loc.fill(0.);

              trilinear_interp(particle.pos(), interp_loc.data(), interp_array,
                               plo, dxi, interp_comp);

              RealVect vel_g(interp_loc[0], interp_loc[1], interp_loc[2]);
              const Real ep_g = interp_loc[3];
              const Real ro_g = interp_loc[4];

              int comp_count = 5;

              Real T_g(0);
              if (adv_enthalpy) {
                T_g = interp_loc[comp_count];
                comp_count += 1;
              }

              Real* X_gk;
              if (fluid_is_a_mixture) {
                X_gk = &interp_loc[comp_count];
              }

              // Indices of cell where particle is located
              int iloc = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
              int jloc = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
              int kloc = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

              Real  mu_g(0);

              if (adv_enthalpy)
                mu_g = fluid_parms.calc_mu_g(T_g);
              else
                mu_g = mu_g0;

              Real rad = p_realarray[SoArealData::radius][ip];
              Real vol = p_realarray[SoArealData::volume][ip];

              int p_id = particle.id();

              RealVect pvel(0.);
              pvel[0] = p_realarray[SoArealData::velx][ip];
              pvel[1] = p_realarray[SoArealData::vely][ip];
              pvel[2] = p_realarray[SoArealData::velz][ip];

              Real rop_g = ro_g * ep_g;

              RealVect vslp(0.);
              vslp[0] = vel_g[0] - pvel[0];
              vslp[1] = vel_g[1] - pvel[1];
              vslp[2] = vel_g[2] - pvel[2];

              Real vrel = sqrt(dot_product(vslp, vslp));
              Real dp = 2.0*rad;

              if (local_cg_dem) {
                 dp = dp/std::cbrt(p_realarray[SoArealData::statwt][ip]);
              }

              Real ep_s = 1.0 - ep_g;

              Real beta = vol*DragFunc(ep_g, mu_g, rop_g, vrel, dp, dp, ep_s,
                 vel_g[0], vel_g[1], vel_g[2], iloc, jloc, kloc, p_id);

              p_realarray[SoArealData::dragcoeff][ip] = beta;

              if(adv_enthalpy){
                Real k_g = fluid_parms.calc_k_g(T_g);

                Real cp_g(0);
                if (!fluid_is_a_mixture)
                  cp_g = run_on_device ?
                    fluid_parms.calc_cp_g<RunOn::Device>(T_g) :
                    fluid_parms.calc_cp_g<RunOn::Host>(T_g);
                else {
                  for (int n_g(0); n_g < nspecies_g; ++n_g) {
                    cp_g += run_on_device ?
                      X_gk[n_g]*fluid_parms.calc_cp_gk<RunOn::Device>(T_g,n_g) :
                      X_gk[n_g]*fluid_parms.calc_cp_gk<RunOn::Host>(T_g,n_g);
                  }
                }

                Real gamma = ConvectionCoeff(ep_g, mu_g, k_g, cp_g, rop_g, vrel,
                                             dp, iloc, jloc, kloc, p_id);

                p_realarray[SoArealData::convection][ip] = 4.0*M_PI*rad*rad*gamma;
              }

            });
          }
          else // FAB not all regular
          {
            // Cell centroids
            const auto& ccent_fab = cellcent->array(pti);
            // Centroid of EB
            const auto& bcent_fab = bndrycent->array(pti);
            // Area fractions
            const auto& apx_fab = areafrac[0]->array(pti);
            const auto& apy_fab = areafrac[1]->array(pti);
            const auto& apz_fab = areafrac[2]->array(pti);

            //const int adv_enthalpy = advect_enthalpy;
            //const Real mu_g0 = fluid.mu_g0;

            amrex::ParallelFor(np,
              [particles_ptr,p_realarray,interp_array,DragFunc,ConvectionCoeff,
               adv_enthalpy,fluid_is_a_mixture,mu_g0,fluid_parms,plo,dx,dxi,flags_array,
               ccent_fab,bcent_fab,apx_fab,apy_fab,apz_fab,nspecies_g,interp_comp,
               local_cg_dem,run_on_device]
              AMREX_GPU_DEVICE (int pid) noexcept
            {
              MFIXParticleContainer::ParticleType& particle = particles_ptr[pid];

              // Cell containing particle centroid
              int ip = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
              int jp = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
              int kp = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

              // No drag force for particles in covered cells.
              if (flags_array(ip,jp,kp).isCovered()) {

                p_realarray[SoArealData::dragcoeff][pid] = 0.;

              // Cut or regular cell and none of the cells in the stencil is
              // covered (Note we can't assume regular cell has no covered
              // cells in the stencil because of the diagonal case)
              } else {

                // Upper cell in trilinear stencil
                int i = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0] + 0.5));
                int j = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1] + 0.5));
                int k = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2] + 0.5));

                // Local array storing interpolated values
                GpuArray< Real, 6+SPECIES::NMAX> interp_loc; // vel_g, ep_g, ro_g, T_g, X_gk
                interp_loc.fill(0.);

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

                  trilinear_interp(particle.pos(), interp_loc.data(),
                                   interp_array, plo, dxi, interp_comp);
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

                    shepard_interp_eb(particle.pos(), ip, jp, kp, dx, dxi, plo,
                                      flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                      interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);

                  }
                  {
                    const int srccomp = 3;
                    const int dstcomp = 3;
                    const int numcomp = interp_comp-3; // ep_g, ro_g, T_g, X_gk

                    shepard_interp(particle.pos(), ip, jp, kp, dx, dxi, plo,
                                   flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                   interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);
                  }
#else
                  const int srccomp = 0;
                  const int dstcomp = 0;
                  const int numcomp = interp_comp; // vel_g, ep_g, ro_g, T_g, X_gk

                  shepard_interp(particle.pos(), ip, jp, kp, dx, dxi, plo,
                                 flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                                 interp_array, interp_loc.data(), srccomp, dstcomp, numcomp);

#endif
                } // Cut cell

                RealVect vel_g(interp_loc[0], interp_loc[1], interp_loc[2]);
                const Real ep_g = interp_loc[3];
                const Real ro_g = interp_loc[4];

                int comp_count(5);

                Real T_g(0);
                if (adv_enthalpy) {
                  T_g = interp_loc[comp_count];
                  comp_count += 1;
                }

                Real* X_gk;
                if (fluid_is_a_mixture) {
                  X_gk = &interp_loc[comp_count];
                }

                Real mu_g(0);

                if (adv_enthalpy)
                  mu_g = fluid_parms.calc_mu_g(T_g);
                else
                  mu_g = mu_g0;

                Real rad = p_realarray[SoArealData::radius][pid];
                Real vol = p_realarray[SoArealData::volume][pid];

                int p_id = particle.id();

                RealVect pvel(0.);
                pvel[0] = p_realarray[SoArealData::velx][pid];
                pvel[1] = p_realarray[SoArealData::vely][pid];
                pvel[2] = p_realarray[SoArealData::velz][pid];

                Real rop_g = ro_g * ep_g;
                RealVect vslp(0.);
                vslp[0] = vel_g[0] - pvel[0];
                vslp[1] = vel_g[1] - pvel[1];
                vslp[2] = vel_g[2] - pvel[2];

                Real vrel = sqrt(dot_product(vslp, vslp));
                Real dp = 2.0*rad;

                if (local_cg_dem) {
                   dp = dp/std::cbrt(p_realarray[SoArealData::statwt][pid]);
                }

                Real ep_s = 1.0 - ep_g;
                Real beta = vol*DragFunc(ep_g, mu_g, rop_g, vrel, dp, dp, ep_s,
                                         vel_g[0], vel_g[1], vel_g[2],
                                         ip, jp, kp, p_id);

                p_realarray[SoArealData::dragcoeff][pid] = beta;

                if(adv_enthalpy) {
                  Real k_g = fluid_parms.calc_k_g(T_g);

                  Real cp_g(0.);
                  if (!fluid_is_a_mixture)
                    cp_g = run_on_device ?
                      fluid_parms.calc_cp_g<RunOn::Device>(T_g) :
                      fluid_parms.calc_cp_g<RunOn::Host>(T_g);
                  else {
                    for (int n_g(0); n_g < nspecies_g; ++n_g) {
                      cp_g += run_on_device ?
                        X_gk[n_g]*fluid_parms.calc_cp_gk<RunOn::Device>(T_g,n_g) :
                        X_gk[n_g]*fluid_parms.calc_cp_gk<RunOn::Host>(T_g,n_g);
                    }
                  }

                  Real gamma = ConvectionCoeff(ep_g, mu_g, k_g, cp_g, rop_g, vrel,
                                               dp, ip, jp, kp, p_id);

                  p_realarray[SoArealData::convection][pid] = 4.0*M_PI*rad*rad*gamma;
                }

              } // Not covered
            }); // pid
          } // type of FAB
        } // if entire FAB not covered
      } // pti
    } // GPU region

    delete interp_ptr;
  } // lev


  // Reset the volume fractions back to the correct values at
  // inflow faces.
  const int dir_bc_out = 1;
  mfix_set_epg_bcs(ep_g_in, dir_bc_out);

}
