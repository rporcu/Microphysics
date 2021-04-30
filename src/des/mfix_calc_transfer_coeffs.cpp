#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_K.H>
#include <mfix_des_drag_K.H>
#include <mfix_des_conv_coeff_K.H>
#include <mfix_mf_helpers.H>

void mfix::mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in)
{
  if (m_drag_type == DragType::WenYu) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in,
                              ComputeDragWenYu(DEM::small_number, DEM::large_number, DEM::eps));
  }
  else if (m_drag_type == DragType::Gidaspow) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in,
                              ComputeDragGidaspow(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else if (m_drag_type == DragType::BVK2) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in,
                              ComputeDragBVK2(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else if (m_drag_type == DragType::UserDrag) {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in,
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
                                      F1 DragFunc)
{
  if (advect_enthalpy)
  {
    if (m_convection_type == ConvectionType::RanzMarshall) {
        mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, DragFunc,
                                  ComputeConvRanzMarshall(DEM::small_number,DEM::large_number,DEM::eps));
    }
    else if (m_convection_type == ConvectionType::Gunn) {
      mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, DragFunc,
                                ComputeConvGunn(DEM::small_number,DEM::large_number,DEM::eps));
    }
    else {
      amrex::Abort("Invalid Convection Type.");
    }
  }
  else {
    mfix_calc_transfer_coeffs(ep_g_in, ro_g_in, vel_g_in, T_g_in, DragFunc,
                              NullConvectionCoeff());
  }
}

template <typename F1, typename F2>
void mfix::mfix_calc_transfer_coeffs (Vector< MultiFab* > const& ep_g_in,
                                      Vector< MultiFab* > const& ro_g_in,
                                      Vector< MultiFab* > const& vel_g_in,
                                      Vector< MultiFab* > const& T_g_in,
                                      F1 DragFunc,
                                      F2 ConvectionCoeff)
{
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  BL_PROFILE("mfix::mfix_calc_transfer_coeff()");

  // We copy the value inside the domain to the outside to avoid
  // unphysical volume fractions.
  const int dir_bc_in = 2;
  mfix_set_epg_bcs(ep_g_in, dir_bc_in);

  // This is just a sanity check to make sure we're not using covered values
  // We can remove these lines once we're confident in the algorithm
  EB_set_covered(*vel_g_in[0], 0, 3, 1, covered_val);
  EB_set_covered(*ep_g_in[0],  0, 1, 1, covered_val);
  EB_set_covered(*ro_g_in[0],  0, 1, 1, covered_val);

  if (advect_enthalpy) {
    EB_set_covered(*T_g_in[0],  0, 1, 1, covered_val);
  }

  for (int lev = 0; lev < nlev; lev++)
  {
    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) &&
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    MultiFab* ro_ptr;
    MultiFab* T_ptr;
    MultiFab* interp_ptr;

    const int interp_ng = 1;    // Only one layer needed for interpolation
    const int interp_comp = 4;  // Four components (3 vel_g + 1 ep_g)

    if (OnSameGrids)
    {
      ro_ptr = ro_g_in[lev];

      if (advect_enthalpy) {
        T_ptr = T_g_in[lev];
      }
      else {
        T_ptr = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
      }

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(grids[lev], dmap[lev], interp_comp, interp_ng, MFInfo(), *ebfactory[lev]);

      // Copy fluid velocity
      MultiFab::Copy(*interp_ptr, *vel_g_in[lev], 0, 0, vel_g_in[lev]->nComp(), interp_ng);

      // Copy volume fraction
      MultiFab::Copy(*interp_ptr, *ep_g_in[lev],  0, 3, ep_g_in[lev]->nComp(), interp_ng);

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
      const int ng_to_copy = 0;

      ro_ptr = new MultiFab(pba, pdm, ro_g_in[lev]->nComp(), 1);
      ro_ptr->copy(*ro_g_in[lev], 0, 0, 1, ng_to_copy, ng_to_copy);

      if (advect_enthalpy) {
        T_ptr = new MultiFab(pba, pdm, T_g_in[lev]->nComp(), 1);
        T_ptr->copy(*T_g_in[lev], 0, 0, 1, ng_to_copy, ng_to_copy);
      }
      else {
        T_ptr = new MultiFab(pba, pdm, 1, 0, MFInfo(), ebfactory_loc);
      }

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_comp, interp_ng, MFInfo(), ebfactory_loc);

      // Copy fluid velocity
      interp_ptr->copy(*vel_g_in[lev], 0, 0, vel_g_in[lev]->nComp(), interp_ng, interp_ng);

      // Copy volume fraction
      interp_ptr->copy(*ep_g_in[lev],  0, 3, ep_g_in[lev]->nComp(), interp_ng, interp_ng);

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
          const auto& interp_array = interp_ptr->array(pti);
          const auto& ro_array     = ro_ptr->array(pti);
          const auto& T_array      = advect_enthalpy ? T_ptr->array(pti) : Array4<const Real>();

          const auto& flags_array  = flags.array();

          auto particles_ptr = particles().dataPtr();

          const int adv_enthalpy = advect_enthalpy;
          const Real mu_g0 = fluid.mu_g0;

          if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
          {
            amrex::ParallelFor(np,
              [particles_ptr,p_realarray,interp_array,ro_array,T_array,
               DragFunc,ConvectionCoeff,plo,dxi,adv_enthalpy,mu_g0,fluid_parms,
               local_cg_dem=DEM::cg_dem]
              AMREX_GPU_DEVICE (int ip) noexcept
            {
              MFIXParticleContainer::ParticleType& particle = particles_ptr[ip];

              GpuArray< Real, interp_comp> interp_loc;
              trilinear_interp(particle.pos(), interp_loc.data(),
                               interp_array, plo, dxi, interp_comp);

              RealVect velfp(0.);
              Real ep(0.);

              velfp[0] = interp_loc[0];
              velfp[1] = interp_loc[1];
              velfp[2] = interp_loc[2];
              ep       = interp_loc[3];

              // Indices of cell where particle is located
              int iloc = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
              int jloc = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
              int kloc = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

              Real  ro = ro_array(iloc,jloc,kloc);
              Real  mu(0);

              if (adv_enthalpy)
                mu = fluid_parms.calc_mu_g(T_array(iloc,jloc,kloc));
              else
                mu = mu_g0;

              Real rad = p_realarray[SoArealData::radius][ip];
              Real vol = p_realarray[SoArealData::volume][ip];

              int p_id = particle.id();

              RealVect pvel(0.);
              pvel[0] = p_realarray[SoArealData::velx][ip];
              pvel[1] = p_realarray[SoArealData::vely][ip];
              pvel[2] = p_realarray[SoArealData::velz][ip];

              Real rop_g = ro * ep;

              RealVect vslp(0.);
              vslp[0] = velfp[0] - pvel[0];
              vslp[1] = velfp[1] - pvel[1];
              vslp[2] = velfp[2] - pvel[2];

              Real vrel = sqrt(dot_product(vslp, vslp));
              Real dp = 2.0*rad;
              if (local_cg_dem)
              {
                 dp = dp/std::cbrt(p_realarray[SoArealData::statwt][ip]);
              }
              Real phis = 1.0 - ep;
              Real beta = vol*DragFunc(ep, mu, rop_g, vrel, dp, dp, phis,
                 velfp[0], velfp[1], velfp[2], iloc, jloc, kloc, p_id);

              p_realarray[SoArealData::dragcoeff][ip] = beta;

              if(adv_enthalpy){
                Real kg = fluid_parms.calc_k_g(T_array(iloc,jloc,kloc));
                Real cp = fluid_parms.calc_cp_g(T_array(iloc,jloc,kloc));
                Real gamma = ConvectionCoeff(ep, mu, kg, cp, rop_g, vrel, dp, iloc, jloc, kloc, p_id);
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

            const int adv_enthalpy = advect_enthalpy;
            const Real mu_g0 = fluid.mu_g0;

            amrex::ParallelFor(np,
              [particles_ptr,p_realarray,interp_array,ro_array,T_array,
               DragFunc,ConvectionCoeff,adv_enthalpy,mu_g0,fluid_parms,
               plo,dx,dxi,flags_array,ccent_fab, bcent_fab,apx_fab,apy_fab,apz_fab,
               local_cg_dem=DEM::cg_dem]
              AMREX_GPU_DEVICE (int pid) noexcept
            {
              MFIXParticleContainer::ParticleType& particle = particles_ptr[pid];

              // Cell containing particle centroid
              int ip = static_cast<int>(amrex::Math::floor((particle.pos(0) - plo[0])*dxi[0]));
              int jp = static_cast<int>(amrex::Math::floor((particle.pos(1) - plo[1])*dxi[1]));
              int kp = static_cast<int>(amrex::Math::floor((particle.pos(2) - plo[2])*dxi[2]));

              // No drag force for particles in covered cells.
              if (flags_array(ip,jp,kp).isCovered() ){

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
                GpuArray< Real, interp_comp> interp_loc;

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

                  const int scomp = 3;
                  fe_interp(particle.pos(), ip, jp, kp, dx, dxi, plo,
                            flags_array, ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                            interp_array, interp_loc.data(), interp_comp, scomp);
                } // Cut cell

                RealVect velfp(0.);
                Real ep(0.);

                velfp[0] = interp_loc[0];
                velfp[1] = interp_loc[1];
                velfp[2] = interp_loc[2];
                ep       = interp_loc[3];

                // Using i/j/k of centroid cell
                Real  ro = ro_array(ip,jp,kp);

                Real mu(0);

                if (adv_enthalpy)
                  mu = fluid_parms.calc_mu_g(T_array(ip,jp,kp));
                else
                  mu = mu_g0;

                Real rad = p_realarray[SoArealData::radius][pid];
                Real vol = p_realarray[SoArealData::volume][pid];

                int p_id = particle.id();

                RealVect pvel(0.);
                pvel[0] = p_realarray[SoArealData::velx][pid];
                pvel[1] = p_realarray[SoArealData::vely][pid];
                pvel[2] = p_realarray[SoArealData::velz][pid];

                Real rop_g = ro * ep;
                RealVect vslp(0.);
                vslp[0] = velfp[0] - pvel[0];
                vslp[1] = velfp[1] - pvel[1];
                vslp[2] = velfp[2] - pvel[2];

                Real vrel = sqrt(dot_product(vslp, vslp));
                Real dp = 2.0*rad;
                if (local_cg_dem)
                {
                   dp = dp/std::cbrt(p_realarray[SoArealData::statwt][pid]);
                }
                Real phis = 1.0 - ep;
                Real beta = vol*DragFunc(ep, mu, rop_g, vrel, dp, dp, phis,
                                         velfp[0], velfp[1], velfp[2],
                                         ip, jp, kp, p_id);

                p_realarray[SoArealData::dragcoeff][pid] = beta;

              if(adv_enthalpy) {
                Real kg = fluid_parms.calc_k_g(T_array(ip,jp,kp));
                Real cp = fluid_parms.calc_cp_g(T_array(ip,jp,kp));
                Real gamma = ConvectionCoeff(ep, mu, kg, cp, rop_g, vrel, dp, ip, jp, kp, p_id);
                p_realarray[SoArealData::convection][pid] = 4.0*M_PI*rad*rad*gamma;
              }

              } // Not covered
            }); // pid
          } // type of FAB
        } // if entire FAB not covered
      } // pti
    } // GPU region

    if (!OnSameGrids) {
      delete ro_ptr;

      if (!advect_enthalpy)
      delete T_ptr;
    }
    else if (!advect_enthalpy) {
      delete T_ptr;
    }

    delete interp_ptr;
  } // lev


  // Reset the volume fractions back to the correct values at
  // inflow faces.
  const int dir_bc_out = 1;
  mfix_set_epg_bcs(ep_g_in, dir_bc_out);

}
