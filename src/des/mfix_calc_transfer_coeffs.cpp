#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_K.H>
#include <mfix_des_drag_K.H>
#include <mfix_des_conv_coeff_K.H>
#include <mfix_mf_helpers.H>

void mfix::mfix_calc_transfer_coeffs ()
{
  if (m_drag_type == DragType::WenYu) {
    mfix_calc_transfer_coeffs(ComputeDragWenYu(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else if (m_drag_type == DragType::Gidaspow) {
    mfix_calc_transfer_coeffs(ComputeDragGidaspow(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else if (m_drag_type == DragType::BVK2) {
    mfix_calc_transfer_coeffs(ComputeDragBVK2(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else if (m_drag_type == DragType::UserDrag) {
    mfix_calc_transfer_coeffs(ComputeDragUser(DEM::small_number,DEM::large_number,DEM::eps));
  }
  else {
    amrex::Abort("Invalid Drag Type.");
  }
}

template <typename F1>
void mfix::mfix_calc_transfer_coeffs (F1 DragFunc)
{
  if (advect_enthalpy)
  {
    if (m_convection_type == ConvectionType::RanzMarshall) {
        mfix_calc_transfer_coeffs(DragFunc, ComputeConvRanzMarshall(DEM::small_number,DEM::large_number,DEM::eps));
    }
    else if (m_convection_type == ConvectionType::Gunn) {
      mfix_calc_transfer_coeffs(DragFunc, ComputeConvGunn(DEM::small_number,DEM::large_number,DEM::eps));
    }
    else {
      amrex::Abort("Invalid Convection Type.");
    }
  }
  else {
    mfix_calc_transfer_coeffs(DragFunc, NullConvectionCoeff());
  }
}

template <typename F1, typename F2>
void mfix::mfix_calc_transfer_coeffs (F1 DragFunc, F2 ConvectionCoeff)
{
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  BL_PROFILE("mfix::mfix_calc_transfer_coeff()");

  // We copy the value inside the domain to the outside to avoid
  // unphysical volume fractions.
  const int dir_bc_in = 2;
  mfix_set_epg_bcs(get_ep_g(), dir_bc_in);

  // This is just a sanity check to make sure we're not using covered values
  // We can remove these lines once we're confident in the algorithm
  EB_set_covered(*m_leveldata[0]->vel_g, 0, 3, 1, covered_val);
  EB_set_covered(*m_leveldata[0]->ep_g,  0, 1, 1, covered_val);
  EB_set_covered(*m_leveldata[0]->mu_g,  0, 1, 1, covered_val);
  EB_set_covered(*m_leveldata[0]->ro_g,  0, 1, 1, covered_val);

  if (advect_enthalpy) {
    EB_set_covered(*m_leveldata[0]->cp_g,  0, 1, 1, covered_val);
    EB_set_covered(*m_leveldata[0]->k_g,   0, 1, 1, covered_val);
  }

  for (int lev = 0; lev < nlev; lev++)
  {
    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) and
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    MultiFab* ro_ptr;
    MultiFab* mu_ptr;
    MultiFab* cp_ptr;
    MultiFab* kg_ptr;
    MultiFab* interp_ptr;

    const int interp_ng = 1;    // Only one layer needed for interpolation
    const int interp_comp = 4;  // Four components (3 vel_g + 1 ep_g)

    if (OnSameGrids)
    {
      ro_ptr = m_leveldata[lev]->ro_g;
      mu_ptr = m_leveldata[lev]->mu_g;

      if (advect_enthalpy) {
        kg_ptr = m_leveldata[lev]->k_g;
        cp_ptr = m_leveldata[lev]->cp_g;
      }
      else {
        kg_ptr = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
        cp_ptr = new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *ebfactory[lev]);
      }

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(grids[lev], dmap[lev], interp_comp, interp_ng, MFInfo(), *ebfactory[lev]);

      // Copy fluid velocity
      interp_ptr->copy(*m_leveldata[lev]->vel_g, 0, 0,
                        m_leveldata[lev]->vel_g->nComp(),
                        interp_ng, interp_ng);

      // Copy volume fraction
      interp_ptr->copy(*m_leveldata[lev]->ep_g,  0, 3,
                        m_leveldata[lev]->ep_g->nComp(),
                        interp_ng, interp_ng);

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

      ro_ptr = new MultiFab(pba, pdm, m_leveldata[lev]->ro_g->nComp(), 1);
      ro_ptr->copy(*m_leveldata[lev]->ro_g, 0, 0, 1, ng_to_copy, ng_to_copy);

      mu_ptr = new MultiFab(pba, pdm, m_leveldata[lev]->mu_g->nComp(), 1);
      mu_ptr->copy(*m_leveldata[lev]->mu_g, 0, 0, 1, ng_to_copy, ng_to_copy);

      if (advect_enthalpy) {
        kg_ptr = new MultiFab(pba, pdm, m_leveldata[lev]->k_g->nComp(), 1);
        kg_ptr->copy(*m_leveldata[lev]->k_g, 0, 0, 1, ng_to_copy, ng_to_copy);

        cp_ptr = new MultiFab(pba, pdm, m_leveldata[lev]->cp_g->nComp(), 1);
        cp_ptr->copy(*m_leveldata[lev]->cp_g, 0, 0, 1, ng_to_copy, ng_to_copy);
      }
      else {
        kg_ptr = new MultiFab(pba, pdm, 1, 0, MFInfo(), ebfactory_loc);
        cp_ptr = new MultiFab(pba, pdm, 1, 0, MFInfo(), ebfactory_loc);
      }

      // Store gas velocity and volume fraction for interpolation
      interp_ptr = new MultiFab(pba, pdm, interp_comp, interp_ng, MFInfo(), ebfactory_loc);

      // Copy fluid velocity
      interp_ptr->copy(*m_leveldata[lev]->vel_g, 0, 0,
                        m_leveldata[lev]->vel_g->nComp(),
                        interp_ng, interp_ng);

      // Copy volume fraction
      interp_ptr->copy(*m_leveldata[lev]->ep_g,  0, 3,
                        m_leveldata[lev]->ep_g->nComp(),
                        interp_ng, interp_ng);

      interp_ptr->FillBoundary(geom[lev].periodicity());
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      const auto dxi_array = geom[lev].InvCellSizeArray();
      const auto dx_array  = geom[lev].CellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const amrex::RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
      const amrex::RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
      const amrex::RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(interp_ptr->Factory());

      const auto cellcent = &(factory.getCentroid());
      const auto bndrycent = &(factory.getBndryCent());
      const auto areafrac = factory.getAreaFrac();

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
          const auto& mu_array     = mu_ptr->array(pti);
          const auto& kg_array     = kg_ptr->array(pti);
          const auto& cp_array     = cp_ptr->array(pti);

          const auto& flags_array  = flags.array();

          auto particles_ptr = particles().dataPtr();

          if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
          {
            amrex::ParallelFor(np,
              [particles_ptr,p_realarray,interp_array,ro_array,mu_array,kg_array,cp_array,
               DragFunc, ConvectionCoeff,plo,dxi,
               local_cg_dem=DEM::cg_dem, local_advect_enthalpy=advect_enthalpy]
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
              Real  mu = mu_array(iloc,jloc,kloc);

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

              if(local_advect_enthalpy){
                Real kg = kg_array(iloc,jloc,kloc);
                Real cp = cp_array(iloc,jloc,kloc);
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

            amrex::ParallelFor(np,
              [particles_ptr,p_realarray,interp_array,ro_array,mu_array,kg_array,cp_array,
               DragFunc, ConvectionCoeff,
               plo,dx,dxi,flags_array,ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
               local_cg_dem=DEM::cg_dem,local_advect_enthalpy=advect_enthalpy]
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
                if (flags_array(i-1,j-1,k-1).isRegular() and
                    flags_array(i  ,j-1,k-1).isRegular() and
                    flags_array(i-1,j  ,k-1).isRegular() and
                    flags_array(i  ,j  ,k-1).isRegular() and
                    flags_array(i-1,j-1,k  ).isRegular() and
                    flags_array(i  ,j-1,k  ).isRegular() and
                    flags_array(i-1,j  ,k  ).isRegular() and
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
                Real  mu = mu_array(ip,jp,kp);

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

              if(local_advect_enthalpy){
                Real kg = kg_array(ip,jp,kp);
                Real cp = cp_array(ip,jp,kp);
                Real gamma = ConvectionCoeff(ep, mu, kg, cp, rop_g, vrel, dp, ip, jp, kp, p_id);
                p_realarray[SoArealData::convection][pid] = 4.0*M_PI*rad*rad*gamma;
              }

              } // Not covered
            }); // pid
          } // type of FAB
        } // if entire FAB not covered
      } // pti
    } // GPU region

    if (not OnSameGrids) {
      delete ro_ptr;
      delete mu_ptr;
      delete kg_ptr;
      delete cp_ptr;
    }
    else if (not advect_enthalpy) {
      delete kg_ptr;
      delete cp_ptr;
    }

    delete interp_ptr;
  } // lev


  // Reset the volume fractions back to the correct values at
  // inflow faces.
  const int dir_bc_out = 1;
  mfix_set_epg_bcs(get_ep_g(), dir_bc_out);

}
