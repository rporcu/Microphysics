#include <mfix.H>
#include <mfix_des_K.H>
#include <mfix_drag_K.H>
#include <des_drag_K.H>

void mfix::mfix_calc_particle_beta (Real time)
{
  if (m_drag_type == DragType::WenYu) {
    mfix_calc_particle_beta(ComputeDragWenYu(), time);
  }
  else if (m_drag_type == DragType::Gidaspow) {
    mfix_calc_particle_beta(ComputeDragGidaspow(), time);
  }
  else if (m_drag_type == DragType::BVK2) {
    mfix_calc_particle_beta(ComputeDragBVK2(), time);
  }
  else if (m_drag_type == DragType::UserDrag) {
    mfix_calc_particle_beta(ComputeDragUser(), time);
  }
    else {
    amrex::Abort("Invalid Drag Type.");
  }
}

template <typename F>
void mfix::mfix_calc_particle_beta (F DragFunc, Real time)
{
  using MFIXParIter = MFIXParticleContainer::MFIXParIter;

  BL_PROFILE("mfix::mfix_calc_particle_beta()");

  // This is just a sanity check to make sure we're not using covered values
  // We can remove these lines once we're confident in the algoirthm 
  EB_set_covered(*m_leveldata[0]->vel_g, 0, 3, 1, covered_val);
  EB_set_covered(*m_leveldata[0]->ep_g, 0, 1, 1, covered_val);
  EB_set_covered(*m_leveldata[0]->mu_g, 0, 1, 1, covered_val);
  EB_set_covered(*m_leveldata[0]->ro_g, 0, 1, 1, covered_val);

  for (int lev = 0; lev < nlev; lev++)
  {
    bool OnSameGrids = ( (dmap[lev] == (pc->ParticleDistributionMap(lev))) and
                         (grids[lev].CellEqual(pc->ParticleBoxArray(lev))) );

    MultiFab* ep_ptr;
    MultiFab* ro_ptr;
    MultiFab* mu_ptr;
    MultiFab* vel_ptr;

    if (OnSameGrids)
    {
      ep_ptr  = m_leveldata[lev]->ep_g;
      ro_ptr  = m_leveldata[lev]->ro_g;
      mu_ptr  = m_leveldata[lev]->mu_g;
      vel_ptr = m_leveldata[lev]->vel_g;
    }
    else
    {
      const BoxArray&            pba = pc->ParticleBoxArray(lev);
      const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

      // Temporary arrays  -- copies with no ghost cells 
      ep_ptr = new MultiFab(pba, pdm, m_leveldata[lev]->ep_g->nComp(), 0);
      ep_ptr->copy(*m_leveldata[lev]->ep_g, 0, 0, 1, 0, 0);

      ro_ptr = new MultiFab(pba, pdm, m_leveldata[lev]->ro_g->nComp(), 0);
      ro_ptr->copy(*m_leveldata[lev]->ro_g, 0, 0, 1, 0, 0);

      mu_ptr = new MultiFab(pba, pdm, m_leveldata[lev]->mu_g->nComp(), 0);
      mu_ptr->copy(*m_leveldata[lev]->mu_g, 0, 0, 1, 0, 0);

      EBFArrayBoxFactory ebfactory_loc(*eb_levels[lev], geom[lev], pba, pdm,
                                       {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                                       EBSupport::basic);

      int ng = m_leveldata[lev]->vel_g->nGrow();
      vel_ptr = new MultiFab(pba, pdm, m_leveldata[lev]->vel_g->nComp(), ng,
                             MFInfo(), ebfactory_loc);
      vel_ptr->copy(*m_leveldata[lev]->vel_g, 0, 0,
                    m_leveldata[lev]->vel_g->nComp(), ng, ng);
      vel_ptr->FillBoundary(geom[lev].periodicity());
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      const auto dxi_array = geom[lev].InvCellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const amrex::RealVect dxi(dxi_array[0], dxi_array[1], dxi_array[2]);
      const amrex::RealVect plo(plo_array[0], plo_array[1], plo_array[2]);

      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        const int np = particles.size();

        Box bx = pti.tilebox();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_ptr)[pti]);
        const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered)
        {
          const auto& vel_array   = vel_ptr->array(pti);
          const auto& ep_array    =  ep_ptr->array(pti);
          const auto& ro_array    =  ro_ptr->array(pti);
          const auto& mu_array    =  mu_ptr->array(pti);
          const auto& flags_array = flags.array();

          auto particles_ptr = particles().dataPtr();

          if (flags.getType(amrex::grow(bx,1)) == FabType::regular)
          {
            amrex::ParallelFor(np,
              [particles_ptr,vel_array,ep_array,ro_array,mu_array,DragFunc,plo,dxi]
              AMREX_GPU_DEVICE (int ip) noexcept
            {
              MFIXParticleContainer::ParticleType& particle = particles_ptr[ip];

              Real velfp[3];
              trilinear_interp(particle.pos(), &velfp[0], vel_array, plo, dxi);

              // Indices of cell where particle is located
              int iloc = floor((particle.pos(0) - plo[0])*dxi[0]);
              int jloc = floor((particle.pos(1) - plo[1])*dxi[1]);
              int kloc = floor((particle.pos(2) - plo[2])*dxi[2]);

              Real  ep = ep_array(iloc,jloc,kloc);
              Real  ro = ro_array(iloc,jloc,kloc);
              Real  mu = mu_array(iloc,jloc,kloc);

              Real rad = particle.rdata(realData::radius);
              Real vol = particle.rdata(realData::volume);

              int p_id = particle.id();

              Real pvel[3]; 
              pvel[0] = particle.rdata(realData::velx);
              pvel[1] = particle.rdata(realData::vely);
              pvel[2] = particle.rdata(realData::velz);

              Real rop_g = ro * ep;

              Real vslp[3];
              vslp[0] = velfp[0] - pvel[0];
              vslp[1] = velfp[1] - pvel[1];
              vslp[2] = velfp[2] - pvel[2];
              
              Real vrel = sqrt(dot_product(vslp, vslp));
              Real dpm = 2.0*rad;
              Real phis = 1.0 - ep;   
              Real beta = vol*DragFunc(ep, mu, rop_g, vrel, dpm, dpm, phis,
                 velfp[0], velfp[1], velfp[2], iloc, jloc, kloc, p_id); 

              particle.rdata(realData::dragx) = beta;
            });
          }
          else // FAB not all regular
          {
            // phi is always on the particle grids -- we do not use the refined
            // level here
            amrex::ParallelFor(np,
              [particles_ptr,vel_array,ro_array,mu_array,ep_array,flags_array,DragFunc,plo,dxi]
              AMREX_GPU_DEVICE (int ip) noexcept
            {
              MFIXParticleContainer::ParticleType& particle = particles_ptr[ip];

              Real velfp[3];

              // This identifies which cell the particle is in
              int iloc = floor((particle.pos(0) - plo[0])*dxi[0]);
              int jloc = floor((particle.pos(1) - plo[1])*dxi[1]); 
              int kloc = floor((particle.pos(2) - plo[2])*dxi[2]);

              // Pick upper cell in the stencil
              Real lx = (particle.pos(0) - plo[0])*dxi[0] + 0.5;
              Real ly = (particle.pos(1) - plo[1])*dxi[1] + 0.5;
              Real lz = (particle.pos(2) - plo[2])*dxi[2] + 0.5;

              int i = std::floor(lx);
              int j = std::floor(ly);
              int k = std::floor(lz);

              // Covered cell
              if (flags_array(iloc,jloc,kloc).isCovered())
              {
                particle.rdata(realData::dragx) = 0.0;
              }
              else
              {
                // Cut or regular cell and none of the cells in the stencil is
                // covered (Note we can't assume regular cell has no covered
                // cells in the stencil because of the diagonal case)
                if (!flags_array(i-1,j-1,k-1).isCovered() and
                    !flags_array(i  ,j-1,k-1).isCovered() and
                    !flags_array(i-1,j  ,k-1).isCovered() and
                    !flags_array(i  ,j  ,k-1).isCovered() and
                    !flags_array(i-1,j-1,k  ).isCovered() and
                    !flags_array(i  ,j-1,k  ).isCovered() and
                    !flags_array(i-1,j  ,k  ).isCovered() and
                    !flags_array(i  ,j  ,k  ).isCovered()) 
                {
                  trilinear_interp(particle.pos(), &velfp[0], vel_array, plo, dxi);
                // At least one of the cells in the stencil is covered
                }
                else
                {
                  // Particle position must be in [-.5:.5] is relative to cell
                  // center and scaled by dx
                  Real gx = particle.pos(0)*dxi[0] - (iloc + 0.5);
                  Real gy = particle.pos(1)*dxi[1] - (jloc + 0.5);
                  Real gz = particle.pos(2)*dxi[2] - (kloc + 0.5);

                  int ii;
                  int jj;
                  int kk;

                  if (not flags_array(iloc-1, jloc, kloc).isCovered())
                  {
                    ii = iloc - 1;
                  }
                  else
                  {
                    ii = iloc + 1; 
                    gx = -gx;
                  }
                  
                  if (not flags_array(iloc, jloc-1, kloc).isCovered() and
                      not flags_array(ii  , jloc-1, kloc).isCovered())
                  {
                    jj = jloc - 1;
                  }
                  else
                  {
                    jj = jloc + 1; 
                    gy = -gy;
                  }
                  
                  if (not flags_array(iloc, jloc, kloc-1).isCovered() and
                      not flags_array(ii  , jloc, kloc-1).isCovered() and
                      not flags_array(iloc, jj  , kloc-1).isCovered() and
                      not flags_array(ii  , jj  , kloc-1).isCovered())
                  {
                    kk = kloc - 1;
                  }
                  else
                  {
                    kk = kloc + 1; 
                    gz = -gz;
                  }

                  Real gxy = gx*gy;
                  Real gxz = gx*gz;
                  Real gyz = gy*gz;
                  Real gxyz = gx*gy*gz;

                  for (int n = 0; n < 3; n++)
                  {
                    velfp[n] = (1.0+gx+gy+gz+gxy+gxz+gyz+gxyz) * vel_array(iloc,jloc,kloc,n)
                             + (-gz - gxz - gyz - gxyz)        * vel_array(iloc,jloc,kk  ,n)
                             + (-gy - gxy - gyz - gxyz)        * vel_array(iloc,jj  ,kloc,n)
                             + (gyz + gxyz)                    * vel_array(iloc,jj  ,kk  ,n) 
                             + (-gx - gxy - gxz - gxyz)        * vel_array(ii  ,jloc,kloc,n)
                             + (gxz + gxyz)                    * vel_array(ii  ,jloc,kk  ,n)
                             + (gxy + gxyz)                    * vel_array(ii  ,jj  ,kloc,n)
                             + (-gxyz)                         * vel_array(ii  ,jj  ,kk  ,n);

                    // Keep the interpolated velocity between the cell value and
                    // the wall value (0)
                    if (velfp[n] > 0.0 && velfp[n] > vel_array(iloc,jloc,kloc,n))
                      velfp[n] = vel_array(iloc,jloc,kloc,n);
                    if (velfp[n] < 0.0 && velfp[n] < vel_array(iloc,jloc,kloc,n))
                      velfp[n] = vel_array(iloc,jloc,kloc,n);
                  }

                } // Cut cell

                Real  ep = ep_array(iloc,jloc,kloc);
                Real  ro = ro_array(iloc,jloc,kloc);
                Real  mu = mu_array(iloc,jloc,kloc);

                Real rad = particle.rdata(realData::radius);
                Real vol = particle.rdata(realData::volume);
     
                int p_id = particle.id();

                Real pvel[3]; 
                pvel[0] = particle.rdata(realData::velx);
                pvel[1] = particle.rdata(realData::vely);
                pvel[2] = particle.rdata(realData::velz);

                Real rop_g = ro * ep;
                        
                Real vslp[3];
                vslp[0] = velfp[0] - pvel[0];
                vslp[1] = velfp[1] - pvel[1];
                vslp[2] = velfp[2] - pvel[2];
                
                Real vrel = sqrt(dot_product(vslp, vslp));
                Real dpm = 2.0*rad;
                Real phis = 1.0 - ep;
            
                Real beta = vol*DragFunc(ep, mu, rop_g, vrel, dpm, dpm, phis,
                                         velfp[0], velfp[1], velfp[2],
                                         iloc, jloc, kloc, p_id); 
                particle.rdata(realData::dragx) = beta;

              } // Not covered
            }); // ip
          } // type of FAB
        } // if entire FAB not covered
      } // pti
    } // GPU region

    if (not OnSameGrids) {
      delete ep_ptr;
      delete ro_ptr;
      delete mu_ptr;
      delete vel_ptr;
    }
  } // lev
}
