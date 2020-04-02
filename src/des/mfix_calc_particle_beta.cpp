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

  // We copy the value inside the domain to the outside to avoid
  // unphysical volume fractions.
  const int dir_bc_in = 2;
  mfix_set_epg_bcs(get_ep_g(), dir_bc_in);

  // This is just a sanity check to make sure we're not using covered values
  // We can remove these lines once we're confident in the algorithm
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

      // Ghost cells needed for interpolation
      int ng = m_leveldata[lev]->ep_g->nGrow();
      ep_ptr = new MultiFab(pba, pdm, m_leveldata[lev]->ep_g->nComp(), ng);
      ep_ptr->copy(*m_leveldata[lev]->ep_g, 0, 0, 1, ng, ng);

      // Temporary arrays  -- copies with no ghost cells
      ro_ptr = new MultiFab(pba, pdm, m_leveldata[lev]->ro_g->nComp(), 1);
      ro_ptr->copy(*m_leveldata[lev]->ro_g, 0, 0, 1, 0, 0);

      mu_ptr = new MultiFab(pba, pdm, m_leveldata[lev]->mu_g->nComp(), 1);
      mu_ptr->copy(*m_leveldata[lev]->mu_g, 0, 0, 1, 0, 0);

      EBFArrayBoxFactory ebfactory_loc(*eb_levels[lev], geom[lev], pba, pdm,
                                       {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                                       EBSupport::volume);

      ng = m_leveldata[lev]->vel_g->nGrow();
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
      const auto dx_array  = geom[lev].CellSizeArray();
      const auto plo_array = geom[lev].ProbLoArray();

      const amrex::RealVect  dx( dx_array[0],  dx_array[1],  dx_array[2]);
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

              Real ep;
              trilinear_interp_scalar(particle.pos(), ep, ep_array, plo, dxi);

              // Indices of cell where particle is located
              int iloc = floor((particle.pos(0) - plo[0])*dxi[0]);
              int jloc = floor((particle.pos(1) - plo[1])*dxi[1]);
              int kloc = floor((particle.pos(2) - plo[2])*dxi[2]);

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

            const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(vel_ptr->Factory());

            const auto cellcent = &(factory.getCentroid());
            const auto& ccent_fab = cellcent->array(pti);

            const auto bndrycent = &(factory.getBndryCent());
            const auto& bcent_fab = bndrycent->array(pti);

            amrex::ParallelFor(np,
              [particles_ptr,vel_array,ro_array,mu_array,ep_array,flags_array,DragFunc,
              plo,dx,dxi,ccent_fab, bcent_fab]
              AMREX_GPU_DEVICE (int pid) noexcept
            {
              MFIXParticleContainer::ParticleType& particle = particles_ptr[pid];

              Real velfp[3];
              Real ep;

              // Cell containing particle centroid
              int ip = floor((particle.pos(0) - plo[0])*dxi[0]);
              int jp = floor((particle.pos(1) - plo[1])*dxi[1]);
              int kp = floor((particle.pos(2) - plo[2])*dxi[2]);

              // No drag force for particles in covered cells.
              if (flags_array(ip,jp,kp).isCovered() ){

                particle.rdata(realData::dragx) = 0.;

              // Cut or regular cell and none of the cells in the stencil is
              // covered (Note we can't assume regular cell has no covered
              // cells in the stencil because of the diagonal case)
              } else {

                // Upper cell in trilinear stencil
                int i = std::floor((particle.pos(0) - plo[0])*dxi[0] + 0.5);
                int j = std::floor((particle.pos(1) - plo[1])*dxi[1] + 0.5);
                int k = std::floor((particle.pos(2) - plo[2])*dxi[2] + 0.5);

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

                  trilinear_interp(particle.pos(), &velfp[0], vel_array, plo, dxi);
                  trilinear_interp_scalar(particle.pos(), ep, ep_array, plo, dxi);

                  // At least one of the cells in the stencil is cut or covered
                } else {

                  // Particle position relative to cell center [-0.5, 0.5]
                  Real gx = particle.pos(0)*dxi[0] - (ip + 0.5);
                  Real gy = particle.pos(1)*dxi[1] - (jp + 0.5);
                  Real gz = particle.pos(2)*dxi[2] - (kp + 0.5);

                  // Use the centoid location of the cell containing the particle
                  // to determine the interpolation stencil. If the particle is
                  // on the low side, then the high side stencil is the index of the
                  // cell contianing the particle, otherwise the particle is in
                  // the low side cell.
                  i = (gx < ccent_fab(ip,jp,kp,0)) ? ip : ip + 1;
                  j = (gy < ccent_fab(ip,jp,kp,1)) ? jp : jp + 1;
                  k = (gz < ccent_fab(ip,jp,kp,2)) ? kp : kp + 1;

                  amrex::Real nodes[8][3];
                  amrex::Real values[8][4];

                  int di = i - ip; // -1 or 0
                  int dj = j - jp; // -1 or 0
                  int dk = k - kp; // -1 or 0

                  // Count the number of non-conntected cells in the stencil
                  int covered = 0;
                  for(int kk(-1); kk<1; kk++){
                    for(int jj(-1); jj<1; jj++){
                      for(int ii(-1); ii<1; ii++){
                        if(not flags_array(ip,jp,kp).isConnected(di+ii,dj+jj,dk+kk))
                          covered += 1;
                      }
                    }
                  }

                  amrex::Print() <<  std::endl <<  std::endl;
                  amrex::Print() << "Particle Index: " << ip << " " << jp << " "<< kp << std::endl;
                  amrex::Print() << "Stencil Index:  " << i  << " " << j  << " "<< k  << std::endl;
                  amrex::Print() << "Covered:        " << covered << std::endl;

                  /*----------------------------------------------------------------------------------*
                   *                                                                                  *
                   *                                    NODE 0                                        *
                   *                                                                                  *
                   *----------------------------------------------------------------------------------*/
                  // Node 0
                  if(flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk-1)) {
                    //  DEBUG START ////////////////////////////////////////////////////////////////////
                    amrex::Print() << std::endl;
                    amrex::Print() << "Node 0 is connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i-1 << " " <<  j-1 << " "<<  k-1 << " " << std::endl;
                    amrex::Print() << "REF:    " << di-1 << " " << dj-1 << " "<< dk-1 << " " << std::endl;
                    amrex::Print() << "Node is regular?  "  << flags_array(i-1,j-1,k-1).isRegular() << std::endl;
                    amrex::Print() << "Node is svalued?  "  << flags_array(i-1,j-1,k-1).isSingleValued() << std::endl;
                    amrex::Print() << "Node is covered?  "  << flags_array(i-1,j-1,k-1).isCovered() << std::endl;
                    amrex::Print() << "Cell Centroids:   "
                                   << (i - 0.5 + ccent_fab(i-1,j-1,k-1,0))*dx[0] << " "
                                   << (j - 0.5 + ccent_fab(i-1,j-1,k-1,1))*dx[1] << " "
                                   << (k - 0.5 + ccent_fab(i-1,j-1,k-1,2))*dx[2] << std::endl;
                    if(flags_array(i-1,j-1,k-1).isSingleValued()){
                      amrex::Print() << "Bndry Centroids:  "
                                   << (i - 0.5 + bcent_fab(i-1,j-1,k-1,0))*dx[0] << " "
                                   << (j - 0.5 + bcent_fab(i-1,j-1,k-1,1))*dx[1] << " "
                                   << (k - 0.5 + bcent_fab(i-1,j-1,k-1,2))*dx[2] << std::endl;
                    }
                    //  DEBUG END   ////////////////////////////////////////////////////////////////////
                    nodes[0][0] = (i - 0.5 + ccent_fab(i-1,j-1,k-1,0))*dx[0];
                    nodes[0][1] = (j - 0.5 + ccent_fab(i-1,j-1,k-1,1))*dx[1];
                    nodes[0][2] = (k - 0.5 + ccent_fab(i-1,j-1,k-1,2))*dx[2];

                    values[0][0] = vel_array(i-1,j-1,k-1,0);
                    values[0][1] = vel_array(i-1,j-1,k-1,1);
                    values[0][2] = vel_array(i-1,j-1,k-1,2);
                    values[0][3] =  ep_array(i-1,j-1,k-1);

                  } else {
                    amrex::Print() << std::endl << "Node 0 is not connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i-1 << " " <<  j-1 << " "<<  k-1 << " " << std::endl;
                    amrex::Print() << "REF:    " << di-1 << " " << dj-1 << " "<< dk-1 << " " << std::endl;

                    int ib = ip;
                    int jb = jp;
                    int kb = kp;

                    if(covered == 2) {
                      if(flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk-1)) {
                        // Node 1 is covered --> use Node 7
                        amrex::Print() << "Using EB in Node 7." << std::endl;
                        ib = i-1;  jb = j  ;  kb = k  ;
                      } else if(flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk-1)) {
                        // Node 3 is covered --> use Node 5
                        amrex::Print() << "Using EB in Node 5." << std::endl;
                        ib = i  ;  jb = j-1;  kb = k  ;
                      } else if(flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk  )) {
                        // Node 4 is covered --> use Node 2
                        amrex::Print() << "Using EB in Node 2." << std::endl;
                        ib = i  ;  jb = j  ;  kb = k-1;
                      }
                    }
                    nodes[0][0] = (ib + 0.5 + bcent_fab(ib,jb,kb,0))*dx[0];
                    nodes[0][1] = (jb + 0.5 + bcent_fab(ib,jb,kb,1))*dx[1];
                    nodes[0][2] = (kb + 0.5 + bcent_fab(ib,jb,kb,2))*dx[2];

                    values[0][0] = 0.;
                    values[0][1] = 0.;
                    values[0][2] = 0.;
                    values[0][3] = ep_array(ib,jb,kb);
                  }

                  /*----------------------------------------------------------------------------------*
                   *                                                                                  *
                   *                                    NODE 1                                        *
                   *                                                                                  *
                   *----------------------------------------------------------------------------------*/
                  if(flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk-1)) {
                    //  DEBUG START ////////////////////////////////////////////////////////////////////
                    amrex::Print() << std::endl;
                    amrex::Print() << "Node 1 is connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i   << " " <<  j-1 << " "<<  k-1 << " " << std::endl;
                    amrex::Print() << "REF:    " << di   << " " << dj-1 << " "<< dk-1 << " " << std::endl;
                    amrex::Print() << "Node is regular?  "  << flags_array(i  ,j-1,k-1).isRegular() << std::endl;
                    amrex::Print() << "Node is svalued?  "  << flags_array(i  ,j-1,k-1).isSingleValued() << std::endl;
                    amrex::Print() << "Node is covered?  "  << flags_array(i  ,j-1,k-1).isCovered() << std::endl;
                    amrex::Print() << "Cell Centroids:   "
                                   << (i + 0.5 + ccent_fab(i  ,j-1,k-1,0))*dx[0] << " "
                                   << (j - 0.5 + ccent_fab(i  ,j-1,k-1,1))*dx[1] << " "
                                   << (k - 0.5 + ccent_fab(i  ,j-1,k-1,2))*dx[2] << std::endl;
                    if(flags_array(i-1,j-1,k-1).isSingleValued()){
                      amrex::Print() << "Bndry Centroids:  "
                                   << (i + 0.5 + bcent_fab(i  ,j-1,k-1,0))*dx[0] << " "
                                   << (j - 0.5 + bcent_fab(i  ,j-1,k-1,1))*dx[1] << " "
                                   << (k - 0.5 + bcent_fab(i  ,j-1,k-1,2))*dx[2] << std::endl;
                    }
                    //  DEBUG END   ////////////////////////////////////////////////////////////////////

                    nodes[1][0] = (i + 0.5 + ccent_fab(i  ,j-1,k-1,0))*dx[0];
                    nodes[1][1] = (j - 0.5 + ccent_fab(i  ,j-1,k-1,1))*dx[1];
                    nodes[1][2] = (k - 0.5 + ccent_fab(i  ,j-1,k-1,2))*dx[2];

                    values[1][0] = vel_array(i  ,j-1,k-1,0);
                    values[1][1] = vel_array(i  ,j-1,k-1,1);
                    values[1][2] = vel_array(i  ,j-1,k-1,2);
                    values[1][3] =  ep_array(i  ,j-1,k-1);
                    
                  } else {
                    amrex::Print() << std::endl << "Node is not connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i   << " " <<  j-1 << " "<<  k-1 << " " << std::endl;
                    amrex::Print() << "REF:    " << di   << " " << dj-1 << " "<< dk-1 << " " << std::endl;

                    int ib = ip;
                    int jb = jp;
                    int kb = kp;

                    if(covered == 2) {
                      if(flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk-1)) {
                        // Node 0 is covered --> use Node 6
                        amrex::Print() << "Using EB in Node 6." << std::endl;
                        ib = i  ;  jb = j  ;  kb = k  ;
                      } else if(flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk-1)) {
                        // Node 2 is covered --> use Node 4
                        amrex::Print() << "Using EB in Node 4." << std::endl;
                        ib = i-1;  jb = j-1;  kb = k  ;
                      } else if(flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk  )) {
                        // Node 5 is covered --> use Node 3
                        amrex::Print() << "Using EB in Node 3." << std::endl;
                        ib = i-1;  jb = j  ;  kb = k-1;
                      }
                    }
                    nodes[0][0] = (ib + 0.5 + bcent_fab(ib,jb,kb,0))*dx[0];
                    nodes[0][1] = (jb + 0.5 + bcent_fab(ib,jb,kb,1))*dx[1];
                    nodes[0][2] = (kb + 0.5 + bcent_fab(ib,jb,kb,2))*dx[2];

                    values[0][0] = 0.;
                    values[0][1] = 0.;
                    values[0][2] = 0.;
                    values[0][3] = ep_array(ib,jb,kb);
                  }

                  /*----------------------------------------------------------------------------------*
                   *                                                                                  *
                   *                                    NODE 2                                        *
                   *                                                                                  *
                   *----------------------------------------------------------------------------------*/
                  if(flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk-1)) {
                    //  DEBUG START ////////////////////////////////////////////////////////////////////
                    amrex::Print() << std::endl;
                    amrex::Print() << "Node 2 is connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i   << " " <<  j   << " "<<  k-1 << " " << std::endl;
                    amrex::Print() << "REF:    " << di   << " " << dj   << " "<< dk-1 << " " << std::endl;
                    amrex::Print() << "Node is regular?  "  << flags_array(i  ,j  ,k-1).isRegular() << std::endl;
                    amrex::Print() << "Node is svalued?  "  << flags_array(i  ,j  ,k-1).isSingleValued() << std::endl;
                    amrex::Print() << "Node is covered?  "  << flags_array(i  ,j  ,k-1).isCovered() << std::endl;
                    amrex::Print() << "Cell Centroids:   "
                                   << (i + 0.5 + ccent_fab(i  ,j  ,k-1,0))*dx[0] << " "
                                   << (j + 0.5 + ccent_fab(i  ,j  ,k-1,1))*dx[1] << " "
                                   << (k - 0.5 + ccent_fab(i  ,j  ,k-1,2))*dx[2] << std::endl;
                    if(flags_array(i  ,j  ,k-1).isSingleValued()){
                      amrex::Print() << "Bndry Centroids:  "
                                   << (i + 0.5 + bcent_fab(i  ,j  ,k-1,0))*dx[0] << " "
                                   << (j + 0.5 + bcent_fab(i  ,j  ,k-1,1))*dx[1] << " "
                                   << (k - 0.5 + bcent_fab(i  ,j  ,k-1,2))*dx[2] << std::endl;
                    }
                    //  DEBUG END   ////////////////////////////////////////////////////////////////////


                    nodes[2][0] = (i + 0.5 + ccent_fab(i  ,j  ,k-1,0))*dx[0];
                    nodes[2][1] = (j + 0.5 + ccent_fab(i  ,j  ,k-1,1))*dx[1];
                    nodes[2][2] = (k - 0.5 + ccent_fab(i  ,j  ,k-1,2))*dx[2];

                    values[2][0] = vel_array(i  ,j  ,k-1,0);
                    values[2][1] = vel_array(i  ,j  ,k-1,1);
                    values[2][2] = vel_array(i  ,j  ,k-1,2);
                    values[2][3] =  ep_array(i  ,j  ,k-1);
                    
                  } else {
                    amrex::Print() << std::endl << "Node 2 is not connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i   << " " <<  j   << " "<<  k-1 << " " << std::endl;
                    amrex::Print() << "REF:    " << di   << " " << dj   << " "<< dk-1 << " " << std::endl;
                    int ib = ip;
                    int jb = jp;
                    int kb = kp;

                    if(covered == 2) {
                      if(flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk-1)) {
                        // Node 1 is covered --> use Node 7
                        amrex::Print() << "Using EB in Node 7." << std::endl;
                        ib = i-1;  jb = j  ;  kb = k  ;
                      } else if(flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk-1)) {
                        // Node 3 is covered --> use Node 5
                        amrex::Print() << "Using EB in Node 5." << std::endl;
                        ib = i  ;  jb = j-1;  kb = k  ;
                      } else if(flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk  )) {
                        // Node 6 is covered --> use Node 0
                        amrex::Print() << "Using EB in Node 0." << std::endl;
                        ib = i-1;  jb = j-1;  kb = k-1;
                      }
                    }
                    nodes[0][0] = (ib + 0.5 + bcent_fab(ib,jb,kb,0))*dx[0];
                    nodes[0][1] = (jb + 0.5 + bcent_fab(ib,jb,kb,1))*dx[1];
                    nodes[0][2] = (kb + 0.5 + bcent_fab(ib,jb,kb,2))*dx[2];

                    values[0][0] = 0.;
                    values[0][1] = 0.;
                    values[0][2] = 0.;
                    values[0][3] = ep_array(ib,jb,kb);
                  }

                  /*----------------------------------------------------------------------------------*
                   *                                                                                  *
                   *                                    NODE 3                                        *
                   *                                                                                  *
                   *----------------------------------------------------------------------------------*/
                  if(flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk-1)) {
                    //  DEBUG START ////////////////////////////////////////////////////////////////////
                    amrex::Print() << std::endl;
                    amrex::Print() << "Node 3 is connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i-1 << " " <<  j   << " "<<  k-1 << " " << std::endl;
                    amrex::Print() << "REF:    " << di-1 << " " << dj   << " "<< dk-1 << " " << std::endl;
                    amrex::Print() << "Node is regular?  "  << flags_array(i-1,j  ,k-1).isRegular() << std::endl;
                    amrex::Print() << "Node is svalued?  "  << flags_array(i-1,j  ,k-1).isSingleValued() << std::endl;
                    amrex::Print() << "Node is covered?  "  << flags_array(i-1,j  ,k-1).isCovered() << std::endl;
                    amrex::Print() << "Cell Centroids:   "
                                   << (i - 0.5 + ccent_fab(i-1,j  ,k-1,0))*dx[0] << " "
                                   << (j + 0.5 + ccent_fab(i-1,j  ,k-1,1))*dx[1] << " "
                                   << (k - 0.5 + ccent_fab(i-1,j  ,k-1,2))*dx[2] << std::endl;
                    if(flags_array(i-1,j  ,k-1).isSingleValued()){
                      amrex::Print() << "Bndry Centroids:  "
                                   << (i - 0.5 + bcent_fab(i-1,j  ,k-1,0))*dx[0] << " "
                                   << (j + 0.5 + bcent_fab(i-1,j  ,k-1,1))*dx[1] << " "
                                   << (k - 0.5 + bcent_fab(i-1,j  ,k-1,2))*dx[2] << std::endl;
                    }
                    //  DEBUG END   ////////////////////////////////////////////////////////////////////
                    nodes[3][0] = (i - 0.5 + ccent_fab(i-1,j  ,k-1,0))*dx[0];
                    nodes[3][1] = (j + 0.5 + ccent_fab(i-1,j  ,k-1,1))*dx[1];
                    nodes[3][2] = (k - 0.5 + ccent_fab(i-1,j  ,k-1,2))*dx[2];

                    values[3][0] = vel_array(i-1,j  ,k-1,0);
                    values[3][1] = vel_array(i-1,j  ,k-1,1);
                    values[3][2] = vel_array(i-1,j  ,k-1,2);
                    values[3][3] =  ep_array(i-1,j  ,k-1);
                    
                  } else {
                    amrex::Print() << std::endl << "Node 3 is not connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i-1 << " " <<  j   << " "<<  k-1 << " " << std::endl;
                    amrex::Print() << "REF:    " << di-1 << " " << dj   << " "<< dk-1 << " " << std::endl;
                    int ib = ip;
                    int jb = jp;
                    int kb = kp;

                    if(covered == 2) {
                      if(flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk-1)) {
                        // Node 0 is covered --> use Node 6
                        amrex::Print() << "Using EB in Node 6." << std::endl;
                        ib = i  ;  jb = j  ;  kb = k  ;
                      } else if(flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk-1)) {
                        // Node 2 is covered --> use Node 4
                        amrex::Print() << "Using EB in Node 4." << std::endl;
                        ib = i-1;  jb = j-1;  kb = k  ;
                      } else if(flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk  )) {
                        // Node 7 is covered --> use Node 1
                        amrex::Print() << "Using EB in Node 1." << std::endl;
                        ib = i  ;  jb = j-1;  kb = k-1;
                      }
                    }
                    nodes[0][0] = (ib + 0.5 + bcent_fab(ib,jb,kb,0))*dx[0];
                    nodes[0][1] = (jb + 0.5 + bcent_fab(ib,jb,kb,1))*dx[1];
                    nodes[0][2] = (kb + 0.5 + bcent_fab(ib,jb,kb,2))*dx[2];

                    values[0][0] = 0.;
                    values[0][1] = 0.;
                    values[0][2] = 0.;
                    values[0][3] = ep_array(ib,jb,kb);
                  }

                  /*----------------------------------------------------------------------------------*
                   *                                                                                  *
                   *                                    NODE 4                                        *
                   *                                                                                  *
                   *----------------------------------------------------------------------------------*/
                  if(flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk  )) {
                    //  DEBUG START ////////////////////////////////////////////////////////////////////
                    amrex::Print() << std::endl;
                    amrex::Print() << "Node 4 is connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i-1 << " " <<  j-1 << " "<<  k   << " " << std::endl;
                    amrex::Print() << "REF:    " << di-1 << " " << dj-1 << " "<< dk   << " " << std::endl;
                    amrex::Print() << "Node is regular?  "  << flags_array(i-1,j-1,k  ).isRegular() << std::endl;
                    amrex::Print() << "Node is svalued?  "  << flags_array(i-1,j-1,k  ).isSingleValued() << std::endl;
                    amrex::Print() << "Node is covered?  "  << flags_array(i-1,j-1,k  ).isCovered() << std::endl;
                    amrex::Print() << "Cell Centroids:   "
                                   << (i - 0.5 + ccent_fab(i-1,j-1,k  ,0))*dx[0] << " "
                                   << (j - 0.5 + ccent_fab(i-1,j-1,k  ,1))*dx[1] << " "
                                   << (k + 0.5 + ccent_fab(i-1,j-1,k  ,2))*dx[2] << std::endl;
                    if(flags_array(i-1,j-1,k  ).isSingleValued()){
                      amrex::Print() << "Bndry Centroids:  "
                                   << (i - 0.5 + bcent_fab(i-1,j-1,k  ,0))*dx[0] << " "
                                   << (j - 0.5 + bcent_fab(i-1,j-1,k  ,1))*dx[1] << " "
                                   << (k + 0.5 + bcent_fab(i-1,j-1,k  ,2))*dx[2] << std::endl;
                    }
                    //  DEBUG END   ////////////////////////////////////////////////////////////////////
                    nodes[4][0] = (i - 0.5 + ccent_fab(i-1,j-1,k  ,0))*dx[0];
                    nodes[4][1] = (j - 0.5 + ccent_fab(i-1,j-1,k  ,1))*dx[1];
                    nodes[4][2] = (k + 0.5 + ccent_fab(i-1,j-1,k  ,2))*dx[2];

                    values[4][0] = vel_array(i-1,j-1,k  ,0);
                    values[4][1] = vel_array(i-1,j-1,k  ,1);
                    values[4][2] = vel_array(i-1,j-1,k  ,2);
                    values[4][3] =  ep_array(i-1,j-1,k  );

                  } else {
                    amrex::Print() << std::endl << "Node 4 is not connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i-1 << " " <<  j-1 << " "<<  k   << " " << std::endl;
                    amrex::Print() << "REF:    " << di-1 << " " << dj-1 << " "<< dk   << " " << std::endl;
                    int ib = ip;
                    int jb = jp;
                    int kb = kp;

                    if(covered == 2) {
                      if(flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk-1)) {
                        // Node 0 is covered --> use Node 6
                        amrex::Print() << "Using EB in Node 6." << std::endl;
                        ib = i  ;  jb = j  ;  kb = k  ;
                      } else if(flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk  )) {
                        // Node 5 is covered --> use Node 3
                        amrex::Print() << "Using EB in Node 3." << std::endl;
                        ib = i-1;  jb = j  ;  kb = k-1;
                      } else if(flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk  )) {
                        // Node 7 is covered --> use Node 1
                        amrex::Print() << "Using EB in Node 1." << std::endl;
                        ib = i  ;  jb = j-1;  kb = k-1;
                      }
                    }
                    nodes[0][0] = (ib + 0.5 + bcent_fab(ib,jb,kb,0))*dx[0];
                    nodes[0][1] = (jb + 0.5 + bcent_fab(ib,jb,kb,1))*dx[1];
                    nodes[0][2] = (kb + 0.5 + bcent_fab(ib,jb,kb,2))*dx[2];

                    values[0][0] = 0.;
                    values[0][1] = 0.;
                    values[0][2] = 0.;
                    values[0][3] = ep_array(ib,jb,kb);
                  }

                  /*----------------------------------------------------------------------------------*
                   *                                                                                  *
                   *                                    NODE 5                                        *
                   *                                                                                  *
                   *----------------------------------------------------------------------------------*/
                  if(flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk  )) {
                    //  DEBUG START ////////////////////////////////////////////////////////////////////
                    amrex::Print() << std::endl;
                    amrex::Print() << "Node 5 is connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i   << " " <<  j-1 << " "<<  k   << " " << std::endl;
                    amrex::Print() << "REF:    " << di   << " " << dj-1 << " "<< dk   << " " << std::endl;
                    amrex::Print() << "Node is regular?  "  << flags_array(i  ,j-1,k  ).isRegular() << std::endl;
                    amrex::Print() << "Node is svalued?  "  << flags_array(i  ,j-1,k  ).isSingleValued() << std::endl;
                    amrex::Print() << "Node is covered?  "  << flags_array(i  ,j-1,k  ).isCovered() << std::endl;
                    amrex::Print() << "Cell Centroids:   "
                                   << (i - 0.5 + ccent_fab(i  ,j-1,k  ,0))*dx[0] << " "
                                   << (j - 0.5 + ccent_fab(i  ,j-1,k  ,1))*dx[1] << " "
                                   << (k + 0.5 + ccent_fab(i  ,j-1,k  ,2))*dx[2] << std::endl;
                    if(flags_array(i  ,j-1,k  ).isSingleValued()){
                      amrex::Print() << "Bndry Centroids:  "
                                   << (i - 0.5 + bcent_fab(i  ,j-1,k  ,0))*dx[0] << " "
                                   << (j - 0.5 + bcent_fab(i  ,j-1,k  ,1))*dx[1] << " "
                                   << (k + 0.5 + bcent_fab(i  ,j-1,k  ,2))*dx[2] << std::endl;
                    }
                    //  DEBUG END   ////////////////////////////////////////////////////////////////////
                    nodes[5][0] = (i + 0.5 + ccent_fab(i  ,j-1,k  ,0))*dx[0];
                    nodes[5][1] = (j - 0.5 + ccent_fab(i  ,j-1,k  ,1))*dx[1];
                    nodes[5][2] = (k + 0.5 + ccent_fab(i  ,j-1,k  ,2))*dx[2];

                    values[5][0] = vel_array(i  ,j-1,k  ,0);
                    values[5][1] = vel_array(i  ,j-1,k  ,1);
                    values[5][2] = vel_array(i  ,j-1,k  ,2);
                    values[5][3] =  ep_array(i  ,j-1,k  );
                    
                  } else {
                    amrex::Print() << std::endl << "Node 5 is not connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i   << " " <<  j-1 << " "<<  k   << " " << std::endl;
                    amrex::Print() << "REF:    " << di   << " " << dj-1 << " "<< dk   << " " << std::endl;
                    int ib = ip;
                    int jb = jp;
                    int kb = kp;

                    if(covered == 2) {
                      if(flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk-1)) {
                        // Node 1 is covered --> use Node 7
                        amrex::Print() << "Using EB in Node 7." << std::endl;
                        ib = i-1;  jb = j  ;  kb = k  ;
                      } else if(flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk  )) {
                        // Node 4 is covered --> use Node 2
                        amrex::Print() << "Using EB in Node 2." << std::endl;
                        ib = i  ;  jb = j  ;  kb = k-1;
                      } else if(flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk  )) {
                        // Node 6 is covered --> use Node 0
                        amrex::Print() << "Using EB in Node 0." << std::endl;
                        ib = i-1;  jb = j-1;  kb = k-1;
                      }
                    }
                    nodes[0][0] = (ib + 0.5 + bcent_fab(ib,jb,kb,0))*dx[0];
                    nodes[0][1] = (jb + 0.5 + bcent_fab(ib,jb,kb,1))*dx[1];
                    nodes[0][2] = (kb + 0.5 + bcent_fab(ib,jb,kb,2))*dx[2];

                    values[0][0] = 0.;
                    values[0][1] = 0.;
                    values[0][2] = 0.;
                    values[0][3] = ep_array(ib,jb,kb);
                  }

                  /*----------------------------------------------------------------------------------*
                   *                                                                                  *
                   *                                    NODE 6                                        *
                   *                                                                                  *
                   *----------------------------------------------------------------------------------*/
                  if(flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk  )) {
                    //  DEBUG START ////////////////////////////////////////////////////////////////////
                    amrex::Print() << std::endl;
                    amrex::Print() << "Node 6 is connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i   << " " <<  j   << " "<<  k   << " " << std::endl;
                    amrex::Print() << "REF:    " << di   << " " << dj   << " "<< dk   << " " << std::endl;
                    amrex::Print() << "Node is regular?  "  << flags_array(i  ,j  ,k  ).isRegular() << std::endl;
                    amrex::Print() << "Node is svalued?  "  << flags_array(i  ,j  ,k  ).isSingleValued() << std::endl;
                    amrex::Print() << "Node is covered?  "  << flags_array(i  ,j  ,k  ).isCovered() << std::endl;
                    amrex::Print() << "Cell Centroids:   "
                                   << (i + 0.5 + ccent_fab(i  ,j  ,k  ,0))*dx[0] << " "
                                   << (j + 0.5 + ccent_fab(i  ,j  ,k  ,1))*dx[1] << " "
                                   << (k + 0.5 + ccent_fab(i  ,j  ,k  ,2))*dx[2] << std::endl;
                    if(flags_array(i  ,j  ,k  ).isSingleValued()){
                      amrex::Print() << "Bndry Centroids:  "
                                   << (i + 0.5 + bcent_fab(i  ,j  ,k  ,0))*dx[0] << " "
                                   << (j + 0.5 + bcent_fab(i  ,j  ,k  ,1))*dx[1] << " "
                                   << (k + 0.5 + bcent_fab(i  ,j  ,k  ,2))*dx[2] << std::endl;
                    }
                    //  DEBUG END   ////////////////////////////////////////////////////////////////////
                    nodes[6][0] = (i + 0.5 + ccent_fab(i  ,j  ,k  ,0))*dx[0];
                    nodes[6][1] = (j + 0.5 + ccent_fab(i  ,j  ,k  ,1))*dx[1];
                    nodes[6][2] = (k + 0.5 + ccent_fab(i  ,j  ,k  ,2))*dx[2];

                    values[6][0] = vel_array(i  ,j  ,k  ,0);
                    values[6][1] = vel_array(i  ,j  ,k  ,1);
                    values[6][2] = vel_array(i  ,j  ,k  ,2);
                    values[6][3] =  ep_array(i  ,j  ,k  );
                    
                  } else {
                    amrex::Print() << std::endl << "Node 6 is not connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i   << " " <<  j   << " "<<  k   << " " << std::endl;
                    amrex::Print() << "REF:    " << di   << " " << dj   << " "<< dk   << " " << std::endl;

                    int ib = ip;
                    int jb = jp;
                    int kb = kp;

                    if(covered == 2) {
                      if(not flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk-1)) {
                        // Node 2 is covered --> use Node 4
                        amrex::Print() << "Using EB in Node 4." << std::endl;
                        ib = i-1;  jb = j-1;  kb = k  ;
                      } else if(not flags_array(ip,jp,kp).isConnected(di  ,dj-1,dk  )) {
                        // Node 5 is covered --> use Node 3
                        amrex::Print() << "Using EB in Node 3." << std::endl;
                        ib = i-1;  jb = j  ;  kb = k-1;
                      } else if(not flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk  )) {
                        // Node 7 is covered --> use Node 1
                        amrex::Print() << "Using EB in Node 1." << std::endl;
                        ib = i  ;  jb = j-1;  kb = k-1;
                      }
                    }

                    amrex::Print() << "Index:  " <<  ib  << " " <<  jb  << " "<<  kb  << " " << std::endl;

                    nodes[6][0] = (ib + 0.5 + bcent_fab(ib,jb,kb,0))*dx[0];
                    nodes[6][1] = (jb + 0.5 + bcent_fab(ib,jb,kb,1))*dx[1];
                    nodes[6][2] = (kb + 0.5 + bcent_fab(ib,jb,kb,2))*dx[2];

                    values[6][0] = 0.;
                    values[6][1] = 0.;
                    values[6][2] = 0.;
                    values[6][3] = ep_array(ib,jb,kb);
                  }

                  /*----------------------------------------------------------------------------------*
                   *                                                                                  *
                   *                                    NODE 7                                        *
                   *                                                                                  *
                   *----------------------------------------------------------------------------------*/
                  if(flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk  )) {
                    //  DEBUG START ////////////////////////////////////////////////////////////////////
                    amrex::Print() << std::endl;
                    amrex::Print() << "Node 7 is connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i-1 << " " <<  j   << " "<<  k   << " " << std::endl;
                    amrex::Print() << "REF:    " << di-1 << " " << dj   << " "<< dk   << " " << std::endl;
                    amrex::Print() << "Node is regular?  "  << flags_array(i-1,j  ,k  ).isRegular() << std::endl;
                    amrex::Print() << "Node is svalued?  "  << flags_array(i-1,j  ,k  ).isSingleValued() << std::endl;
                    amrex::Print() << "Node is covered?  "  << flags_array(i-1,j  ,k  ).isCovered() << std::endl;
                    amrex::Print() << "Cell Centroids:   "
                                   << (i - 0.5 + ccent_fab(i-1,j  ,k  ,0))*dx[0] << " "
                                   << (j + 0.5 + ccent_fab(i-1,j  ,k  ,1))*dx[1] << " "
                                   << (k + 0.5 + ccent_fab(i-1,j  ,k  ,2))*dx[2] << std::endl;
                    if(flags_array(i-1,j  ,k  ).isSingleValued()){
                      amrex::Print() << "Bndry Centroids:  "
                                   << (i - 0.5 + bcent_fab(i-1,j  ,k  ,0))*dx[0] << " "
                                   << (j + 0.5 + bcent_fab(i-1,j  ,k  ,1))*dx[1] << " "
                                   << (k + 0.5 + bcent_fab(i-1,j  ,k  ,2))*dx[2] << std::endl;
                    }
                    //  DEBUG END   ////////////////////////////////////////////////////////////////////
                    nodes[7][0] = (i - 0.5 + ccent_fab(i-1,j  ,k  ,0))*dx[0];
                    nodes[7][1] = (j + 0.5 + ccent_fab(i-1,j  ,k  ,1))*dx[1];
                    nodes[7][2] = (k + 0.5 + ccent_fab(i-1,j  ,k  ,2))*dx[2];

                    values[7][0] = vel_array(i-1,j  ,k  ,0);
                    values[7][1] = vel_array(i-1,j  ,k  ,1);
                    values[7][2] = vel_array(i-1,j  ,k  ,2);
                    values[7][3] =  ep_array(i-1,j  ,k  );
                    
                  } else {
                    amrex::Print() << std::endl << "Node 7 is not connected." << std::endl;
                    amrex::Print() << "Index:  " <<  i-1 << " " <<  j   << " "<<  k   << " " << std::endl;
                    amrex::Print() << "REF:    " << di-1 << " " << dj   << " "<< dk   << " " << std::endl;

                    int ib = ip;
                    int jb = jp;
                    int kb = kp;

                    if(covered == 2) {
                      if(flags_array(ip,jp,kp).isConnected(di-1,dj  ,dk-1)) {
                        // Node 3 is covered --> use Node 5
                        amrex::Print() << "Using EB in Node 5." << std::endl;
                        ib = i  ;  jb = j-1;  kb = k  ;
                      } else if(flags_array(ip,jp,kp).isConnected(di-1,dj-1,dk  )) {
                        // Node 4 is covered --> use Node 2
                        amrex::Print() << "Using EB in Node 2." << std::endl;
                        ib = i  ;  jb = j  ;  kb = k-1;
                      } else if(flags_array(ip,jp,kp).isConnected(di  ,dj  ,dk  )) {
                        // Node 6 is covered --> use Node 0
                        amrex::Print() << "Using EB in Node 6." << std::endl;
                        ib = i-1;  jb = j-1;  kb = k-1;
                      }
                    }
                    nodes[7][0] = (ib + 0.5 + bcent_fab(ib,jb,kb,0))*dx[0];
                    nodes[7][1] = (jb + 0.5 + bcent_fab(ib,jb,kb,1))*dx[1];
                    nodes[7][2] = (kb + 0.5 + bcent_fab(ib,jb,kb,2))*dx[2];

                    values[7][0] = 0.;
                    values[7][1] = 0.;
                    values[7][2] = 0.;
                    values[7][3] = ep_array(ib,jb,kb);

                  }

                    amrex::Print() << std::endl << std::endl;
                      for(int n(0); n<8; n++){
                      amrex::Print() << " node " << n << ": " 
                        << nodes[n][0] << " "
                        << nodes[n][1] << " "
                        << nodes[n][2] << std::endl;
                      }

                    amrex::Print() << std::endl << std::endl;
                      for(int n(0); n<8; n++){
                      amrex::Print() << " values " << n << ": " 
                        << values[n][0] << " "
                        << values[n][1] << " "
                        << values[n][2] << " "
                        << values[n][3] << std::endl;
                      }




                velfp[0] = vel_array(ip,jp,kp,0);
                velfp[1] = vel_array(ip,jp,kp,1);
                velfp[2] = vel_array(ip,jp,kp,2);

                ep = ep_array(ip,jp,kp);

                } // Cut cell

                // Using i/j/k of centroid cell
                Real  ro = ro_array(ip,jp,kp);
                Real  mu = mu_array(ip,jp,kp);

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
                                         ip, jp, kp, p_id);

                particle.rdata(realData::dragx) = beta;

              } // Not covered
            }); // pid
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


  // Reset the volume fractions back to the correct values at
  // inflow faces.
  const int dir_bc_out = 1;
  mfix_set_epg_bcs(get_ep_g(), dir_bc_out);

}
