#include <mfix.H>

#include <mfix_bc_parms.H>
#include <mfix_fluid_parms.H>


using namespace BC;

void
mfix::mfix_set_bc_type (int lev, int nghost_bc)
{
    if (ooo_debug) amrex::Print() << "mfix_set_bc_type" << std::endl;

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    const GpuArray<Real, 3> plo = geom[lev].ProbLoArray();

    const int und_  = bc_list.get_undefined();
    const int ig_   = bc_list.get_ig();

    const int l_species = fluid.nspecies;
    const int l_ntrac = ntrac;
    const int l_force = amrex::max(AMREX_SPACEDIM, l_ntrac, l_species);

    // Set the defaults for BCRecs
    m_bcrec_velocity.resize(AMREX_SPACEDIM);
    m_bcrec_hydro_velocity.resize(AMREX_SPACEDIM);
    m_bcrec_density.resize(1);
    m_bcrec_enthalpy.resize(1);
    m_bcrec_tracer.resize(l_ntrac);
    m_bcrec_species.resize(l_species);
    m_bcrec_force.resize(l_force);

    { // begin x direction

      const int dir = 0;

      Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
      Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();

      const int init_x = geom[lev].isPeriodic(0) ? und_ : ig_;

      Box domainx(geom[lev].Domain());
      domainx.grow(1,nghost_bc);  // Add ghost cells to y
      domainx.grow(2,nghost_bc);  // Add ghost cells to z

      { // x-lo side of the domain

        Box box_ilo = amrex::adjCellLo(domainx,0,1);

        IntVect ibx_lo(box_ilo.loVect());
        IntVect ibx_hi(box_ilo.hiVect());

        int xlo_type = init_x;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_ilo, [bc_ilo_type, init_x]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_ilo_type(i,j,k,0) = init_x;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_xlo.size(); ++lc) {

          const int bcv  = bc_xlo[lc];
          const int type = bc[bcv].type;

          xlo_type = type;

          if (lc > 0){
            ibx_lo[1] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(1)-plo[1])/dy + 0.5));
            ibx_lo[2] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(2)-plo[2])/dz + 0.5));

            ibx_hi[1] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(1)-plo[1])/dy + 0.5))-1;
            ibx_hi[2] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(2)-plo[2])/dz + 0.5))-1;
          }

          const Box lo_box(ibx_lo, ibx_hi);

          amrex::ParallelFor(lo_box, [bc_ilo_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_ilo_type(i,j,k,0) = type;
               bc_ilo_type(i,j,k,1) = bcv;
             });
        }

        set_bcrec_lo(lev, dir, xlo_type);

      } // end x-lo side of the domain


      { // x-hi side of the domain

        Box box_ihi = amrex::adjCellHi(domainx,0,1);

        IntVect ibx_lo(box_ihi.loVect());
        IntVect ibx_hi(box_ihi.hiVect());

        int xhi_type = init_x;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_ihi, [bc_ihi_type, init_x]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_ihi_type(i,j,k,0) = init_x;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_xhi.size(); ++lc) {

          const int bcv  = bc_xhi[lc];
          const int type = bc[bcv].type;

          xhi_type = type;

          if (lc > 0){
            ibx_lo[1] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(1)-plo[1])/dy + 0.5));
            ibx_lo[2] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(2)-plo[2])/dz + 0.5));

            ibx_hi[1] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(1)-plo[1])/dy + 0.5))-1;
            ibx_hi[2] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(2)-plo[2])/dz + 0.5))-1;
          }

          const Box hi_box(ibx_lo, ibx_hi);

          amrex::ParallelFor(hi_box, [bc_ihi_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_ihi_type(i,j,k,0) = type;
               bc_ihi_type(i,j,k,1) = bcv;
             });

        }

        set_bcrec_hi(lev, dir, xhi_type);
      } // end x-hi side of the domain
    } // end x-direction



    { // begin y-direction

      const int dir = 1;

      const int init_y = geom[lev].isPeriodic(1) ? und_ : ig_;

      Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
      Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();

      Box domainy(geom[lev].Domain());
      domainy.grow(0,nghost_bc);  // Add ghost cells to x
      domainy.grow(2,nghost_bc);  // Add ghost cells to z

      { // y-lo side of the domain

        Box box_jlo = amrex::adjCellLo(domainy,1,1);

        IntVect jbx_lo(box_jlo.loVect());
        IntVect jbx_hi(box_jlo.hiVect());

        int ylo_type = init_y;

        // Initialize y-lo domain extent.
        amrex::ParallelFor(box_jlo, [bc_jlo_type, init_y]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_jlo_type(i,j,k,0) = init_y;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_ylo.size(); ++lc) {

          const int bcv  = bc_ylo[lc];
          const int type = bc[bcv].type;

          ylo_type = type;

          if (lc > 0){
            jbx_lo[0] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(0)-plo[0])/dx + 0.5));
            jbx_lo[2] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(2)-plo[2])/dz + 0.5));

            jbx_hi[0] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(0)-plo[0])/dx + 0.5))-1;
            jbx_hi[2] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(2)-plo[2])/dz + 0.5))-1;
          }

          const Box lo_box(jbx_lo, jbx_hi);

          amrex::ParallelFor(lo_box, [bc_jlo_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_jlo_type(i,j,k,0) = type;
               bc_jlo_type(i,j,k,1) = bcv;
             });

        }

        set_bcrec_lo(lev, dir, ylo_type);

      }// end y-lo side of the domain


      { // y-hi side of the domain

        Box box_jhi = amrex::adjCellHi(domainy,1,1);

        IntVect jbx_lo(box_jhi.loVect());
        IntVect jbx_hi(box_jhi.hiVect());

        int yhi_type = init_y;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_jhi, [bc_jhi_type, init_y]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_jhi_type(i,j,k,0) = init_y;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_yhi.size(); ++lc) {

          const int bcv  = bc_yhi[lc];
          const int type = bc[bcv].type;

          yhi_type = type;

          if (lc > 0){
            jbx_lo[0] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(0)-plo[0])/dx + 0.5));
            jbx_lo[2] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(2)-plo[2])/dz + 0.5));

            jbx_hi[0] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(0)-plo[0])/dx + 0.5))-1;
            jbx_hi[2] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(2)-plo[2])/dz + 0.5))-1;
          }

          const Box hi_box(jbx_lo, jbx_hi);

          amrex::ParallelFor(hi_box, [bc_jhi_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_jhi_type(i,j,k,0) = type;
               bc_jhi_type(i,j,k,1) = bcv;
             });
        }

        set_bcrec_hi(lev, dir, yhi_type);

      } // end y-hi side of the domain
    } // end y-direction

    { // begin z-direction

      const int dir = 2;

      const int init_z = geom[lev].isPeriodic(2) ? und_ : ig_;

      Array4<int> const& bc_klo_type = bc_klo[lev]->array();
      Array4<int> const& bc_khi_type = bc_khi[lev]->array();

      Box domainz(geom[lev].Domain());
      domainz.grow(0,nghost_bc);  // Add ghost cells to x
      domainz.grow(1,nghost_bc);  // Add ghost cells to y

      { // z-lo side of the domain

        Box box_klo = amrex::adjCellLo(domainz,2,1);

        IntVect kbx_lo(box_klo.loVect());
        IntVect kbx_hi(box_klo.hiVect());

        int zlo_type = init_z;

        // Initialize y-lo domain extent.
        amrex::ParallelFor(box_klo, [bc_klo_type, init_z]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_klo_type(i,j,k,0) = init_z;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_zlo.size(); ++lc) {

          const int bcv  = bc_zlo[lc];
          const int type = bc[bcv].type;

          zlo_type = type;

          if (lc > 0){
            kbx_lo[0] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(0)-plo[0])/dx + 0.5));
            kbx_lo[1] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(1)-plo[1])/dy + 0.5));

            kbx_hi[0] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(0)-plo[0])/dx + 0.5))-1;
            kbx_hi[1] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(1)-plo[1])/dy + 0.5))-1;
          }

          const Box lo_box(kbx_lo, kbx_hi);

          amrex::ParallelFor(lo_box, [bc_klo_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_klo_type(i,j,k,0) = type;
               bc_klo_type(i,j,k,1) = bcv;
             });

        }

        set_bcrec_lo(lev, dir, zlo_type);

      } // end z-lo side of the domain


      { // z-hi side of the domain

        Box box_khi = amrex::adjCellHi(domainz,2,1);

        IntVect kbx_lo(box_khi.loVect());
        IntVect kbx_hi(box_khi.hiVect());

        int zhi_type = init_z;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_khi, [bc_khi_type, init_z]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_khi_type(i,j,k,0) = init_z;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_zhi.size(); ++lc) {

          const int bcv  = bc_zhi[lc];
          const int type = bc[bcv].type;

          zhi_type = type;

          if (lc > 0){
            kbx_lo[0] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(0)-plo[0])/dx + 0.5));
            kbx_lo[1] = static_cast<int>(amrex::Math::floor((bc[bcv].region->lo(1)-plo[1])/dy + 0.5));

            kbx_hi[0] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(0)-plo[0])/dx + 0.5))-1;
            kbx_hi[1] = static_cast<int>(amrex::Math::floor((bc[bcv].region->hi(1)-plo[1])/dy + 0.5))-1;
          }

          const Box hi_box(kbx_lo, kbx_hi);

          amrex::ParallelFor(hi_box, [bc_khi_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_khi_type(i,j,k,0) = type;
               bc_khi_type(i,j,k,1) = bcv;
             });

        }

        set_bcrec_hi(lev, dir, zhi_type);

      } // end z-hi side of the domain
    }// end z-direction


    {
      m_bcrec_velocity_d.resize(AMREX_SPACEDIM);
      m_bcrec_hydro_velocity_d.resize(AMREX_SPACEDIM);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_velocity_d.data(), m_bcrec_velocity.data(), sizeof(BCRec)*AMREX_SPACEDIM);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_hydro_velocity_d.data(), m_bcrec_hydro_velocity.data(), sizeof(BCRec)*AMREX_SPACEDIM);
    }


    if (advect_density) {
      m_bcrec_density_d.resize(1);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_density_d.data(), m_bcrec_density.data(), sizeof(BCRec));
    }

    if (advect_enthalpy) {
      m_bcrec_enthalpy_d.resize(1);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_enthalpy_d.data(), m_bcrec_enthalpy.data(), sizeof(BCRec));
    }

    if (l_ntrac > 0) {
      m_bcrec_tracer_d.resize(l_ntrac);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_tracer_d.data(), m_bcrec_tracer.data(), sizeof(BCRec)*l_ntrac);
    }

    if (l_species > 0) {
      m_bcrec_species_d.resize(l_species);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_species_d.data(), m_bcrec_species.data(), sizeof(BCRec)*l_species);
    }


    {
      m_bcrec_force_d.resize(l_force);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_force_d.data(), m_bcrec_force.data(), sizeof(BCRec)*l_force);
    }


    if ( fluid.solve ) {
      Real ltime(0.);

      set_velocity_bc_values (ltime);
      set_density_bc_values (ltime);

      if (advect_tracer ) {
        set_tracer_bc_values (ltime);
      }

      if (advect_fluid_species) {
        set_species_bc_values (ltime);
      }

      // Species
      if (advect_enthalpy ) {
        set_temperature_bc_values (ltime);
      }

    }


    m_h_bc_ep_g.resize(bc.size());
    m_h_bc_p_g.resize(bc.size());


    for(unsigned bcv(0); bcv < bc.size(); ++bcv)
    {

      if ( fluid.solve ) {
        m_h_bc_ep_g[bcv] = bc[bcv].fluid.volfrac;
        m_h_bc_p_g[bcv]  = bc[bcv].fluid.pressure;

      } else {
        m_h_bc_ep_g[bcv] = 1e50;
        m_h_bc_p_g[bcv]  = 1e50;
      }

    }

    Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_p_g.begin(), m_h_bc_p_g.end(), m_bc_p_g.begin());
    Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_ep_g.begin(), m_h_bc_ep_g.end(), m_bc_ep_g.begin());


    Gpu::synchronize();
}


void mfix::set_bcrec_lo(const int lev, const int dir, const int l_type)
{

  const int minf_ = bc_list.get_minf();
  const int pinf_ = bc_list.get_pinf();
  const int pout_ = bc_list.get_pout();
  const int nsw_  = bc_list.get_nsw();

  // Velocity BC Recs
  if (l_type == pinf_) {

    m_bcrec_velocity[0].setLo(dir, BCType::foextrap);
    m_bcrec_velocity[1].setLo(dir, BCType::foextrap);
    m_bcrec_velocity[2].setLo(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[0].setLo(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[1].setLo(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[2].setLo(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[dir].setLo(dir, BCType::ext_dir);

  } else if (l_type == pout_) {

    m_bcrec_velocity[0].setLo(dir, BCType::foextrap);
    m_bcrec_velocity[1].setLo(dir, BCType::foextrap);
    m_bcrec_velocity[2].setLo(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[0].setLo(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[1].setLo(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[2].setLo(dir, BCType::foextrap);

  } else if (l_type == minf_ || l_type == nsw_) {

    m_bcrec_velocity[0].setLo(dir, BCType::ext_dir);
    m_bcrec_velocity[1].setLo(dir, BCType::ext_dir);
    m_bcrec_velocity[2].setLo(dir, BCType::ext_dir);

    m_bcrec_hydro_velocity[0].setLo(dir, BCType::ext_dir);
    m_bcrec_hydro_velocity[1].setLo(dir, BCType::ext_dir);
    m_bcrec_hydro_velocity[2].setLo(dir, BCType::ext_dir);

  } else if (geom[lev].isPeriodic(dir)) {

    m_bcrec_velocity[0].setLo(dir, BCType::int_dir);
    m_bcrec_velocity[1].setLo(dir, BCType::int_dir);
    m_bcrec_velocity[2].setLo(dir, BCType::int_dir);

    m_bcrec_hydro_velocity[0].setLo(dir, BCType::int_dir);
    m_bcrec_hydro_velocity[1].setLo(dir, BCType::int_dir);
    m_bcrec_hydro_velocity[2].setLo(dir, BCType::int_dir);
  }

  // Scalar BC Recs
  if (l_type == pinf_ || l_type == pout_ || l_type == nsw_) {

    m_bcrec_density[0].setLo(dir, BCType::foextrap);
    m_bcrec_enthalpy[0].setLo(dir, BCType::foextrap);
    for (auto& b : m_bcrec_tracer) b.setLo(dir, BCType::foextrap);
    for (auto& b : m_bcrec_species) b.setLo(dir, BCType::foextrap);

  } else if (l_type == minf_) {

    m_bcrec_density[0].setLo(dir, BCType::ext_dir);
    m_bcrec_enthalpy[0].setLo(dir, BCType::ext_dir);
    for (auto& b : m_bcrec_tracer) b.setLo(dir, BCType::ext_dir);
    for (auto& b : m_bcrec_species) b.setLo(dir, BCType::ext_dir);

  } else if (geom[lev].isPeriodic(dir)) {

    m_bcrec_density[0].setLo(dir, BCType::int_dir);
    m_bcrec_enthalpy[0].setLo(dir, BCType::int_dir);
    for (auto& b : m_bcrec_tracer) b.setLo(dir, BCType::int_dir);
    for (auto& b : m_bcrec_species) b.setLo(dir, BCType::int_dir);
  }

  // Force BC Recs
  if (geom[lev].isPeriodic(dir)) {

    for (auto& b : m_bcrec_force) b.setLo(dir, BCType::int_dir);

  } else {

    for (auto& b : m_bcrec_force) b.setLo(dir, BCType::foextrap);
  }
}



void mfix::set_bcrec_hi(const int lev, const int dir, const int l_type)
{

  const int minf_ = bc_list.get_minf();
  const int pinf_ = bc_list.get_pinf();
  const int pout_ = bc_list.get_pout();
  const int nsw_  = bc_list.get_nsw();

  // Velocity BC Recs
  if (l_type == pinf_) {

    m_bcrec_velocity[0].setHi(dir, BCType::foextrap);
    m_bcrec_velocity[1].setHi(dir, BCType::foextrap);
    m_bcrec_velocity[2].setHi(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[0].setHi(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[1].setHi(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[2].setHi(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[dir].setHi(dir, BCType::ext_dir);

  } else if (l_type == pout_) {

    m_bcrec_velocity[0].setHi(dir, BCType::foextrap);
    m_bcrec_velocity[1].setHi(dir, BCType::foextrap);
    m_bcrec_velocity[2].setHi(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[0].setHi(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[1].setHi(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[2].setHi(dir, BCType::foextrap);

  } else if (l_type == minf_ || l_type == nsw_) {

    m_bcrec_velocity[0].setHi(dir, BCType::ext_dir);
    m_bcrec_velocity[1].setHi(dir, BCType::ext_dir);
    m_bcrec_velocity[2].setHi(dir, BCType::ext_dir);

    m_bcrec_hydro_velocity[0].setHi(dir, BCType::ext_dir);
    m_bcrec_hydro_velocity[1].setHi(dir, BCType::ext_dir);
    m_bcrec_hydro_velocity[2].setHi(dir, BCType::ext_dir);

  } else if (geom[lev].isPeriodic(dir)) {

    m_bcrec_velocity[0].setHi(dir, BCType::int_dir);
    m_bcrec_velocity[1].setHi(dir, BCType::int_dir);
    m_bcrec_velocity[2].setHi(dir, BCType::int_dir);

    m_bcrec_hydro_velocity[0].setHi(dir, BCType::int_dir);
    m_bcrec_hydro_velocity[1].setHi(dir, BCType::int_dir);
    m_bcrec_hydro_velocity[2].setHi(dir, BCType::int_dir);
  }

  // Scalar BC Recs
  if (l_type == pinf_ || l_type == pout_ || l_type == nsw_) {

    m_bcrec_density[0].setHi(dir, BCType::foextrap);
    m_bcrec_enthalpy[0].setHi(dir, BCType::foextrap);
    for (auto& b : m_bcrec_tracer) b.setHi(dir, BCType::foextrap);
    for (auto& b : m_bcrec_species) b.setHi(dir, BCType::foextrap);

  } else if (l_type == minf_) {

    m_bcrec_density[0].setHi(dir, BCType::ext_dir);
    m_bcrec_enthalpy[0].setHi(dir, BCType::ext_dir);
    for (auto& b : m_bcrec_tracer) b.setHi(dir, BCType::ext_dir);
    for (auto& b : m_bcrec_species) b.setHi(dir, BCType::ext_dir);

  } else if (geom[lev].isPeriodic(dir)) {

    m_bcrec_density[0].setHi(dir, BCType::int_dir);
    m_bcrec_enthalpy[0].setHi(dir, BCType::int_dir);
    for (auto& b : m_bcrec_tracer) b.setHi(dir, BCType::int_dir);
    for (auto& b : m_bcrec_species) b.setHi(dir, BCType::int_dir);
  }

  // Force BC Recs
  if (geom[lev].isPeriodic(dir)) {

    for (auto& b : m_bcrec_force) b.setHi(dir, BCType::int_dir);

  } else {

    for (auto& b : m_bcrec_force) b.setHi(dir, BCType::foextrap);
  }
}

void
mfix::set_velocity_bc_values (Real time_in) const
{

  m_h_bc_u_g.resize(bc.size());
  m_h_bc_v_g.resize(bc.size());
  m_h_bc_w_g.resize(bc.size());

  const int minf_ = bc_list.get_minf();

  for(unsigned bcv(0); bcv < BC::bc.size(); ++bcv) {

    if ( bc[bcv].type == minf_ ) {

      const auto& bc_vels = BC::bc[bcv].fluid.get_velocity(time_in);
      m_h_bc_u_g[bcv] = bc_vels[0];
      m_h_bc_v_g[bcv] = bc_vels[1];
      m_h_bc_w_g[bcv] = bc_vels[2];

    } else {

      m_h_bc_u_g[bcv] = 1e50;
      m_h_bc_v_g[bcv] = 1e50;
      m_h_bc_w_g[bcv] = 1e50;

    }

  }

  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_u_g.begin(), m_h_bc_u_g.end(), m_bc_u_g.begin());
  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_v_g.begin(), m_h_bc_v_g.end(), m_bc_v_g.begin());
  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_w_g.begin(), m_h_bc_w_g.end(), m_bc_w_g.begin());

  Gpu::synchronize();
}

void
mfix::set_temperature_bc_values (Real time_in) const
{
  m_h_bc_t_g.resize(bc.size());
  m_h_bc_h_g.resize(bc.size());

  const int minf_ = bc_list.get_minf();
  const int pinf_ = bc_list.get_pinf();

  // Flag to understand if fluid is a mixture
  const int fluid_is_a_mixture = fluid.is_a_mixture;

  auto& fluid_parms = *fluid.parameters;

  for(unsigned bcv(0); bcv < bc.size(); ++bcv) {
    if ( bc[bcv].type == minf_ || bc[bcv].type == pinf_ ) {
      const Real Tg = bc[bcv].fluid.get_temperature(time_in);
      m_h_bc_t_g[bcv] = Tg;
      if (!fluid_is_a_mixture) {
        m_h_bc_h_g[bcv] = fluid_parms.calc_h_g(m_h_bc_t_g[bcv]);
      } else {
        m_h_bc_h_g[bcv] = 0.0;
        for (int n(0); n < fluid.nspecies; n++) {
          const Real X_gk = bc[bcv].fluid.species[n].mass_fraction;
          m_h_bc_h_g[bcv] += X_gk*fluid_parms.calc_h_gk(Tg,n);
        }
      }

    } else {
      m_h_bc_t_g[bcv] = 1e50;
      m_h_bc_h_g[bcv] = 1e50;
    }
  }

  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_t_g.begin(), m_h_bc_t_g.end(), m_bc_t_g.begin());
  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_h_g.begin(), m_h_bc_h_g.end(), m_bc_h_g.begin());

  Gpu::synchronize();
}

void
mfix::set_density_bc_values (Real /*time_in*/) const
{

  m_h_bc_ro_g.resize(bc.size());

  const int minf_ = bc_list.get_minf();
  const int pinf_ = bc_list.get_pinf();

  // HACK -- BC density is constant given current implementation.
  // This was copied over from the mfix_set_density_bcs routine.
  const Real ro_g0 = fluid.ro_g0;

  for(unsigned bcv(0); bcv < BC::bc.size(); ++bcv) {
    if ( bc[bcv].type == minf_ || bc[bcv].type == pinf_ ) {
      m_h_bc_ro_g[bcv] = ro_g0;
    } else {
      m_h_bc_ro_g[bcv] = 1e50;
    }
  }

  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_ro_g.begin(), m_h_bc_ro_g.end(), m_bc_ro_g.begin());

  Gpu::synchronize();
}

void
mfix::set_tracer_bc_values (Real /*time_in*/) const
{
  m_h_bc_tracer.resize(bc.size());

  const int minf_ = bc_list.get_minf();
  const int pinf_ = bc_list.get_pinf();

  // HACK -- BC tracer is constant given current implementation.
  // This was copied over from the mfix_set_tracer_bcs routine.
  const Real trac0 = fluid.trac_0;

  for(unsigned bcv(0); bcv < BC::bc.size(); ++bcv) {
    if ( bc[bcv].type == minf_ || bc[bcv].type == pinf_ ) {
      m_h_bc_tracer[bcv] = trac0;
    } else {
      m_h_bc_tracer[bcv] = 1e50;
    }
  }

  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_tracer.begin(), m_h_bc_tracer.end(), m_bc_tracer.begin());

  Gpu::synchronize();
}

void
mfix::set_species_bc_values (Real /*time_in*/) const
{
  m_h_bc_X_gk.resize(fluid.nspecies, Gpu::HostVector<Real>(bc.size()));
  m_bc_X_gk.resize(fluid.nspecies, Gpu::DeviceVector<Real>(bc.size()));

  // Important! Resize the bc vector for the fluid species mass fractions
  // We have to do it here because the size has to match the number of fluid
  // species
  m_bc_X_gk_ptr.resize(fluid.nspecies, nullptr);
  m_h_bc_X_gk_ptr.resize(fluid.nspecies, nullptr);

  const int minf_ = bc_list.get_minf();
  const int pinf_ = bc_list.get_pinf();

  for(unsigned bcv(0); bcv < BC::bc.size(); ++bcv) {
    if ( bc[bcv].type == minf_ || bc[bcv].type == pinf_ ) {
      for (int n(0); n < fluid.nspecies; n++) {
        m_h_bc_X_gk[n][bcv] = bc[bcv].fluid.species[n].mass_fraction;
      }
    } else {
      for (int n(0); n < fluid.nspecies; n++) {
        m_h_bc_X_gk[n][bcv] = 1e50;
      }
    }
  }

  for (int n = 0; n < fluid.nspecies; ++n) {
    Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_X_gk[n].begin(), m_h_bc_X_gk[n].end(), m_bc_X_gk[n].begin());
    m_h_bc_X_gk_ptr[n] = m_bc_X_gk[n].dataPtr();
  }
  Gpu::copyAsync(Gpu::hostToDevice, m_h_bc_X_gk_ptr.begin(), m_h_bc_X_gk_ptr.end(), m_bc_X_gk_ptr.begin());

  Gpu::synchronize();
}
