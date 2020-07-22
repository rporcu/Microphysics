#include <mfix.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

#include <mfix_bc_parms.H>
#include <mfix_fluid_parms.H>

#include <math.h>

using namespace BC;

void
mfix::mfix_set_bc_type (int lev)
{
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    const GpuArray<Real, 3> plo = geom[lev].ProbLoArray();

    // Extract the lower and upper boundaries of Box Domain
    const int und_  = bc_list.get_undefined();
    const int ig_   = bc_list.get_ig();
    const int minf_ = bc_list.get_minf();
    const int pinf_ = bc_list.get_pinf();

    {
      Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
      Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();

      const int init_x = geom[lev].isPeriodic(0) ? und_ : ig_;

      Box domainx(geom[lev].Domain());
      domainx.grow(1,nghost);  // Add ghost cells to y
      domainx.grow(2,nghost);  // Add ghost cells to z

      { // x-lo side of the domain

        Box box_ilo = amrex::adjCellLo(domainx,0,1);

        IntVect ibx_lo(box_ilo.loVect());
        IntVect ibx_hi(box_ilo.hiVect());

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_ilo, [bc_ilo_type, init_x]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_ilo_type(i,j,k,0) = init_x;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_xlo.size(); ++lc) {

          const int bcv  = bc_xlo[lc];
          const int type = bc[bcv].type;

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
      }


      { // x-hi side of the domain

        Box box_ihi = amrex::adjCellHi(domainx,0,1);

        IntVect ibx_lo(box_ihi.loVect());
        IntVect ibx_hi(box_ihi.hiVect());

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_ihi, [bc_ihi_type, init_x]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_ihi_type(i,j,k,0) = init_x;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_xhi.size(); ++lc) {

          const int bcv  = bc_xhi[lc];
          const int type = bc[bcv].type;

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
      }
    }



    {

      const int init_y = geom[lev].isPeriodic(1) ? und_ : ig_;

      Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
      Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();

      Box domainy(geom[lev].Domain());
      domainy.grow(0,nghost);  // Add ghost cells to x
      domainy.grow(2,nghost);  // Add ghost cells to z

      { // y-lo side of the domain

        Box box_jlo = amrex::adjCellLo(domainy,1,1);

        IntVect jbx_lo(box_jlo.loVect());
        IntVect jbx_hi(box_jlo.hiVect());

        // Initialize y-lo domain extent.
        amrex::ParallelFor(box_jlo, [bc_jlo_type, init_y]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_jlo_type(i,j,k,0) = init_y;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_ylo.size(); ++lc) {

          const int bcv  = bc_ylo[lc];
          const int type = bc[bcv].type;

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
      }


      { // y-hi side of the domain

        Box box_jhi = amrex::adjCellHi(domainy,1,1);

        IntVect jbx_lo(box_jhi.loVect());
        IntVect jbx_hi(box_jhi.hiVect());

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_jhi, [bc_jhi_type, init_y]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_jhi_type(i,j,k,0) = init_y;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_yhi.size(); ++lc) {

          const int bcv  = bc_yhi[lc];
          const int type = bc[bcv].type;

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
      }
    }




    {
      const int init_z = geom[lev].isPeriodic(2) ? und_ : ig_;

      Array4<int> const& bc_klo_type = bc_klo[lev]->array();
      Array4<int> const& bc_khi_type = bc_khi[lev]->array();

      Box domainz(geom[lev].Domain());
      domainz.grow(0,nghost);  // Add ghost cells to x
      domainz.grow(1,nghost);  // Add ghost cells to y

      { // z-lo side of the domain

        Box box_klo = amrex::adjCellLo(domainz,2,1);

        IntVect kbx_lo(box_klo.loVect());
        IntVect kbx_hi(box_klo.hiVect());

        // Initialize y-lo domain extent.
        amrex::ParallelFor(box_klo, [bc_klo_type, init_z]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_klo_type(i,j,k,0) = init_z;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_zlo.size(); ++lc) {

          const int bcv  = bc_zlo[lc];
          const int type = bc[bcv].type;

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
      }


      { // z-hi side of the domain

        Box box_khi = amrex::adjCellHi(domainz,2,1);

        IntVect kbx_lo(box_khi.loVect());
        IntVect kbx_hi(box_khi.hiVect());

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_khi, [bc_khi_type, init_z]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_khi_type(i,j,k,0) = init_z;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_zhi.size(); ++lc) {

          const int bcv  = bc_zhi[lc];
          const int type = bc[bcv].type;

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
      }
    }

    for(unsigned bcv(0); bcv < bc.size(); ++bcv)
    {

      if ( FLUID::solve and bc[bcv].type == minf_ ) {

        m_bc_u_g[bcv] = bc[bcv].fluid.velocity[0];
        m_bc_v_g[bcv] = bc[bcv].fluid.velocity[1];
        m_bc_w_g[bcv] = bc[bcv].fluid.velocity[2];

      } else {

        m_bc_u_g[bcv] = 1e50;
        m_bc_v_g[bcv] = 1e50;
        m_bc_w_g[bcv] = 1e50;

      }

      if ( FLUID::solve ) {
        m_bc_ep_g[bcv] = bc[bcv].fluid.volfrac;
        m_bc_p_g[bcv]  = bc[bcv].fluid.pressure;

      } else {
        m_bc_ep_g[bcv] = 1e50;
        m_bc_p_g[bcv]  = 1e50;
      }

      // Fluid temperature
      if ( FLUID::solve and advect_enthalpy ) {
        if ( bc[bcv].type == minf_ or bc[bcv].type == pinf_ ) {
          m_bc_t_g[bcv]  = bc[bcv].fluid.temperature;

        } else {
          m_bc_t_g[bcv]  = 1e50;
        }
      }

      // Fluid species mass fractions
      if ( FLUID::solve and advect_fluid_species) {
        if ( bc[bcv].type == minf_ or bc[bcv].type == pinf_ ) {
          for (int n(0); n < FLUID::nspecies_g; n++) {
            m_bc_X_g[n][bcv] = bc[bcv].fluid.species.mass_fractions[n];
          }
        }
        else {
          for (int n(0); n < FLUID::nspecies_g; n++)
            m_bc_X_g[n][bcv] = 1e50;
        }
      }

    }
}
