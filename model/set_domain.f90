module set_domain_module

   use check_boundary_conditions_module, only: check_boundary_conditions
   use check_initial_conditions_module, only: check_initial_conditions
   use check_point_sources_module, only: check_point_sources

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: SET_DOMAIN                                              !
!  Purpose: read and verify input data, open files                     !
!                                                                      !
!  Author: P. Nicoletti                               Date: 04-DEC-91  !
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine set_domain(dx,dy,dz) &
          bind(C, name="set_domain")

      ! Cyclic domain flags.
      use geometry, only: cyclic
      use geometry, only: cyclic_x, cyclic_x_pd, cyclic_x_mf
      use geometry, only: cyclic_y, cyclic_y_pd, cyclic_y_mf
      use geometry, only: cyclic_z, cyclic_z_pd, cyclic_z_mf

      ! Flag for specificed constant mass flux.
      use bc, only: Flux_g

      use get_bc_area_module, only: get_bc_area
      use set_bc_flow_module, only: set_bc_flow

      use param1, only: is_defined

      implicit none

      real(c_real)  , intent(in   ) :: dx,dy,dz

      ! This used to be in set_geometry

      !  Determine the cyclic direction with a specified mass flux
      cyclic_x_mf = (IS_DEFINED(FLUX_G) .AND. cyclic_x_pd)
      cyclic_y_mf = (IS_DEFINED(FLUX_G) .AND. cyclic_y_pd)
      cyclic_z_mf = (IS_DEFINED(FLUX_G) .AND. cyclic_z_pd)

      ! Force the cyclic flag if cyclic with pressure drop.
      IF (cyclic_x_pd) cyclic_x = .TRUE.
      IF (cyclic_y_pd) cyclic_y = .TRUE.
      IF (cyclic_z_pd) cyclic_z = .TRUE.
      cyclic = cyclic_x .OR. cyclic_y .OR. cyclic_z

      ! End of what used to be in set_geometry

      call check_initial_conditions(dx,dy,dz)
      call check_boundary_conditions(dx,dy,dz)
      call check_point_sources(dx,dy,dz)

      ! Compute area of boundary surfaces.
      call get_bc_area(dx,dy,dz)

      ! Convert (mass, volume) flows to velocities.
      call set_bc_flow

      end subroutine set_domain
end module set_domain_module
