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
      subroutine set_domain(slo,shi,lo,hi,flag,dx,dy,dz) &
          bind(C, name="set_domain")

      ! Cyclic domain flags.
      use geometry, only: CYCLIC
      use geometry, only: CYCLIC_X, CYCLIC_X_PD, CYCLIC_X_MF
      use geometry, only: CYCLIC_Y, CYCLIC_Y_PD, CYCLIC_Y_MF
      use geometry, only: CYCLIC_Z, CYCLIC_Z_PD, CYCLIC_Z_MF

      ! Flag for specificed constant mass flux.
      use bc, only: Flux_g

      use set_icbc_flags_module, only: set_icbc_flag
      use get_bc_area_module, only: get_bc_area
      use set_bc_flow_module, only: set_bc_flow
      use set_flags_module, only: set_flags
!     use set_flags_module, only: set_flags1

      use param1, only: is_defined

      use geometry, only: flag_mod

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(inout) :: flag(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),1)
      real(c_real)  , intent(in   ) :: dx,dy,dz

      ! This used to be in set_geometry

      !  Determine the cyclic direction with a specified mass flux
      CYCLIC_X_MF = (IS_DEFINED(FLUX_G) .AND. CYCLIC_X_PD)
      CYCLIC_Y_MF = (IS_DEFINED(FLUX_G) .AND. CYCLIC_Y_PD)
      CYCLIC_Z_MF = (IS_DEFINED(FLUX_G) .AND. CYCLIC_Z_PD)

      ! Force the cyclic flag if cyclic with pressure drop.
      IF (CYCLIC_X_PD) CYCLIC_X = .TRUE.
      IF (CYCLIC_Y_PD) CYCLIC_Y = .TRUE.
      IF (CYCLIC_Z_PD) CYCLIC_Z = .TRUE.
      CYCLIC = CYCLIC_X .OR. CYCLIC_Y .OR. CYCLIC_Z

      ! End of what used to be in set_geometry

      call check_initial_conditions(dx,dy,dz)
      call check_boundary_conditions(dx,dy,dz)
      call check_point_sources(dx,dy,dz)

      ! This call needs to occur before any of the IC/BC checks.
      call set_icbc_flag(slo,shi,flag)

      ! Compute area of boundary surfaces.
      call get_bc_area(dx,dy,dz)

      ! Convert (mass, volume) flows to velocities.
      call set_bc_flow

      ! Set the flags for identifying computational cells
      call set_flags(slo,shi,lo,hi,flag)
      flag_mod = flag

      ! Set the flags for wall surfaces impermeable and identify flow
      ! boundaries using FLAG_E, FLAG_N, and FLAG_T
!     call set_flags1(slo,shi,lo,hi,flag)
      flag_mod = flag

      end subroutine set_domain
end module set_domain_module
