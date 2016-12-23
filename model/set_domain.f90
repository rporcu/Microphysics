MODULE set_domain_module

   use check_boundary_conditions_module, only: check_boundary_conditions
   use check_initial_conditions_module, only: check_initial_conditions
   use check_point_sources_module, only: check_point_sources

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
      SUBROUTINE SET_DOMAIN(flag)

      use compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

      ! Domain decomposition and dimensions
      use geometry, only: DX, oDX
      use geometry, only: DY, oDZ
      use geometry, only: DZ, oDY

      ! Cyclic domain flags.
      use geometry, only: CYCLIC
      use geometry, only: CYCLIC_X, CYCLIC_X_PD, CYCLIC_X_MF
      use geometry, only: CYCLIC_Y, CYCLIC_Y_PD, CYCLIC_Y_MF
      use geometry, only: CYCLIC_Z, CYCLIC_Z_PD, CYCLIC_Z_MF

      use set_geometry_des_module, only: set_geometry_des

      ! Flag for specificed constant mass flux.
      use bc, only: Flux_g

      use set_icbc_flags_module, only: set_icbc_flag
      use get_bc_area_module, only: get_bc_area
      use set_bc_flow_module, only: set_bc_flow
      use set_flags_module, only: set_flags


      use param1, only: one, is_defined

      use run, only: dem_solids

      implicit none

      integer, intent(inout) :: flag(istart3:iend3,jstart3:jend3,kstart3:kend3,4)

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

      ODX = ONE/DX
      ODY = ONE/DY
      ODZ = ONE/DZ

      ! End of what used to be in set_geometry

      IF(DEM_SOLIDS) CALL SET_GEOMETRY_DES

      CALL CHECK_INITIAL_CONDITIONS
      CALL CHECK_BOUNDARY_CONDITIONS
      CALL CHECK_POINT_SOURCES

! This call needs to occur before any of the IC/BC checks.
      CALL SET_ICBC_FLAG(flag)

! Compute area of boundary surfaces.
      CALL GET_BC_AREA

! Convert (mass, volume) flows to velocities.
      CALL SET_BC_FLOW

! Set the flags for identifying computational cells
      CALL SET_FLAGS(flag)


      end subroutine set_domain
end module set_domain_module
