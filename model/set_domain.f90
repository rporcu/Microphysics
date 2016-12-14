MODULE SET_DOMAIN_MODULE
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

      use set_geo_mod, only: set_geometry
      use set_geometry_des_module, only: set_geometry_des
      use set_icbc_flags_module, only: set_icbc_flag
      use get_bc_area_module, only: get_bc_area
      use set_bc_flow_module, only: set_bc_flow
      use set_flags_module, only: set_flags
      use desgrid        , only: desgrid_init
      use mpi_init_des   , only: desmpi_init
      use stl_preproc_des, only: des_stl_preprocessing

      use run, only: dem_solids

      implicit none

      INTEGER, allocatable, intent(INOUT) :: FLAG(:,:,:,:)

! Set grid spacing variables.
      CALL SET_GEOMETRY(flag)
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

      IF(DEM_SOLIDS) THEN
         CALL DESGRID_INIT
         CALL DESMPI_INIT
         CALL DES_STL_PREPROCESSING
      ENDIF

      END SUBROUTINE SET_DOMAIN
END MODULE SET_DOMAIN_MODULE
