MODULE SET_BC_DEM_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM                                              !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose: Check the data provided for the des mass inflow boundary   !
!  condition and flag errors if the data is improper.  This module is  !
!  also used to convert the proveded information into the format       !
!  necessary for the dependent subrountines to function properly.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_DEM(flag)

      USE des_bc, only: dem_bcmo, dem_mio, dem_bcmi
      USE discretelement, only: particles
      USE error_manager, only: err_msg, init_err_msg, flush_err_msg, finl_err_msg
      USE layout_mi_dem_module, only: layout_mi_dem
      USE param1, only: undefined_i
      USE set_bc_dem_mi_module, only: set_bc_dem_mi
      USE set_bc_dem_mo_module, only: set_bc_dem_mo

      IMPLICIT NONE

      INTEGER, DIMENSION(:,:,:,:), INTENT(IN) :: FLAG

      CALL INIT_ERR_MSG("SET_BC_DEM")

! The variable PARTICLES should already be set by this point if using
! gener_part_config option
      IF(PARTICLES == UNDEFINED_I) THEN
         PARTICLES = 0
      ENDIF

! If the system is started without any particles and an inlet is not
! specified, the run is likely aborted.
      IF(PARTICLES == 0 .AND. DEM_BCMI == 0) THEN
         WRITE(ERR_MSG, 1202)
         CALL FLUSH_ERR_MSG
      ENDIF

 1202 FORMAT('WARNING 1202: The system is initiated with no particles',&
         ' and no',/'solids inlet was detected.')

      IF(DEM_BCMI > 0) CALL SET_BC_DEM_MI(flag)
      IF(DEM_BCMO > 0) CALL SET_BC_DEM_MO

! Set the flag that one or more DEM MI/MO exists.
      DEM_MIO = (DEM_BCMI /= 0 .OR. DEM_BCMO /= 0)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_BC_DEM
END MODULE SET_BC_DEM_MODULE
