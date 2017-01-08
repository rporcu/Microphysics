MODULE CHECK_BC_WALLS_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS                                           !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Driver routine to call checks for WALL BCs.                 !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS(BCV)

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases
      use param, only: DIM_M

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
! Index of BC being checked.
      INTEGER, INTENT(in) :: BCV
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS")

! Input checks for gas phase.
      CALL CHECK_BC_WALLS_GAS(BCV)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_BC_WALLS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS_GAS                                       !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Check user-input for gas phase WALL BC parameters.          !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS_GAS(BCV)

! Global Variables:
!---------------------------------------------------------------------//
! User-input: type of BC
      use bc, only: BC_TYPE
! User-Input: gas velocity at wall BCs.
      use bc, only: BC_UW_G, BC_VW_G, BC_WW_G

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants.
      use param1, only: IS_UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: BCV

!......................................................................!


! Initialize the error manger.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS_GAS")

! The wall velocities are not needed for no-slip or free-slip
      IF(BC_TYPE(BCV) == 'PAR_SLIP_WALL') THEN
         IF(IS_UNDEFINED(BC_UW_G(BCV))) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_Uw_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IS_UNDEFINED(BC_VW_G(BCV))) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_Vw_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IS_UNDEFINED(BC_WW_G(BCV))) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_Ww_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF


! Clear the error manager.
      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_WALLS_GAS
END MODULE CHECK_BC_WALLS_MODULE
