!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS                                           !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Driver routine to call checks for WALL BCs.                 !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS(M_TOT, SKIP, BCV)


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Identifies solids model (TFM,DEM,PIC)
      use run, only: SOLIDS_MODEL

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases
      use param, only: DIM_M

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
! Index of BC being checked.
      INTEGER, INTENT(in) :: BCV
! Total number of solids phases.
      INTEGER, INTENT(in) :: M_TOT
! Flag. Solids not present at this BC (used for flow BCs).
      LOGICAL, INTENT(in) :: SKIP(DIM_M)

! Local Variables:
!---------------------------------------------------------------------//
! Local total number of solids phases
      INTEGER :: MTOT_L
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
! Flag: Solve K-th direction (3D)
      use geometry, only: DO_K

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants.
      use param1, only: ZERO, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: BCV

!......................................................................!


! Initialize the error manger.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS_GAS")

! The wall velocities are not needed for no-slip or free-slip
      IF(BC_TYPE(BCV) == 'PAR_SLIP_WALL') THEN
         IF(BC_UW_G(BCV) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_Uw_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(BC_VW_G(BCV) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_Vw_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(BC_WW_G(BCV) == UNDEFINED) THEN
            IF(DO_K)THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_Ww_g',BCV))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSE
               BC_WW_G(BCV) = ZERO
            ENDIF
         ENDIF
      ENDIF


! Clear the error manager.
      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_WALLS_GAS
