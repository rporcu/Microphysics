MODULE CHECK_SOLIDS_COMMON_ALL_MODULE
   use param1, only: is_undefined, is_defined
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_ALL                                 !
!  Purpose: Check the solid phase input that is common to all solids   !
!  phase models.                                                       !
!                                                                      !
!    ****** DO NOT PUT MODEL SPECIFIC CHECKS IN THIS ROUTINE ******    !
!                                                                      !
!  Use the companion routines for checks specific to a particular      !
!  solids phase model:                                                 !
!                                                                      !
!    > CHECK_SOLIDS_DEM       :: DEM solids phase model                !
!                                                                      !
!  Author: J.Musser                                  Date: 03-FEB-14   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_ALL


! Global Variables:
!---------------------------------------------------------------------//
! Number of continuum solids phases.
      use constant, only: MMAX
! User specified: Initial solids diameter.
      use constant, only: D_P0


! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases.
      use param, only: DIM_M
! Parameter constants
      use param1, only: ZERO


! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg

      implicit none


! Local Variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER :: M
! Total number of all solids
      INTEGER :: MMAX_L

!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_ALL")

! Set the number of solids phases to be checked.
      MMAX_L = MMAX

! Check D_p0
      DO M = 1, MMAX_L
         IF(IS_UNDEFINED(D_P0(M))) THEN
            WRITE(ERR_MSG, 1000) trim(iVar('D_p0',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(D_P0(M) <= ZERO)THEN
            WRITE(ERR_MSG, 1001) trim(iVar('D_p0',M)), iVal(D_P0(M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

      DO M = MMAX_L+1, DIM_M
         IF(IS_DEFINED(D_P0(M)))THEN
            WRITE(ERR_MSG,1002) trim(iVar('D_p0',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

! Check solids drag model selection.
      CALL CHECK_SOLIDS_DRAG

! Check the solids density input parameters.
      CALL CHECK_SOLIDS_DENSITY(MMAX_L)

! Finalize the error messges
      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal input: ',A,' specified out of range.', &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_COMMON_ALL


!----------------------------------------------------------------------!
! Subroutine: CHECK_SOLIDS_DRAG                                        !
! Purpose: Check solids species input.                                 !
!                                                                      !
! Author: J. Musser                                  Date: 07-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_DRAG

! Global Variables:
!---------------------------------------------------------------------//
! User specifed drag type, as string and enum
      use run, only: DRAG_TYPE
      use run, only: DRAG_TYPE_ENUM
! Possible DRAG_TYPE_ENUM values:
      use run, only: SYAM_OBRIEN
      use run, only: GIDASPOW
      use run, only: GIDASPOW_PCF
      use run, only: GIDASPOW_BLEND
      use run, only: GIDASPOW_BLEND_PCF
      use run, only: WEN_YU
      use run, only: WEN_YU_PCF
      use run, only: KOCH_HILL
      use run, only: KOCH_HILL_PCF
      use run, only: BVK
      use run, only: HYS
      use run, only: USER_DRAG

! Global Parameters:
!---------------------------------------------------------------------//

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg

      implicit none


! Local Variables:
!---------------------------------------------------------------------//
!     NONE

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DRAG")

      SELECT CASE(trim(adjustl(DRAG_TYPE)))

      CASE ('SYAM_OBRIEN'); DRAG_TYPE_ENUM = SYAM_OBRIEN
      CASE ('GIDASPOW'); DRAG_TYPE_ENUM = GIDASPOW
      CASE ('GIDASPOW_PCF'); DRAG_TYPE_ENUM = GIDASPOW_PCF
      CASE ('GIDASPOW_BLEND'); DRAG_TYPE_ENUM = GIDASPOW_BLEND
      CASE ('GIDASPOW_BLEND_PCF'); DRAG_TYPE_ENUM = GIDASPOW_BLEND_PCF
      CASE ('WEN_YU'); DRAG_TYPE_ENUM = WEN_YU
      CASE ('WEN_YU_PCF'); DRAG_TYPE_ENUM = WEN_YU_PCF
      CASE ('KOCH_HILL'); DRAG_TYPE_ENUM = KOCH_HILL
      CASE ('KOCH_HILL_PCF'); DRAG_TYPE_ENUM = KOCH_HILL_PCF
      CASE ('BVK'); DRAG_TYPE_ENUM = BVK
      CASE ('HYS'); DRAG_TYPE_ENUM = HYS
      CASE ('USER_DRAG','USR_DRAG'); DRAG_TYPE_ENUM = USER_DRAG

      CASE DEFAULT
         WRITE(ERR_MSG,1001)'DRAG_TYPE', trim(adjustl(DRAG_TYPE))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      END SELECT

      CALL FINL_ERR_MSG

      RETURN

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_DRAG



!----------------------------------------------------------------------!
!  Subroutine: CHECK_SOLIDS_DENSITY                                    !
!  Purpose: check the solid phase density input                        !
!                                                                      !
!  Author: J.Musser                                  Date: 03-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_DENSITY(MMAX_LL)


! Global Variables:
!---------------------------------------------------------------------//
! User specified: constant solids density
      use constant, only: RO_s0

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases.
      use param, only: DIM_M

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg


      implicit none


! Arguments:
!---------------------------------------------------------------------//
! Total number of solids phases
      INTEGER, intent(in) :: MMAX_LL

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER :: M


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DENSITY")

! Check each solids phase.
      DO M = 1, MMAX_LL

! Verify that one -and only one- solids density model is in use.
         IF(IS_UNDEFINED(RO_S0(M))) THEN
            WRITE(ERR_MSG, 1100) M
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error 1101: No solids density information for phase ',  &
         I2,'.',/'Please correct the mfix.dat file.')

         ENDIF

      ENDDO

! Check for input overflow.
      DO M = MMAX_LL+1, DIM_M
         IF(IS_DEFINED(RO_S0(M))) THEN
            WRITE(ERR_MSG,1002) trim(iVar('RO_s0',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

! Finalize the error messges
      CALL FINL_ERR_MSG

      RETURN

 1002 FORMAT('Error 1002: Illegal input: ',A,' specified out of range.',&
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_DENSITY
END MODULE CHECK_SOLIDS_COMMON_ALL_MODULE
