!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_DRAG                                               !
!                                                                      !
!  Purpose: Provide a hook for user defined drag law implementation.   !
!                                                                      !
!  This routine is called from inside fluid (TFM) and particle (DES)   !
!  loops. The fluid cell index (IJK) and phase (TFM) or particle index !
!  (DES) is passed.                                                    !
!                                                                      !
!  ***************************   WARNING   **************************  !
!  *----------------------------------------------------------------*  !
!  * The dummy arguments changed in the 2015-1 MFIX Release.        *  !
!  *                                                                *  !
!  *   1) Phase index (M) is now particle index (NP) for DES. This  *  !
!  *      is reflected in the name change M --> M_NP.               *  !
!  *                                                                *  !
!  *   2) The fluid velocity was added as a dummy argument. This    *  !
!  *      provides access to the interpolated gas velocity for      *  !
!  *      coupled DES simulations.                                  *  !
!  *                                                                *  !
!  ******************************************************************  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_USR(I,J,K, M_NP, lDgA, EPg, Mug, ROg, VREL, DPM, &
         ROs, lUg, lVg, lWg)

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

! Index of fluid cell:
      INTEGER, INTENT(IN) :: I,J,K
! TFM SOLIDS --> Index of phase (M)
! DES SOLIDS --> Index of particle (NP); M = particle_phase(NP,5)
      INTEGER, INTENT(IN) :: M_NP

! drag coefficient
      real(c_real), INTENT(OUT) :: lDgA
! gas volume fraction
      real(c_real), INTENT(IN) :: EPg
! gas laminar viscosity
      real(c_real), INTENT(IN) :: Mug
! gas density
      real(c_real), INTENT(IN) :: ROg
! Magnitude of gas-solids relative velocity
      real(c_real), INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      real(c_real), INTENT(IN) :: DPM
! particle density of solids phase M
      real(c_real), INTENT(IN) :: ROs
! fluid velocity components:
! o TFM: Averaged from faces to cell center
! o DES: Interpolated to the particle's position
      real(c_real), INTENT(IN) :: lUg, lVg, lWg


! The following error message is used to make sure that if a user
! defined drag law is invoked, that this routine has been modified.


!- REMOVE THE FOLLOWING ---------------------------------------------->>

      lDgA = 0.0

      CALL INIT_ERR_MSG('USR_DRAG')
      WRITE(ERR_MSG,9999)
      CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 9999 FORMAT('ERROR 9999: The user-defined drag routine was invoked ', &
         'but this',/'generic error message exits. Either choose a ',  &
         'different drag law',/'or correct mfix/model/usr_drag.f')

!- END REMOVE --------------------------------------------------------<<

      RETURN
      END SUBROUTINE DRAG_USR
