!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!                                                                      C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR0

      use usr
      use discretelement, only: PARTICLES
      use discretelement, only: DES_VEL_NEW
      use constant, only: PI
      use compar, only: myPE, PE_IO

      IMPLICIT NONE

      INTEGER :: NP
      double precision :: lTMP(62,3)

      IF(PARTICLES /= 93) THEN
         write(*,"(3x, 'invalid setup for test case')")
         call mfix_exit(0)
      ENDIF

      INIT_VEL_T = 0.0d0
      INIT_ANGLE = 0.0d0

! Collect initial translational velocity to IO processor.
      CALL COLLECT_DEM05_DATA(DES_VEL_NEW(:,1), lTMP(:,1))
      CALL COLLECT_DEM05_DATA(DES_VEL_NEW(:,2), lTMP(:,2))
      CALL COLLECT_DEM05_DATA(DES_VEL_NEW(:,3), lTMP(:,3))

! Store the collision angle and initial tangential velocity
      IF(myPE == PE_IO) THEN
         DO NP=1, 62
            INIT_VEL_T(NP) = sqrt(lTMP(NP,1)**2 + lTMP(NP,3)**2)
            INIT_ANGLE(NP) = abs(atan(INIT_VEL_T(NP)/lTMP(NP,2)))*180.0/PI
         ENDDO
      ENDIF

      return
      END SUBROUTINE USR0
