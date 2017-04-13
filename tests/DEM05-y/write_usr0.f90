!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_USR0                                             !
!  Purpose: Write initial part of user-defined output                  !
!                                                                      !
!  Author:                                            Date: dd-mmm-yy  !
!  Reviewer:                                          Date: dd-mmm-yy  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine write_usr0() &
        bind(C, name="write_usr0")

      implicit none

      CALL WRITE_DAT_HEADER('POST_ALPHA.dat','ALPHA')
      CALL WRITE_DAT_HEADER('POST_OMEGA.dat','OMEGA')
      CALL WRITE_DAT_HEADER('POST_COEFF.dat','COEFF')

      contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_DAT_HEADER(FNAME, VAR)

      use run, only: DESCRIPTION
      use discretelement, only: des_coll_model
      use constant, only: MMAX

      IMPLICIT NONE

      INTEGER :: M, N

      CHARACTER(len=*) :: FNAME
      CHARACTER(len=*) :: VAR

! logical used for testing is the data file already exists
      LOGICAL :: EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: fUNIT = 2030

      INQUIRE(FILE=FNAME,EXIST=EXISTS)
      IF (.NOT.EXISTS) THEN
         OPEN(UNIT=fUNIT,FILE=FNAME,STATUS='NEW')
         WRITE(fUNIT, 1000) trim(DESCRIPTION)

         IF(VAR=='ALPHA')THEN
            CALL WRITE_EXP_ALPHA(fUNIT)
         ELSEIF(VAR=='OMEGA') THEN
            CALL WRITE_EXP_OMEGA(fUNIT)
         ELSEIF(VAR=='COEFF') THEN
            CALL WRITE_EXP_COEFF(fUNIT)
         ENDIF

      ELSE
         OPEN(UNIT=fUNIT,FILE=FNAME,POSITION="APPEND",STATUS='OLD')
      ENDIF

      select case (trim(des_coll_model))
      case('LSD');      write(funit,1450) 'LINEAR SPRING-DASHPOT'
      case('HERTZIAN'); write(funit,1450) 'HERTZIAN SPRING-DASHPOT'
      end select
1450  FORMAT(3/4X,'Collision model: ',A)


      WRITE(fUNIT, 1200)


 1000 FORMAT(25x,A)


 1200 FORMAT(/6X,'Impact Angle',2X,'Part-Part',5x,'Part-Wall')

      CLOSE(fUNIT)
      RETURN
      END SUBROUTINE WRITE_DAT_HEADER

!----------------------------------------------------------------------!
!                                                                      !
!  Experimentally measured rebound angle. Data taken from Figure 7(b)  !
!  of Di Renzo and Di Maio (2004).                                     !
!                                                                      !
!  REF: Di Renzo, A. and Di Maio F.P. "Comparison of contact-force     !
!       models for the simulation of collisions in DEM-based granular  !
!       flow codes," Chemical Engineering Science, 59(3), pg 525-541.  !
!                                                                      !
!  REF: Kharaz, A.H., Gorham, D.A., and Salman, A.D. "An experimental  !
!       study of the elastic rebound of spheres," Powder Technology,   !
!       120(3), pg 281-291.                                            !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_EXP_ALPHA(lUNIT)

      INTEGER, INTENT(IN) :: lUNIT

! Number of data points.
      INTEGER, PARAMETER :: ANGLES = 9
! Impact angle.
      DOUBLE PRECISION, PARAMETER :: IMPT_ANGLE(ANGLES) =              &
         (/  0.2397,  2.2945,  4.2922,  9.9429, 20.5594,               &
            30.3767, 39.6233, 49.7831, 60.0571 /)
! Experimentally measured rebound angles.
      DOUBLE PRECISION, PARAMETER :: EXP_ANGLE(ANGLES) =               &
         (/  0.1194,  1.7630,  2.9577,  6.6917, 12.3645,               &
            21.7809, 32.5455, 45.1051, 57.3652 /)

! Loop counter
      INTEGER :: LC

      WRITE(lUNIT,1000)
 1000 FORMAT(/4x,'Impact Angle',4x,'Experiment')

      DO LC=1, ANGLES
! Write the results to a file.
         WRITE(lUNIT,1100) IMPT_ANGLE(LC), EXP_ANGLE(LC)
      ENDDO

 1100 FORMAT(2(3x,F11.4))
      END SUBROUTINE WRITE_EXP_ALPHA


!----------------------------------------------------------------------!
!                                                                      !
!  Experimentally measured post-collision angular velocity. Data taken !
!  from Figure 8 of Di Renzo and Di Maio (2004).                       !
!                                                                      !
!  REF: Di Renzo, A. and Di Maio F.P. "Comparison of contact-force     !
!       models for the simulation of collisions in DEM-based granular  !
!       flow codes," Chemical Engineering Science, 59(3), pg 525-541.  !
!                                                                      !
!  REF: Kharaz, A.H., Gorham, D.A., and Salman, A.D. "An experimental  !
!       study of the elastic rebound of spheres," Powder Technology,   !
!       120(3), pg 281-291.                                            !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_EXP_OMEGA(lUNIT)

      INTEGER, INTENT(IN) :: lUNIT

! Number of data points.
      INTEGER, PARAMETER :: OMEGAS = 7
! Impact Angle.
      DOUBLE PRECISION, PARAMETER :: IMPT_ANGLE(OMEGAS) =              &
         (/  6.0025, 10.9840, 21.1931, 31.0947,                        &
            39.0283, 50.2214, 60.2460 /)
! Experimentally measured post-collision angluar velocities
      DOUBLE PRECISION, PARAMETER :: EXP_OMEGA(OMEGAS) =               &
         (/ 135.5703, 264.5493, 577.0096, 612.6481,                    &
            595.4843, 452.9444, 381.3185 /)

! Loop counter
      INTEGER :: LC

      WRITE(lUNIT,1000)
 1000 FORMAT(/4x,'Impact Angle',4x,'Experiment')

      DO LC=1, OMEGAS
! Write the results to a file.
         WRITE(lUNIT,1100) IMPT_ANGLE(LC), EXP_OMEGA(LC)
      ENDDO

 1100 FORMAT(2(3x,F11.4))

      END SUBROUTINE WRITE_EXP_OMEGA

!----------------------------------------------------------------------!
!                                                                      !
!  Experimentally measured tangetial coefficient of restitution as a   !
!  function of impact angle. Data taken from Figure 6 of Di Renzo and  !
!  Di Maio (2004).                                                     !
!                                                                      !
!  REF: Di Renzo, A. and Di Maio F.P. "Comparison of contact-force     !
!       models for the simulation of collisions in DEM-based granular  !
!       flow codes," Chemical Engineering Science, 59(3), pg 525-541.  !
!                                                                      !
!  REF: Kharaz, A.H., Gorham, D.A., and Salman, A.D. "An experimental  !
!       study of the elastic rebound of spheres," Powder Technology,   !
!       120(3), pg 281-291.                                            !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_EXP_COEFF(lUNIT)

      INTEGER, INTENT(IN) :: lUNIT

! Number of data points
      INTEGER, PARAMETER :: COEFFS = 8
! Impact Angle.
      DOUBLE PRECISION, PARAMETER :: IMPT_ANGLE(COEFFS) =              &
         (/  2.0275,  3.9282,  9.6304, 20.0845,                        &
            30.2218, 39.3031, 49.5459, 59.8416 /)
! Experimentally measured post-collision
      DOUBLE PRECISION, PARAMETER :: EXP_COEFF(COEFFS) =               &
         (/ 0.7938, 0.7450, 0.6728, 0.5964,                            &
            0.6883, 0.7719, 0.8469, 0.9019 /)

      INTEGER :: LC

      WRITE(lUNIT,1000)
 1000 FORMAT(/4x,'Impact Angle',4x,'Experiment')

      DO LC=1, COEFFS
! Write the results to a file.
         WRITE(lUNIT,1100) IMPT_ANGLE(LC), EXP_COEFF(LC)
      ENDDO

 1100 FORMAT(2(3x,F11.4))

      END SUBROUTINE WRITE_EXP_COEFF

      END SUBROUTINE WRITE_USR0
