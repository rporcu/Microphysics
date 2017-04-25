!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS3_DES                                               !
!  Author: J.Musser                                   Date: 17-Nov-14  !
!                                                                      !
!  REF: Di Renzo, A. and Di Maio F.P. "Comparison of contact-force     !
!       models for the simulation of collisions in DEM-based granular  !
!       flow codes," Chemical Engineering Science, 59(3), pg 525-541.  !
!                                                                      !
!  REF: Kharaz, A.H., Gorham, D.A., and Salman, A.D. "An experimental  !
!       study of the elastic rebound of spheres," Powder Technology,   !
!       120(3), pg 281-291.                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR3_DES(max_pip, des_pos_new, des_vel_new, omega_new)


      use constant, only: PI

      use usr, only: init_angle, init_vel_t
      use amrex_fort_module, only : c_real => amrex_real

      IMPLICIT NONE

      integer     , intent(in   ) :: max_pip
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)
      real(c_real), intent(in   ) :: omega_new(max_pip,3)

! Local variables
!---------------------------------------------------------------------//
! file name
      CHARACTER(LEN=64) :: FNAME
! file unit for heat transfer data
      INTEGER, PARAMETER :: UDF_UNIT = 2030
! Calculated tangential restitution coefficient.
      REAL(C_REAL) :: RST_COEFF(2)
! Rebound Angle (degrees)
      REAL(C_REAL) :: RBND_ANGLE(2)
! Particle loop counter
      INTEGER :: LC, NP
! YZ Velocity Mag
      REAL(C_REAL) :: VEL_YZ(2), ROT_YZ(2)

! Particle-Wall Rebound Angle (degrees)
!---------------------------------------------------------------------//
! Open the files.
      FNAME = 'POST_ALPHA.dat'
      OPEN(UNIT=UDF_UNIT,FILE=FNAME, POSITION="APPEND",STATUS='OLD')

         DO LC=1,31
! Calculate the particle-wall rebound angle.
            NP = LC+31
            VEL_YZ(1) = sqrt(des_vel_new(NP,2)**2 + des_vel_new(NP,3)**2)
            RBND_ANGLE(1) = atan(VEL_YZ(1)/des_vel_new(NP,1))*180.0/PI

! Calculate the particle-wall rebound angle.
            NP = LC
            VEL_YZ(2) = sqrt(des_vel_new(NP,2)**2 + des_vel_new(NP,3)**2)
            RBND_ANGLE(2) = atan(VEL_YZ(2)/des_vel_new(NP,1))*180.0/PI

! Write the results to a file.
            WRITE(UDF_UNIT,"(3(3x,F11.4))") INIT_ANGLE(NP), RBND_ANGLE(1:2)
         ENDDO
         CLOSE(UDF_UNIT)

! Particle-wall Tangential Restitution Coefficient
!---------------------------------------------------------------------//
! Open the files.
         FNAME = 'POST_COEFF.dat'
         OPEN(UNIT=UDF_UNIT,FILE=FNAME, POSITION="APPEND",STATUS='OLD')

         DO LC=2, 31
! Calculate the particle-wall restitution coefficient.
            NP=LC+31
            VEL_YZ(1) = sqrt(des_vel_new(NP,2)**2 + des_vel_new(NP,3)**2)
            RST_COEFF(1) =  VEL_YZ(1)/INIT_VEL_T(NP)
! Calculate the particle-wall restitution coefficient.
            NP=LC
            VEL_YZ(2) = sqrt(des_vel_new(NP,2)**2 + des_vel_new(NP,3)**2)
            RST_COEFF(2) =  VEL_YZ(2)/INIT_VEL_T(NP)
! Write the results to a file.
            WRITE(UDF_UNIT,"(3(3x,F11.4))") INIT_ANGLE(NP), RST_COEFF(1:2)
         ENDDO
         CLOSE(UDF_UNIT)

! Particle-wall Angular velocity (rad/sec)
!---------------------------------------------------------------------//
! Open the files.
         FNAME = 'POST_OMEGA.dat'
         OPEN(UNIT=UDF_UNIT,FILE=FNAME, POSITION="APPEND",STATUS='OLD')

         DO LC=1,31

! Calculate particle-particle angular velocity
            NP=LC+31
            ROT_YZ(1) = sqrt(omega_new(NP,2)**2 + omega_new(NP,3)**2)

! Calculate particle-wall angular velocity
            NP=LC
            ROT_YZ(2) = sqrt(omega_new(NP,2)**2 + omega_new(NP,3)**2)

! Write the results to a file.
            WRITE(UDF_UNIT,"(3(3x,F11.4))") INIT_ANGLE(NP), ROT_YZ(1:2)
         ENDDO
         CLOSE(UDF_UNIT)

      RETURN
      END SUBROUTINE USR3_DES
