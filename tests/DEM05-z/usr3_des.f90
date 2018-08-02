!$vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
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
subroutine usr3_des ( np, particles )

   use constant,          only: PI
   use usr,               only: init_angle, init_vel_t
   use amrex_fort_module, only: c_real => amrex_real
   use particle_mod,      only: particle_t
   use run,               only: des_tstart, des_dt, tstop

   implicit none

   integer,          intent(in   ) :: np
   type(particle_t), intent(in   ) :: particles(np)

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
   INTEGER :: LC, P
   ! YZ Velocity Mag
   REAL(C_REAL) :: VEL_YZ(2), ROT_YZ(2)

   if ( des_tstart + 1.5*des_dt .ge. tstop ) then
       ! Particle-Wall Rebound Angle (degrees)
       !---------------------------------------------------------------------//
       ! Open the files.
       FNAME = 'POST_ALPHA.dat'
       open(UNIT=UDF_UNIT,FILE=FNAME, POSITION="APPEND",STATUS='OLD')

       do LC = 1, 31
           ! Calculate the particle-wall rebound angle.
           P = LC+31
           VEL_YZ(1) = sqrt(particles(p) % vel(1)**2 + particles(p) % vel(2)**2)
           RBND_ANGLE(1) = atan(VEL_YZ(1)/particles(p) % vel(3))*180.0/PI

           ! Calculate the particle-wall rebound angle.
           P = LC
           VEL_YZ(2) = sqrt(particles(p) % vel(1)**2 + particles(p) % vel(2)**2)
           RBND_ANGLE(2) = atan(VEL_YZ(2)/particles(p) % vel(3))*180.0/PI

          ! Write the results to a file.
          write(UDF_UNIT,"(3(3x,F11.4))") INIT_ANGLE(P), RBND_ANGLE(1:2)
       end do 

       close(UDF_UNIT)

       ! Particle-wall Tangential Restitution Coefficient
       !---------------------------------------------------------------------//
       ! Open the files.
       FNAME = 'POST_COEFF.dat'
       open(UNIT=UDF_UNIT,FILE=FNAME, POSITION="APPEND",STATUS='OLD')

       do LC = 2, 31
           ! Calculate the particle-wall restitution coefficient.
           P = LC+31
           VEL_YZ(1) = sqrt(particles(p) % vel(1)**2 + particles(p) % vel(2)**2)
           RST_COEFF(1) =  VEL_YZ(1)/INIT_VEL_T(p)

           ! Calculate the particle-wall restitution coefficient.
           P=LC
           VEL_YZ(2) = sqrt(particles(p) % vel(1)**2 + particles(p) % vel(2)**2)
           RST_COEFF(2) =  VEL_YZ(2)/INIT_VEL_T(p)

           ! Write the results to a file.
           write(UDF_UNIT,"(3(3x,F11.4))") INIT_ANGLE(P), RST_COEFF(1:2)
       end do

       close(UDF_UNIT)

       ! Particle-wall Angular velocity (rad/sec)
       !---------------------------------------------------------------------//
       ! Open the files.
       FNAME = 'POST_OMEGA.dat'
       open(UNIT=UDF_UNIT,FILE=FNAME, POSITION="APPEND",STATUS='OLD')

       do LC = 1, 31
           ! Calculate particle-particle angular velocity
           P=LC+31
           ROT_YZ(1) = sqrt(particles(p) % omega(1)**2 + particles(p) % omega(2)**2)

           ! Calculate particle-wall angular velocity
           P=LC
           ROT_YZ(2) = sqrt(particles(p) % omega(1)**2 + particles(p) % omega(2)**2)

           ! Write the results to a file.
           write(UDF_UNIT,"(3(3x,F11.4))") INIT_ANGLE(P), ROT_YZ(1:2)
       enddo
       close(UDF_UNIT)
   end if
   return
end subroutine USR3_DES
