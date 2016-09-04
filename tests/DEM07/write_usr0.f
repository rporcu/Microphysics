!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR0                                             C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Purpose: This routine is called before the time loop starts and is  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_USR0

      use constant
      use desgrid
      use des_allocate
      use discretelement
      use functions
      use physprop
      use mpi_utility

      IMPLICIT NONE

      integer, parameter :: pCount = 1222

      double precision :: ru1, ru2, ru3, ru4, rand3(3)

      double precision :: x0, y0, z0, xLen, yLen, zLen

      double precision :: lRad, lDp2
      double precision :: lPos(3), lVel(3), meanVel(3), vScale

      integer :: lc1

      double precision :: gTemp, ldist(3), ldmag

! Initial granular energy (m^2/sec^2)
      double precision, parameter :: T0 = 1.0d-1
      double precision, parameter :: one_third = 1.0d0/3.0d0


      lRad = 0.5d0*D_p0(1)
      lDp2 = D_p0(1)**2

      x0 = dg_xStart + lRad
      y0 = dg_yStart + lRad
      z0 = dg_zStart + lRad

      xLen = dg_xEnd - dg_xStart - D_p0(1)
      yLen = dg_yEnd - dg_yStart - D_p0(1)
      zLen = dg_zEnd - dg_zStart - D_p0(1)

      seed_lp: do while(pip < pCount)
         call random_number(rand3)
         lPos = x0 + xLen*rand3
         do lc1=1, pip
            ldist = des_pos_new(:,lc1) - lPos
            ldmag = dot_product(ldist, ldist)
            if(ldmag + 0.001 < lDp2) cycle seed_lp
         enddo

         pip = pip+1
         call particle_grow(pip)
         call set_normal(pip)

         call random_number(ru1)
         call random_number(ru2)
         call random_number(ru3)
         call random_number(ru4)

         lVel(1) = DSQRT(-2.0d0*DLOG(DBLE(ru1)))*COS(2.0d0*PI*ru2)
         lVel(2) = DSQRT(-2.0d0*DLOG(DBLE(ru1)))*SIN(2.0d0*PI*ru2)
         lVel(3) = DSQRT(-2.0d0*DLOG(DBLE(ru3)))*COS(2.0d0*PI*ru4)

         des_pos_new(:,pip) = lPos
         des_vel_new(:,pip) = lVel

         omega_new(:,pip) = 0.0d0

         des_radius(pip) = lRad
         ro_sol(pip) = RO_s0(1)

         meanVel = meanVel + lVel

      enddo seed_lp

! Calc the average mean velocity in each direction 
      meanVel = meanVel/dble(pip)

! Subtract mean velocity from the random velocities to get a zero
! mean velocity. Also, calculate the mean granular temperature.
      gTemp = 0.0d0
      do lc1 = 1, pip
        des_vel_new(:,lc1) = des_vel_new(:,lc1) - meanVel
        gTemp = gTemp + one_third * dot_product &
           (des_vel_new(:,lc1), des_vel_new(:,lc1))
      enddo
      gTemp = gTemp/DBLE(pip)

! Scale velocities so the mean granular temperature is equal to
! the targeted valued.
      gTemp = dsqrt(T0/gTemp)

      do lc1 = 1, pip
        des_vel_new(:,lc1) = des_vel_new(:,lc1)*gTemp
      enddo

      particles = pip
      call global_all_sum(particles)

      RETURN
      END SUBROUTINE WRITE_USR0
