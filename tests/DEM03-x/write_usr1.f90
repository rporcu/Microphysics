
!                                                                      C
!  Module name: write_usr1 (L)                                         C
!  Purpose: Write user-defined output                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine write_usr1(l, time, dt, max_pip, des_pos_new, des_vel_new, omega_new,
                            xlength, ylength, zlength )

      use amrex_fort_module, only : c_real => amrex_real

      IMPLICIT NONE

      integer,      intent(in   ) :: l, max_pip
      real(c_real), intent(in   ) :: time, dt, xlength, ylength, zlength 
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)
      real(c_real), intent(in   ) :: omega_new(max_pip,3)

      SELECT CASE(L)
      CASE(1); call write_des_out(TIME, max_pip, des_pos_new, des_vel_new, &
                                  xlength)
      end SELECT

      end subroutine write_usr1

!......................................................................!
!  Subroutine: write_des_out                                           !
!                                                                      !
!  Purpose: Calculate the position and velocity (1D Y-axis) of a free  !
!  falling particle. Compare the results to the MFIX-DEM solultion.    !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!......................................................................!
      subroutine write_des_out(lTime, max_pip, des_pos_new, des_vel_new, length)

      use amrex_fort_module, only : c_real => amrex_real
      Use usr, only: gz1, gz2, gx1, gx2, rk4_v4
      Use param1, only: undefined, is_defined

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      integer,      intent(in   ) ::  max_pip
      real(c_real), intent(in   ) :: ltime, length
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)

! Local variables
!---------------------------------------------------------------------//
! file unit for heat transfer data
      integer, parameter :: upos1 = 2030
      integer, parameter :: upos2 = 2031

      real(c_real), save :: rk4_time = 0.0d0
      real(c_real) :: rk4_dt, rk4_dt_last
      real(c_real) time_interval
      real(c_real), parameter :: rk4_dt_default = 1.0d-6
      integer :: i, rk4_steps

! Open the files.
      OPEN(UNIT=uPos1,FILE='POST_POS1.dat', &
         POSITION="APPend",STATUS='OLD')

      OPEN(UNIT=uPos2,FILE='POST_POS2.dat', &
         POSITION="APPend",STATUS='OLD')

! Calculate the value for the rk4 solutions.
      TIME_INTERVAL = lTime - rk4_TIME
      IF(TIME_INTERVAL .LE. rk4_DT) THEN
         rk4_STEPS = 1
         rk4_DT = TIME_INTERVAL
         rk4_DT_LAST = UNDEFINED
      ELSE
         rk4_STEPS = floor(real(TIME_INTERVAL/rk4_DT_DEFAULT))
         rk4_DT = rk4_DT_DEFAULT
         rk4_DT_LAST = lTime - (rk4_TIME + rk4_STEPS*rk4_DT)
      endif

      do I=1, rk4_STEPS
         CALL rk4_V4(rk4_DT, gx1, gz1, gx2, gz2)
         rk4_TIME = rk4_TIME + rk4_DT
      enddo

      IF(IS_DEFINED(rk4_DT_LAST)) THEN
         CALL rk4_V4(rk4_DT_LAST, gx1, gz1, gx2, gz2, length)
         rk4_TIME = rk4_TIME + rk4_DT_LAST
      endif

      ! Write the results to a file.
      write(uPos1,"(3x,F15.8,5X,F15.8,2(3x,F15.8))") lTime, gx1,   &
         DES_POS_new(1,1), (ABS(gx1 - DES_POS_new(1,1))/ABS(gx1))*100

      write(uPos2,"(3x,F15.8,5X,F15.8,2(3x,F15.8))") lTime, gx2,   &
         DES_POS_new(2,1), (ABS(gx2 - DES_POS_new(2,1))/ABS(gx1))*100

      CLOSE(uPos1)
      CLOSE(uPos2)

      end subroutine write_des_out
