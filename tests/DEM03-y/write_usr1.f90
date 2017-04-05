!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
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
                                  ylength)
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
      Use usr, only: gx1, gx2, gy1, gy2, rk4_v4
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
         POSITION="APPEND",STATUS='OLD')

      OPEN(UNIT=uPos2,FILE='POST_POS2.dat', &
         POSITION="APPEND",STATUS='OLD')

! Calculate the value for the RK4 solutions.
      TIME_INTERVAL = lTime - RK4_TIME
      IF(TIME_INTERVAL .LE. RK4_DT) THEN
         RK4_STEPS = 1
         RK4_DT = TIME_INTERVAL
         RK4_DT_LAST = UNDEFINED
      ELSE
         RK4_STEPS = floor(real(TIME_INTERVAL/RK4_DT_DEFAULT))
         RK4_DT = RK4_DT_DEFAULT
         RK4_DT_LAST = lTime - (RK4_TIME + RK4_STEPS*RK4_DT)
      ENDIF

      DO I=1, RK4_STEPS
         CALL RK4_V4(RK4_DT, gY1, gX1, gY2, gX2)
         RK4_TIME = RK4_TIME + RK4_DT
      ENDDO

      IF(IS_DEFINED(RK4_DT_LAST)) THEN
         CALL RK4_V4(RK4_DT_LAST, gY1, gX1, gY2, gX2, length)
         RK4_TIME = RK4_TIME + RK4_DT_LAST
      ENDIF

! Write the results to a file.
      WRITE(uPos1,"(3x,F15.8,5X,F15.8,2(3x,F15.8))") lTime, gY1,   &
         DES_POS_new(1,2), (ABS(gY1 - DES_POS_new(1,2))/ABS(gY1))*100

      WRITE(uPos2,"(3x,F15.8,5X,F15.8,2(3x,F15.8))") lTime, gY2,   &
         DES_POS_new(2,2), (ABS(gY2 - DES_POS_new(2,2))/ABS(gY1))*100

      CLOSE(uPos1)
      CLOSE(uPos2)

      end subroutine write_des_out
