!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR1 (L)                                         C
!  Author:                                            Date: dd-mmm-yy  C
!                                                                      C
!  Purpose: Write user-defined output                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine write_usr1(l, time, dt, max_pip, des_pos_new, des_vel_new, omega_new)

      use bl_fort_module, only : c_real
      use usr, only: UPDATE_RK4_SOL

      IMPLICIT NONE

      integer,      intent(in   ) :: l, max_pip
      real(c_real), intent(in   ) :: time, dt
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)
      real(c_real), intent(in   ) :: omega_new(max_pip,3)

      SELECT CASE(L)
      CASE(1); CALL WRITE_TEST_DATA(max_pip, des_vel_new)
      END SELECT

      RETURN
      END SUBROUTINE WRITE_USR1

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_TEST_DATA(max_pip, des_vel_new)

      use run, only: time
      use bl_fort_module, only : c_real

      implicit none

      integer,      intent(in   ) :: max_pip
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)
      integer :: lc1, lc2
      real(c_real) :: ltime, gTemp

      integer, parameter :: fUnit=2030

      OPEN(UNIT=fUnit,FILE='POST_GRAN_TEMP.dat',&
         POSITION="APPEND",STATUS='OLD')

! This is a test of the setup
      lc2 = 0
      gTemp = 0.0d0
      do lc1=1, max_pip
      !   if(normal_particle==particle_state(lc1)) then
            gTemp = gTemp + dot_product &
               (des_vel_new(lc1,:),des_vel_new(lc1,:))
            lc2 = lc2 + 1
      !   endif
      enddo

      gTemp = gTemp/(3.0d0*DBLE(lc2))
      ltime = time*sqrt(0.1d0)/100.0e-6

      write(fUnit,"(3(2x,es13.6))") ltime, ToT0(time), gTemp/0.1d0

      close(fUnit)

      contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      real(c_real) function ToT0(ptime)

      real(c_real), intent(in) :: pTime

      real(c_real), parameter :: A = 222.511014d0
      real(c_real), parameter :: B =  91.382708d0

      real(c_real) :: expBto2

      expBto2 = exp(B*pTime/2.0d0)

      ToT0 = 1.0d0/(expBto2 + (A/B)*sqrt(0.10d0)*(expBto2-1.0d0))
      ToT0 = ToT0**2

      end function ToT0

      END SUBROUTINE WRITE_TEST_DATA
