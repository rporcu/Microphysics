!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR1 (L)                                         C
!  Author:                                            Date: dd-mmm-yy  C
!                                                                      C
!  Purpose: Write user-defined output                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine write_usr1(l, np, time, dt, particles)

  use amrex_fort_module, only : c_real => amrex_real
  use particle_mod,      only: particle_t

  implicit none

  integer,          intent(in   ) :: l, np
  real(c_real),     intent(in   ) :: time, dt
  type(particle_t), intent(in   ) :: particles(np)

  select case(l)
  case(1); call write_test_data(np, particles, time)
  end select

  return

contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine write_test_data(np, particles, time)

  implicit none

  integer,          intent(in   ) :: np
  type(particle_t), intent(in   ) :: particles(np)
  real(c_real),     intent(in   ) :: time

  integer :: lc1, lc2
  real(c_real) :: ltime, gTemp

  integer, parameter :: fUnit=2030

  OPEN(UNIT=fUnit,FILE='POST_GRAN_TEMP.dat',&
       POSITION="APPEND",STATUS='UNKNOWN')

! This is a test of the setup
  lc2 = 0
  gTemp = 0.0d0
  do lc1=1, np
     gTemp = gTemp + dot_product &
          (particles(lc1) % vel, particles(lc1) % vel)
     lc2 = lc2 + 1
  enddo

  gtemp = gtemp/(3.0d0*dble(lc2))
  ltime = time*sqrt(0.1d0)/100.0e-6

  write(fUnit,"(3(2x,es13.6))") ltime, ToT0(time), gTemp/0.1d0

  close(fUnit)

end subroutine write_test_data

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

end function tot0

end subroutine write_usr1
