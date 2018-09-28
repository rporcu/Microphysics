!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DISPLAY_RESID(NIT, IER)                                !
!  Purpose: Display residuals                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine display_resid(time, dt, nit, resid)&
     bind(C, name="display_resid")

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding, only: c_int

  use error_manager, only: err_msg, flush_err_msg

  implicit none

  ! iteration number
  integer(c_int), intent(in) :: nit
  real(rt),   intent(in) :: time, dt, resid(8,2)

  if(nit == 1) then
     write (err_msg, '(/" Time = ",g12.5,"  Dt = ",g12.5)') time, dt
     call flush_err_msg(header=.false., footer=.false., log=.false.)
     write(err_msg,"(3x,'Nit',5x,'Pg',8x,'Ug',8x,'Vg',8x,'Wg')")
     call flush_err_msg(header=.false., footer=.false., log=.false.)
  endif

  write(err_msg,"(i6,4(2x,1pg8.1))") nit, resid(1:4,1)
  call flush_err_msg(header=.false., footer=.false., log=.false.)

  return

end subroutine display_resid
