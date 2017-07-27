module calc_trd_g_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_trD_g                                             !
!  Purpose: Calculate the trace of gas phase rate of strain tensor     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine calc_trd_g(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi,&
       trd_g, u_g, v_g, w_g, dx, dy, dz) &
       bind(C, name="calc_trd_g")

    implicit none

    integer(c_int), intent(in   ) :: slo(3),shi(3)
    integer(c_int), intent(in   ) ::  lo(3), hi(3)
    integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
    real(c_real),   intent(in   ) :: dx, dy, dz

! Dummy Arguments
!-----------------------------------------------
    real(c_real), intent(  out) :: &
         trd_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(c_real), intent(in   ) :: &
         u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
         v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
         w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

! Local variables
!-----------------------------------------------
    integer :: i, j, k
    real(c_real) :: odx, ody, odz

    odx = 1.d0 / dx
    ody = 1.d0 / dy
    odz = 1.d0 / dz

    print *,"VG IN TRD ",v_g(0,49,0)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             trd_g(i,j,k) = &
                  (u_g(i+1,j,k)-u_g(i,j,k))*odx + &
                  (v_g(i,j+1,k)-v_g(i,j,k))*ody + &
                  (w_g(i,j,k+1)-w_g(i,j,k))*odz
          end do
       end do
    end do

  end subroutine calc_trd_g
end module calc_trd_g_module
