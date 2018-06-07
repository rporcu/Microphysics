   module adjust_a

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: adjust_a                                               !
!                                                                      !
!  Purpose: Handle the special case of the center coefficient in       !
!           momentum eqs becoming zero.                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

contains

  subroutine adjust_a_u(slo, shi, alo, ahi, lo, hi, A_m, b_m, rop_g, dy, dz)

    use param, only: one, zero, small_number

    implicit none

    integer(c_int), intent(in   ) :: lo(3),hi(3)
    integer(c_int), intent(in   ) :: alo(3),ahi(3)
    integer(c_int), intent(in   ) :: slo(3),shi(3)

    real(rt), intent(inout) :: &
         A_m (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3), &
         b_m (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

    real(rt), intent(in   ) :: &
         rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(rt), intent(in   ) :: dy, dz
!---------------------------------------------------------------------//

    integer      :: i, j, k
    real(rt) :: half_ayz

    half_ayz = 0.5d0*dy*dz

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (abs(A_m(i,j,k,0)) < small_number) then

                A_m(i,j,k,:) = zero
                A_m(i,j,k,0) = -one

                if (b_m(i,j,k) < zero) then

                   b_m(i,j,k) = sqrt(abs(b_m(i,j,k))/(half_ayz*rop_g(i+1,j,k)))

                else if (b_m(i,j,k) > zero) then

                   b_m(i,j,k) = sqrt(abs(b_m(i,j,k))/(half_ayz*rop_g(i  ,j,k)))

                end if

             end if
          end do
       end do
    end do

  end subroutine adjust_a_u


  subroutine adjust_a_v(slo, shi, alo, ahi, lo, hi, A_m, b_m, rop_g, dx, dz)

    use param, only: one, zero, small_number

    implicit none

    integer(c_int), intent(in   ) :: lo(3),hi(3)
    integer(c_int), intent(in   ) :: alo(3),ahi(3)
    integer(c_int), intent(in   ) :: slo(3),shi(3)

    real(rt), intent(inout) :: &
         A_m (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3), &
         b_m (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

    real(rt), intent(in   ) :: &
         rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(rt), intent(in   ) :: dx, dz
!---------------------------------------------------------------------//

    integer      :: i, j, k
    real(rt) :: half_axz

    half_axz = 0.5d0*dx*dz

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (abs(A_m(i,j,k,0)) < small_number) then

                A_m(i,j,k,:) = zero
                A_m(i,j,k,0) = -one

                if (b_m(i,j,k) < zero) then

                   b_m(i,j,k) = sqrt(abs(b_m(i,j,k))/(half_axz*rop_g(i,j+1,k)))

                else if (b_m(i,j,k) > zero) then

                   b_m(i,j,k) = sqrt(abs(b_m(i,j,k))/(half_axz*rop_g(i,j  ,k)))

                end if

             end if
          end do
       end do
    end do

  end subroutine adjust_a_v


  subroutine adjust_a_w(slo, shi, alo, ahi, lo, hi, A_m, b_m, rop_g, dx, dy)

    use param, only: one, zero, small_number

    implicit none

    integer(c_int), intent(in   ) :: lo(3),hi(3)
    integer(c_int), intent(in   ) :: alo(3),ahi(3)
    integer(c_int), intent(in   ) :: slo(3),shi(3)

    real(rt), intent(inout) :: &
         A_m (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3), &
         b_m (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

    real(rt), intent(in   ) :: &
         rop_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(rt), intent(in   ) :: dx, dy
!---------------------------------------------------------------------//

    integer      :: i, j, k
    real(rt) :: half_axy

    half_axy = 0.5d0*dx*dy

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (abs(A_m(i,j,k,0)) < small_number) then

                A_m(i,j,k,:) = zero
                A_m(i,j,k,0) = -one

                if (b_m(i,j,k) < zero) then

                   b_m(i,j,k) = sqrt(abs(b_m(i,j,k))/(half_axy*rop_g(i,j,k+1)))

                else if (b_m(i,j,k) > zero) then

                   b_m(i,j,k) = sqrt(abs(b_m(i,j,k))/(half_axy*rop_g(i,j,k  )))

                end if

             end if
          end do
       end do
    end do

  end subroutine adjust_a_w
   end module adjust_a
