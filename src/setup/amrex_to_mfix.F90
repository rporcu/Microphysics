module amrex_to_mfix_module
! _________________________________________________________________

  use amrex_fort_module, only: rt => amrex_real
  use iso_c_binding,     only: c_int, c_char

  implicit none

contains

!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_usr3 (vel_g, ulo, uhi, &
                        p_g, slo, shi, dx, dy, dz) bind(C, name="mfix_usr3")

    integer(c_int), intent(in   ) :: ulo(3),uhi(3)
    integer(c_int), intent(in   ) :: slo(3),shi(3)

    real(rt), intent(inout) :: vel_g&
        (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)
    real(rt), intent(inout) :: p_g&
        (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(rt), intent(in   ) :: dx, dy, dz

    call usr3(vel_g(ulo(1):,ulo(2):,ulo(3):,1), ulo, uhi, &
              vel_g(ulo(1):,ulo(2):,ulo(3):,2), ulo, uhi, &
              vel_g(ulo(1):,ulo(2):,ulo(3):,3), ulo, uhi, &
              p_g, slo, shi, dx, dy, dz)

  end subroutine mfix_usr3

end module amrex_to_mfix_module
