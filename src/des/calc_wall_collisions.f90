module wall_collisions
    use amrex_fort_module, only : rt => amrex_real
    use iso_c_binding    , only: c_int

    implicit none
contains

  subroutine ls_has_walls (has_wall, phi, phlo, phhi, tol) bind(c, name="ls_has_walls")

    integer,  intent(  out) :: has_wall
    integer,  intent(in   ) :: phlo(3), phhi(3)
    real(rt), intent(in   ) :: phi( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
    real(rt), intent(in   ) :: tol

    if (any(phi(phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3)) .le. tol)) then
       has_wall = 1
    else
       has_wall = 0
    end if

  end subroutine ls_has_walls

end module wall_collisions
