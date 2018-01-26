subroutine compute_levelset(lo, hi,             &
                            flag,  fglo, fghi,  &
                            valid, vlo,  vhi,   &
                            phi,   phlo, phhi ) &
           bind(C, name="compute_levelset")

    use amrex_fort_module, only : c_real => amrex_real
    use iso_c_binding    , only: c_int

    use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell

    implicit none

    integer,      dimension(3), intent(in   ) :: lo, hi, fglo, fghi, vlo, vhi, phlo, phhi
    integer,                    intent(in   ) :: flag(  fglo(1):fghi(1), fglo(2):fghi(2), fglo(3):fghi(3) )
    integer,                    intent(in   ) :: valid(  vlo(1): vhi(1),  vlo(2): vhi(2),  vlo(3): vhi(3) )
    real(c_real),               intent(  out) :: phi(   phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )

    integer :: i, j, k

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)

            end do
        end do
    end do

end subroutine compute_levelset
