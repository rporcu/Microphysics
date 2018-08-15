!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: mfix_fixup_eb_geom                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine mfix_fixup_eb_geom(lo, hi, flag, flo, fhi, vfrac, vlo, vhi, domlo, domhi) &
       bind(c,name='mfix_fixup_eb_geom')

    use amrex_fort_module, only : rt=>amrex_real

    use amrex_ebcellflag_module, only : is_regular_cell, is_single_valued_cell,  &
                                        get_neighbor_cells, set_regular_cell
    integer, dimension(3), intent(in) :: lo, hi, flo, fhi, vlo, vhi, domlo, domhi
    integer, intent(inout) :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(rt), intent(in)  :: vfrac(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

    integer :: i,j,k, nbr(-1:1,-1:1,-1:1)
    real(rt), parameter :: almostone = 1.d0 - 1.d-13

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (is_single_valued_cell(flag(i,j,k)) .and. vfrac(i,j,k) .gt. almostone) then
                call get_neighbor_cells(flag(i,j,k),nbr)
                if ( (nbr(-1,0,0).eq.1 .or. i.eq.domlo(1)) .and. &
                     (nbr( 1,0,0).eq.1 .or. i.eq.domhi(1)) .and. &
                     (nbr(0,-1,0).eq.1 .or. j.eq.domlo(2)) .and. &
                     (nbr(0, 1,0).eq.1 .or. j.eq.domhi(2)) .and. &
                     (nbr(0,0,-1).eq.1 .or. k.eq.domlo(3)) .and. &
                     (nbr(0,0, 1).eq.1 .or. k.eq.domhi(3)) ) then
                   call set_regular_cell(flag(i,j,k))
                end if
             end if

          end do
       end do
    end do

end subroutine mfix_fixup_eb_geom

