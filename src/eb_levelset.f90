module eb_levelset
    use amrex_fort_module, only : c_real => amrex_real
    use iso_c_binding    , only: c_int
    implicit none

contains

subroutine compute_levelset(lo, hi,             &
                            flag,  fglo, fghi,  &
                            valid, vlo,  vhi,   &
                            phi,   phlo, phhi,  &
                            dx                ) &
           bind(C, name="compute_levelset")

    implicit none

    integer,      dimension(3), intent(in   ) :: lo, hi, fglo, fghi, vlo, vhi, phlo, phhi
    integer,                    intent(in   ) :: flag(  fglo(1):fghi(1), fglo(2):fghi(2), fglo(3):fghi(3) )
    integer,                    intent(  out) :: valid(  vlo(1): vhi(1),  vlo(2): vhi(2),  vlo(3): vhi(3) )
    real(c_real),               intent(  out) :: phi(   phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
    real(c_real), dimension(3), intent(in   ) :: dx

    real(c_real), dimension(3) :: pos_node, xy_centre
    real(c_real)               :: r_edge, h_bottom, levelset_node

    integer :: i, j, k, ii, jj, kk, klo, khi, jlo, jhi, ilo, ihi
    logical :: valid_cell

    xy_centre = (/ 0.0016, 0.0016, 0. /)
    r_edge    = 0.0016
    h_bottom  = 1e-6;
    
    klo = lo(3)-1
    khi = hi(3)+1
    jlo = lo(2)-1
    jhi = hi(2)+1
    ilo = lo(1)-1
    ihi = hi(1)+1

    do kk = klo, khi
        do jj = jlo, jhi
            do ii = ilo, ihi
                pos_node        = (/ ii * dx(1), jj * dx(2), kk * dx(3) /)
                levelset_node   = levelset_cylinder(pos_node, xy_centre, r_edge, h_bottom)
                phi(ii, jj, kk) = levelset_node
                if ( levelset_node .le. 0 ) then
                    valid(ii, jj, kk) = 1
                else
                    valid(ii, jj, kk) = 0
                end if
                end do
        end do
    end do

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                
                klo = k-1
                khi = k+1
                jlo = j-1
                jhi = j+1
                ilo = i-1
                ihi = i+1
                
                do kk = klo, khi
                    do jj = jlo, jhi
                        do ii = ilo, ihi
                            if ( phi(ii, jj, kk) .le. 0 ) then
                                valid_cell = .true.
                            end if
                        end do
                    end do
                end do

                if ( valid_cell ) then
                    valid(i, j, k) = 1
                    ! write(*,*) i, j, k, "phi=", phi(i, j, k)
                end if

                valid_cell = .false.
            end do
        end do
    end do

contains

    function value_wall(pos, centre, r_edge)
        implicit none
        
        ! output type
        real(c_real) :: value_wall
        
        ! input types
        real(c_real), dimension(3), intent(in) :: pos, centre
        real(c_real),               intent(in) :: r_edge
        
        ! internal type
        real(c_real), dimension(2) :: r
        
        ! vector from centre (cylinder cylinder independent of z-value)
        r(1:2) = pos(1:2) - centre(1:2)
        
        ! write(*,*) "pos=", pos, "r=", r, "levelset=", r_edge - sqrt( r(1)**2 + r(2)**2 )

        value_wall = r_edge - sqrt( r(1)**2 + r(2)**2 )

    end function value_wall

    pure function value_bottom(pos, bottom_h)
        implicit none

        ! output type
        real(c_real) :: value_bottom

        ! input types
        real(c_real), dimension(3), intent(in) :: pos
        real(c_real),               intent(in) :: bottom_h

        value_bottom = pos(3) - bottom_h

    end function value_bottom

    function levelset_cylinder(pos, centre, r_edge, h_bottom)
        implicit none

        ! output type
        real(c_real) :: levelset_cylinder

        ! input types
        real(c_real), dimension(3), intent(in) :: pos, centre
        real(c_real),               intent(in) :: r_edge, h_bottom

        levelset_cylinder = min( value_wall(pos, centre, r_edge), &
                                 value_bottom(pos, h_bottom)      )

    end function levelset_cylinder

end subroutine compute_levelset

subroutine normal_levelset(pos,   plo,         &
                           lo,    hi,          &
                           valid, vlo,  vhi,   &
                           phi,   phlo, phhi,  &
                           dx,    normal     ) &
           bind(C, name="normal_levelset")

    implicit none

    real(c_real), dimension(3), intent(in   ) :: pos, plo
    integer,      dimension(3), intent(in   ) :: lo, hi, fglo, fghi, vlo, vhi, phlo, phhi
    integer,                    intent(in   ) :: valid(  vlo(1): vhi(1),  vlo(2): vhi(2),  vlo(3): vhi(3) )
    real(c_real),               intent(in   ) :: phi(   phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
    real(c_real), dimension(3), intent(in   ) :: dx
    real(c_real), dimension(3), intent(  out) :: normal

    integer                    :: i, j, k
    real(c_real)               :: xp, yp, zp, lx, ly, lz, wx_lo, wx_hi, wy_lo, wy_hi, wz_lo, wz_hi
    real(c_real), dimension(3) :: inv_dx

    inv_dx = 1.0d0 / dx

    xp = pos(1) - plo(1)
    yp = pos(2) - plo(2)
    zp = pos(3) - plo(3)

    lx = xp * inv_dx(1)
    ly = yp * inv_dx(2)
    lz = zp * inv_dx(3)

    i = floor(lx)
    j = floor(ly)
    z = floor(lz)

    wx_hi = lx - i
    wy_hi = ly - j
    wz_hi = lz - k

    wx_lo = 1.0d0 - wx_hi
    wy_lo = 1.0d0 - wy_hi
    wz_lo = 1.0d0 - wz_hi

end subroutine normal_levelset


end module eb_levelset
