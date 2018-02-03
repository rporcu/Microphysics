module eb_levelset
    use amrex_fort_module, only : c_real => amrex_real
    use iso_c_binding    , only: c_int

    implicit none

contains

subroutine fill_levelset(lo,    hi,   n_refine, &
                         valid, vlo,  vhi,      &
                         phi,   phlo, phhi,     &
                         dx                   ) &
           bind(C, name="fill_levelset")

    implicit none

    integer,      dimension(3), intent(in   ) :: lo, hi, vlo, vhi, phlo, phhi
    integer,                    intent(in   ) :: n_refine
    integer,                    intent(  out) :: valid(vlo(1):vhi(1),   vlo(2):vhi(2),   vlo(3):vhi(3)  )
    real(c_real),               intent(  out) :: phi( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
    real(c_real), dimension(3), intent(in   ) :: dx

    real(c_real), dimension(3) :: pos_node, xy_centre, normal, plo
    real(c_real)               :: r_edge, h_bottom, levelset_node, phi_interp

    integer :: i, j, k, ii, jj, kk, klo, khi, jlo, jhi, ilo, ihi
    logical :: valid_cell

    xy_centre = (/ 0.0016, 0.0016, 0. /)
    r_edge    = 0.0016
    h_bottom  = 1e-8;

    klo = lo(3) - n_refine
    khi = hi(3) + n_refine
    jlo = lo(2) - n_refine
    jhi = hi(2) + n_refine
    ilo = lo(1) - n_refine
    ihi = hi(1) + n_refine
    
    do kk = klo, khi
        do jj = jlo, jhi
            do ii = ilo, ihi
                pos_node        = (/ ii*dx(1)/n_refine, jj*dx(2)/n_refine, kk*dx(3)/n_refine /)
                levelset_node   = levelset_cylinder(pos_node, xy_centre, r_edge, h_bottom)
                phi(ii, jj, kk) = levelset_node
                if ( levelset_node .le. 0 ) then
                    valid(ii, jj, kk) = 1
                end if
            end do
        end do
    end do

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                valid_cell = neighbour_is_valid(phi, phlo, phhi, i, j, k, n_refine)

                if ( valid_cell ) then
                    valid(i, j, k) = 1
                end if
            end do
        end do
    end do

contains

    pure function neighbour_is_valid(phi, phlo, phhi, i, j, k, n_refine)
        implicit none

        ! ** output type
        logical :: neighbour_is_valid

        ! ** input types
        integer,      dimension(3), intent(in) :: phlo, phhi
        real(c_real),               intent(in) :: phi( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
        integer,                    intent(in) :: i, j, k, n_refine

        ! ** declare local variables
        ! ii, jj, kk : loop variables itterating over neighbour stencil
        ! klo ... ihi: boundaries of stencil which will be checked for valid cells
        !              a cell is valid if phi <= 0
        integer :: ii, jj, kk, klo, khi, jlo, jhi, ilo, ihi

        ! neighbour stencil of size n_refine
        klo = k - 2 * n_refine
        if ( klo .lt. phlo(3) ) then
            klo = klo + n_refine
        end if

        khi = k + 2 * n_refine
        if ( khi .gt. phhi(3) ) then
            khi = khi - n_refine
        end if

        jlo = j - 2 * n_refine
        if ( jlo .lt. phlo(2) ) then
            jlo = jlo + n_refine
        end if

        jhi = j + 2 * n_refine
        if ( jhi .gt. phhi(2) ) then
            jhi = jhi - n_refine
        end if

        ilo = i - 2 * n_refine
        if ( ilo .lt. phlo(1) ) then
            ilo = ilo + n_refine
        end if

        ihi = i + 2 * n_refine
        if ( ihi .gt. phhi(1) ) then
            ihi = ihi - n_refine
        end if

        neighbour_is_valid = .false.
        
        do kk = klo, khi
            do jj = jlo, jhi
                do ii = ilo, ihi
                    if ( phi(ii, jj, kk) .le. 0 ) then
                        neighbour_is_valid = .true.
                        return
                    end if
                end do
            end do
        end do

    end function neighbour_is_valid

    pure function value_wall(pos, centre, r_edge)
        implicit none

        ! output type
        real(c_real) :: value_wall

        ! input types
        real(c_real), dimension(3), intent(in) :: pos, centre
        real(c_real),               intent(in) :: r_edge

        ! internal types
        real(c_real), dimension(2) :: r

        ! vector from centre (cylinder cylinder independent of z-value)
        r(1:2) = pos(1:2) - centre(1:2)

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

    pure function levelset_cylinder(pos, centre, r_edge, h_bottom)
        implicit none

        ! output type
        real(c_real) :: levelset_cylinder

        ! input types
        real(c_real), dimension(3), intent(in) :: pos, centre
        real(c_real),               intent(in) :: r_edge, h_bottom

        levelset_cylinder = min( value_wall(pos, centre, r_edge), &
                                 value_bottom(pos, h_bottom)      )

    end function levelset_cylinder

end subroutine fill_levelset

subroutine interp_levelset(pos, plo,  n_refine, &
                           phi, phlo, phhi,     &
                           dx,  phi_interp    ) &
           bind(C, name="interp_levelset")


    implicit none

    real(c_real), dimension(3), intent(in   ) :: pos, plo
    integer,      dimension(3), intent(in   ) :: phlo, phhi
    integer,                    intent(in   ) :: n_refine
    real(c_real),               intent(in   ) :: phi( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
    real(c_real), dimension(3), intent(in   ) :: dx
    real(c_real),               intent(  out) :: phi_interp

    integer                    :: i, j, k
    real(c_real)               :: x, y, z, xp, yp, zp, lx, ly, lz, wx_lo, wx_hi, wy_lo, wy_hi, wz_lo, wz_hi
    real(c_real), dimension(3) :: inv_dx

    real(c_real) :: inv_norm

    inv_dx(:) = n_refine / dx(:)

    xp = pos(1) - plo(1)
    yp = pos(2) - plo(2)
    zp = pos(3) - plo(3)

    lx = xp * inv_dx(1)
    ly = yp * inv_dx(2)
    lz = zp * inv_dx(3)

    i = floor(lx)
    j = floor(ly)
    k = floor(lz)

    wx_hi = lx - i
    wy_hi = ly - j
    wz_hi = lz - k

    wx_lo = 1.0d0 - wx_hi
    wy_lo = 1.0d0 - wy_hi
    wz_lo = 1.0d0 - wz_hi

    phi_interp = phi(i,   j,   k  ) * wx_lo * wy_lo * wz_lo &
               + phi(i+1, j,   k  ) * wx_hi * wy_lo * wz_lo &
               + phi(i,   j+1, k  ) * wx_lo * wy_hi * wz_lo &
               + phi(i,   j,   k+1) * wx_lo * wy_lo * wz_hi &
               + phi(i+1, j+1, k  ) * wx_hi * wy_hi * wz_lo &
               + phi(i,   j+1, k+1) * wx_lo * wy_hi * wz_hi &
               + phi(i+1, j,   k+1) * wx_hi * wy_lo * wz_hi &
               + phi(i+1, j+1, k+1) * wx_hi * wy_hi * wz_hi

    !phi_interp = -phi_interp
    !if( phi_interp .le. 0.) then
    !    write(*,*) i,j,k,  phi_interp
    !end if 
end subroutine interp_levelset

subroutine normal_levelset(pos, plo,   n_refine, &
                           phi, phlo,  phhi,     &
                           dx,  normal         ) &
           bind(C, name="normal_levelset")

    implicit none

    real(c_real), dimension(3), intent(in   ) :: pos, plo
    integer,      dimension(3), intent(in   ) :: phlo, phhi
    integer,                    intent(in   ) :: n_refine
    real(c_real),               intent(in   ) :: phi( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
    real(c_real), dimension(3), intent(in   ) :: dx
    real(c_real), dimension(3), intent(  out) :: normal

    integer                    :: i, j, k
    real(c_real)               :: x, y, z, xp, yp, zp, lx, ly, lz, wx_lo, wx_hi, wy_lo, wy_hi, wz_lo, wz_hi
    real(c_real), dimension(3) :: inv_dx

    real(c_real) :: inv_norm

    inv_dx = n_refine / dx

    xp = pos(1) - plo(1)
    yp = pos(2) - plo(2)
    zp = pos(3) - plo(3)

    lx = xp * inv_dx(1)
    ly = yp * inv_dx(2)
    lz = zp * inv_dx(3)

    i = floor(lx)
    j = floor(ly)
    k = floor(lz)

    wx_hi = lx - i
    wy_hi = ly - j
    wz_hi = lz - k

    wx_lo = 1.0d0 - wx_hi
    wy_lo = 1.0d0 - wy_hi
    wz_lo = 1.0d0 - wz_hi

    normal(1) = - phi(i,   j,   k  )*inv_dx(1) * wy_lo * wz_lo  &
                + phi(i+1, j,   k  )*inv_dx(1) * wy_lo * wz_lo  &
                - phi(i,   j+1, k  )*inv_dx(1) * wy_hi * wz_lo  &
                + phi(i+1, j+1, k  )*inv_dx(1) * wy_hi * wz_lo  &
                - phi(i,   j,   k+1)*inv_dx(1) * wy_lo * wz_hi  &
                + phi(i+1, j,   k+1)*inv_dx(1) * wy_lo * wz_hi  &
                - phi(i,   j+1, k+1)*inv_dx(1) * wy_hi * wz_hi  &
                + phi(i+1, j+1, k+1)*inv_dx(1) * wy_hi * wz_hi

    normal(2) = - phi(i,   j,   k  )*inv_dx(2) * wx_lo * wz_lo  &
                + phi(i,   j+1, k  )*inv_dx(2) * wx_lo * wz_lo  &
                - phi(i+1, j,   k  )*inv_dx(2) * wx_hi * wz_lo  &
                + phi(i+1, j+1, k  )*inv_dx(2) * wx_hi * wz_lo  &
                - phi(i,   j,   k+1)*inv_dx(2) * wx_lo * wz_hi  &
                + phi(i,   j+1, k+1)*inv_dx(2) * wx_lo * wz_hi  &
                - phi(i+1, j,   k+1)*inv_dx(2) * wx_hi * wz_hi  &
                + phi(i+1, j+1, k+1)*inv_dx(2) * wx_hi * wz_hi

    normal(3) = - phi(i,   j,   k  )*inv_dx(3) * wx_lo * wy_lo  &
                + phi(i,   j,   k+1)*inv_dx(3) * wx_lo * wy_lo  &
                - phi(i+1, j,   k  )*inv_dx(3) * wx_hi * wy_lo  &
                + phi(i+1, j,   k+1)*inv_dx(3) * wx_hi * wy_lo  &
                - phi(i,   j+1, k  )*inv_dx(3) * wx_lo * wy_hi  &
                + phi(i,   j+1, k+1)*inv_dx(3) * wx_lo * wy_hi  &
                - phi(i+1, j+1, k  )*inv_dx(3) * wx_hi * wy_hi  &
                + phi(i+1, j+1, k+1)*inv_dx(3) * wx_hi * wy_hi

    ! this might not be necessary if the phi grid is dense enough...
    !inv_norm = 1.0d0 / sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)
    !normal(:) = normal(:) * inv_norm

end subroutine normal_levelset


end module eb_levelset
