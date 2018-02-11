module eb_levelset
    use amrex_fort_module, only : c_real => amrex_real
    use iso_c_binding,    only: c_int


    use param,                   only: small_number, zero, one
    use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, is_single_valued_cell, &
                                        get_neighbor_cells

    implicit none

contains

subroutine init_levelset(lo,    hi,   n_refine, &
                         valid, vlo,  vhi,      &
                         phi,   phlo, phhi    ) &
           bind(C, name="init_levelset")

    implicit none

    integer,      dimension(3), intent(in   ) :: lo, hi, vlo, vhi, phlo, phhi
    integer,                    intent(in   ) :: n_refine
    integer,                    intent(  out) :: valid ( vlo(1):vhi(1),  vlo(2):vhi(2),  vlo(3):vhi(3) )
    real(c_real),               intent(  out) :: phi   (phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))


    integer :: i, j, k, ii, jj, kk, klo, khi, jlo, jhi, ilo, ihi
    
    
    klo = lo(3) - n_refine
    khi = hi(3) + n_refine
    jlo = lo(2) - n_refine
    jhi = hi(2) + n_refine
    ilo = lo(1) - n_refine
    ihi = hi(1) + n_refine

    do kk = klo, khi
        do jj = jlo, jhi
            do ii = ilo, ihi
                phi(ii, jj, kk)   = huge(phi)
                valid(ii, jj, kk) = -1
            end do
        end do
    end do

end subroutine init_levelset

subroutine fill_levelset(lo,    hi,   n_refine, &
                         norm,  nlo,  nhi,      &
                         valid, vlo,  vhi,      &
                         phi,   phlo, phhi,     &
                         dx                   ) &
           bind(C, name="fill_levelset")

    implicit none

    integer,      dimension(3), intent(in   ) :: lo, hi, nlo, nhi, vlo, vhi, phlo, phhi
    integer,                    intent(in   ) :: n_refine
    integer,                    intent(  out) :: valid ( vlo(1):vhi(1),  vlo(2):vhi(2),  vlo(3):vhi(3) )
    real(c_real),               intent(inout) :: phi   (phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))
    real(c_real),               intent(in   ) :: norm  ( nlo(1):nhi(1),  nlo(2):nhi(2),  nlo(3):nhi(3) )
    real(c_real), dimension(3), intent(in   ) :: dx

    real(c_real), dimension(3) :: pos_node, xy_centre, normal, plo
    real(c_real)               :: r_edge, h_bottom, levelset_node, phi_interp

    integer :: i, j, k, ii, jj, kk, klo, khi, jlo, jhi, ilo, ihi
    logical :: valid_cell

    write(*,*) lo, hi, dx

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
                levelset_node   = closest_dist(norm, nlo, nhi, pos_node, n_refine, dx)
                !if ( dabs(levelset_node) < 5e-5 ) then
                !    write(*,*) "pt=", pos_node, "levelset=", levelset_node
                !end if
                if ( levelset_node .lt. phi(ii, jj, kk) ) then
                    phi(ii, jj, kk) = levelset_node
                    if ( levelset_node .le. 0 ) then
                        valid(ii, jj, kk) = 1
                    end if
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

    function closest_dist(eb_data, nlo,  nhi,  &
                          pos,     n_refine, dx)

        implicit none

        real(c_real) :: closest_dist

        integer,                    intent(in) :: n_refine
        integer,      dimension(3), intent(in) :: nlo, nhi
        real(c_real), dimension(3), intent(in) :: pos, dx
        real(c_real),               intent(in) :: eb_data ( nlo(1):nhi(1), nlo(2):nhi(2), nlo(3):nhi(3), 6)

        integer                    :: ii, jj, kk, klo, khi, jlo, jhi, ilo, ihi
        integer,      dimension(3) :: vindex_pt, vindex_loop
        real(c_real)               :: fudge, dist_proj, dist2, min_dist2
        real(c_real), dimension(3) :: inv_dx, eb_norm, eb_cent, eb_min_pt
        logical                    :: proj_valid

        inv_dx(:)    = n_refine / dx(:)
        fudge        = one - 1.d-8
        closest_dist = huge(closest_dist)
        min_dist2    = huge(min_dist2)

        klo = nlo(3)
        khi = nhi(3)
        jlo = nlo(2)
        jhi = nhi(2)
        ilo = nlo(1)
        ihi = nhi(1)
        
        proj_valid = .false.

        do kk = klo, khi
            do jj = jlo, jhi
                do ii = ilo, ihi
                        
                        eb_cent(:) = eb_data(ii, jj, kk, 1:3)
                        eb_norm(:) = eb_data(ii, jj, kk, 3:5)

                        dist_proj  = dot_product( pos(:) - eb_cent(:), -eb_norm(:) )
                        eb_min_pt(:) = pos(:) + eb_norm(:) * dist_proj * fudge
                        vindex_pt = floor( eb_min_pt(:) * inv_dx(:) )
                        vindex_loop = floor( eb_cent(:) * inv_dx(:) )

                        !write(*,*) vindex_loop(:), eb_cent(:), sqrt( &
                        !    dot_product ( ( eb_cent(1:2) - (/0.0016, 0.0016/) ), &
                        !                  ( eb_cent(1:2) - (/0.0016, 0.0016/) ) )&
                        !)


                        if ( any( vindex_loop /= vindex_pt ) ) then
                            dist2 = dot_product( pos(:) - eb_cent(:), pos(:) - eb_cent(:) )
                            if ( dist2 < min_dist2) then
                                min_dist2 = dist2
                            end if
                        else
                            if ( dist_proj < closest_dist ) then
                                closest_dist = dist_proj
                                proj_valid   = .true.
                            end if
                        end if
                end do
            end do
        end do

        if ( .not. proj_valid ) then
            !write(*,*) "Level-set fallback."
            closest_dist = -sqrt( min_dist2)
        end if
    end function closest_dist

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


        !----------------------------------------------------------------------------------------------------!
        ! build neighbour stencil of size 2 * n_refine                                                       !
        !                                ^^^                                                                 !
        !                  *** fudge factor: finite-sized particle can overlap with neighbouring cells       !
        ! note: stencil could be out-of-bounds (due to fudge factor) => bounds-checking                      !
        !----------------------------------------------------------------------------------------------------!

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


        !---------------------------------------------------------------------------------------------------!
        ! check members of neighbour stencil:                                                               !
        !       cell is "valid" whenever at least one cell in the neighbour stencil has a level-set phi     !
        !       less than, or equal to, 0                                                                   !
        !---------------------------------------------------------------------------------------------------!

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

end subroutine fill_levelset

subroutine count_eb_facets(lo, hi, flag, flo, fhi, n_facets) &
           bind(C, name="count_eb_facets")

    implicit none

    integer, dimension(3), intent(in   ) :: lo, hi, flo, fhi
    integer,               intent(in   ) :: flag ( flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3) )
    integer,               intent(inout) :: n_facets

    integer :: i, j, k

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                if ( is_single_valued_cell( flag(i, j, k) ) ) then
                    n_facets = n_facets + 1
                end if
            end do
        end do
    end do

end subroutine count_eb_facets

subroutine eb_as_list(lo,       hi,   n_refine, &
                      flag,     flo,  fhi,      &
                      norm,     nlo,  nhi,      &
                      bcent,    blo,  bhi,      &
                      list_out, llo,  lhi, dx,  &
                      c_facets)&
           bind(C, name="eb_as_list")

    implicit none

    integer,               intent(in   ) :: n_refine
    integer, dimension(3), intent(in   ) :: lo, hi, flo, fhi, nlo, nhi, blo, bhi, llo, lhi
    real(c_real),          intent(in   ) :: dx(3)
    integer,               intent(in   ) :: flag  ( flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3) )
    real(c_real),          intent(in   ) :: norm  ( nlo(1):nhi(1), nlo(2):nhi(2), nlo(3):nhi(3), 3)
    real(c_real),          intent(in   ) :: bcent ( blo(1):bhi(1), blo(2):bhi(2), blo(3):bhi(3), 3)
    real(c_real),          intent(  out) :: list_out ( llo(1):lhi(1), llo(2):lhi(2), llo(3):lhi(3), 6)
    integer,               intent(inout) :: c_facets

    integer                    :: i, j, k
    real(c_real), dimension(3) :: eb_cent

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                if ( is_single_valued_cell( flag(i, j, k) ) ) then
                    
                    eb_cent(:) = ( bcent(i, j, k, :)                           &
                                   + (/ dble(i), dble(j), dble(k) /)           &
                                   + (/ 0.5d0, 0.5d0, 0.5d0 /) ) * dx(:)/n_refine

                    c_facets = c_facets + 1

                    list_out(c_facets, llo(2), llo(3), 1:3) = eb_cent(:)
                    list_out(c_facets, llo(2), llo(3), 4:6) = norm(i, j, k, :)
                end if
            end do
        end do
    end do

end subroutine eb_as_list

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
