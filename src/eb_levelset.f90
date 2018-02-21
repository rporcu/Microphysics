module eb_levelset
    use amrex_fort_module, only : c_real => amrex_real
    use iso_c_binding,    only: c_int


    use param,                   only: small_number, zero, one
    use amrex_ebcellflag_module, only: is_regular_cell, is_covered_cell, is_single_valued_cell, &
                                       get_neighbor_cells

    implicit none

contains

    pure subroutine init_levelset(lo,  hi,          &
                                  phi, phlo, phhi ) &
               bind(C, name="init_levelset")

        implicit none

        integer,      dimension(3), intent(in   ) :: lo, hi, phlo, phhi
        real(c_real),               intent(  out) :: phi ( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )

        integer :: i, j, k

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    phi(i, j, k) = huge(phi)
                end do
            end do
        end do

    end subroutine init_levelset


subroutine fill_levelset_eb(lo,      hi,          &
                            eb_list, l_eb,        &
                            valid,   vlo,  vhi,   &
                            phi,     phlo, phhi,  &
                            dx,      n_refine   ) &
            bind(C, name="fill_levelset_eb")

    implicit none

    integer,                       intent(in   ) :: l_eb
    integer,      dimension(3),    intent(in   ) :: lo, hi, vlo, vhi, phlo, phhi
    real(c_real), dimension(l_eb), intent(in   ) :: eb_list
    real(c_real),                  intent(  out) :: phi     (phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3))
    integer,                       intent(  out) :: valid   ( vlo(1):vhi(1),   vlo(2):vhi(2),   vlo(3):vhi(3) )
    real(c_real), dimension(3),    intent(in   ) :: dx
    integer,                       intent(in   ) :: n_refine

    real(c_real), dimension(3) :: pos_node, c_vec
    real(c_real)               :: levelset_node

    integer :: ii, jj, kk
    logical :: valid_cell

    do kk = lo(3), hi(3)
        do jj = lo(2), hi(2)
            do ii = lo(1), hi(1)
                pos_node      = (/ ii*dx(1)/n_refine, jj*dx(2)/n_refine, kk*dx(3)/n_refine /)
                levelset_node = closest_dist(eb_list, l_eb, pos_node, valid_cell, dx, n_refine)
                !if ( dabs(levelset_node) < 5e-5 ) then
                !    write(*,*) ii, jj, kk, "pt=", pos_node, "valid=", valid_cell, "levelset=", levelset_node
                !end if 
                phi(ii, jj, kk) = levelset_node;
                if ( valid_cell ) then
                    valid(ii, jj, kk) = 1
                else
                    valid(ii, jj, kk) = 0
                end if
            end do
        end do
    end do

contains

    function closest_dist(eb_data, l_eb,       &
                          pos,     proj_valid, &
                          dx,      n_refine   )
        implicit none

        real(c_real) :: closest_dist

        integer,                       intent(in)  :: n_refine, l_eb
        logical,                       intent(out) :: proj_valid
        real(c_real), dimension(3),    intent(in)  :: pos, dx
        real(c_real), dimension(l_eb), intent(in)  :: eb_data

        integer                    :: i
        integer,      dimension(3) :: vi_pt_closest, vi_loop_closest
        real(c_real)               :: dist_proj, dist2, min_dist2, min_edge_dist2
        real(c_real), dimension(3) :: inv_dx, eb_norm, eb_cent, eb_min_pt, eb_cent_closest, eb_norm_closest, eb_pt_closest

        inv_dx(:)      = n_refine / dx(:)
        closest_dist   = huge(closest_dist)
        min_dist2      = huge(min_dist2)
        min_edge_dist2 = huge(min_edge_dist2)

        proj_valid = .false.
        
        do i = 1, l_eb, 6
            eb_cent(:)   = eb_data(i     : i + 2)
            eb_norm(:)   = eb_data(i + 3 : i + 5)

            dist2        = dot_product( pos(:) - eb_cent(:), pos(:) - eb_cent(:) )
            dist_proj    = dot_product( pos(:) - eb_cent(:), -eb_norm(:) )
            
            eb_min_pt(:) = pos(:) + eb_norm(:) * dist_proj
            
            if ( dist2 < min_dist2 ) then
                min_dist2          = dist2
                closest_dist       = dist_proj
                vi_loop_closest(:) = floor( eb_cent(:) * inv_dx(:) / n_refine)
                vi_pt_closest(:)   = floor( eb_min_pt(:) * inv_dx(:) / n_refine)
                eb_cent_closest(:) = eb_cent(:)
                eb_norm_closest(:) = eb_norm(:)
                eb_pt_closest(:)   = eb_min_pt(:)
            end if
        end do

        if ( all( vi_pt_closest == vi_loop_closest ) ) then
            proj_valid = .true.
        end if
        
        if ( .not. proj_valid ) then
            c_vec = facets_nearest_pt(vi_pt_closest, vi_loop_closest, pos, eb_norm_closest, eb_cent_closest)
            min_edge_dist2 = dot_product( c_vec(:) - pos(:), c_vec(:) - pos(:))
            !if ( min_dist2 < min_edge_dist2 ) then
            !    write(*,*) "ha!"
            !end if
            closest_dist = -sqrt( min(min_dist2, min_edge_dist2) )
        end if

    end function closest_dist

   !----------------------------------------------------------------------!
   !                                                                      !
   !  Pure Function: DOT_3D_REAL                                          !
   !                                                                      !
   !  Purpose: Returns the cartesian dot product for two vectors in three !
   !  dimensions.                                                         !
   !                                                                      !
   !  Comments: Vectors are represented as one-dimensional arrays of type !
   !  real(c_real) and of dimension(3) (i.e. indices range from 1..3).    !
   !----------------------------------------------------------------------!
    pure function dot_3d_real (v1, v2)
        implicit none

        real(c_real)                           :: dot_3d_real
        real(c_real), dimension(3), intent(in) :: v1, v2

        ! really naive implementation
        dot_3d_real = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

    end function dot_3d_real

   !----------------------------------------------------------------------!
   !                                                                      !
   !  Pure Function: CROSS_3D_REAL                                        !
   !                                                                      !
   !  Purpose: Returns the cartesian cross product for two vectors in     !
   !  three dimensions.                                                   !
   !                                                                      !
   !  Comments: Vectors are represented as one-dimensional arrays of type !
   !  real(c_real) and of dimension(3) (i.e. indices range from 1..3).    !
   !----------------------------------------------------------------------!
    pure function cross_3d_real (v1, v2)
        implicit none

        real(c_real), dimension(3)             :: cross_3d_real
        real(c_real), dimension(3), intent(in) :: v1, v2

        cross_3d_real(1) = v1(2)*v2(3) - v1(3)*v2(2)
        cross_3d_real(2) = v1(3)*v2(1) - v1(1)*v2(3)
        cross_3d_real(3) = v1(1)*v2(2) - v1(2)*v2(1)

    end function cross_3d_real

   !----------------------------------------------------------------------!
   !                                                                      !
   !  Pure Function: PT_IN_BOX                                            !
   !                                                                      !
   !  Purpose: Returns true if the coordinate vector represents a point   !
   !  inside a three-dimensional cube of size dx and origin specified by  !
   !  integer vector. As a final input, the user specifies a dimension    !
   !  (axis) to ignore. This allows IN-BOX checking on a box face.        !
   !                                                                      !
   !  Comments: Position vectors are represented as one-dimensional arrays!
   !  of type real(c_real) and of dimension(3) (i.e. indices range        !
   !  from 1..3). Cells are enumerated using arrays of type integer and   !
   !  dimension(3).
   !----------------------------------------------------------------------!
   ! pure function pt_in_box (pt, id, id_ignore)
   !     implicit none

   !     logical                                :: pt_in_box
   !     real(c_real), dimension(3), intent(in) :: pt
   !     integer,      dimension(3), intent(in) :: id
   !     integer,                    intent(in) :: id_ignore

   !     real(c_real), dimension(3) :: box_max, box_min
   !     integer                    :: i

   !     ! Determine box boundaries
   !     box_min(:) = id(:) * dx(:)
   !     box_max(:) = (id(:) + 1.) * dx(:)

   !     pt_in_box = .true.

   !     ! Check each coordinate. Skip ignored coordinate.
   !     do i = 1, 3
   !         if (.not. (i .eq. id_ignore) ) then
   !             if ( pt(i) .le. box_min(i) ) then
   !                 pt_in_box = .false.
   !                 exit
   !             end if

   !             if ( pt(i) .ge. box_max(i) ) then
   !                 pt_in_box = .false.
   !                 exit
   !             end if
   !         end if
   !     end do

   ! end function pt_in_box

   !----------------------------------------------------------------------!
   !                                                                      !
   !  Pure Subroutine: CALC_FACET_EDGE                                    !
   !                                                                      !
   !  Purpose: Calculates the line (represented by a position and a       !
   !  direction vector) given by the intersection of two planes (defined  !
   !  by two normal (n1, n2) and two positions (h1 = n1.p1, h2 = n2.p2).  !
   !                                                                      !
   !  When one plane is the EB surface, and the other is a face of the    !
   !  cell. Then this line represents the edge of the EB facet.           !
   !                                                                      !
   !  Comments: Vectors are represented as one-dimensional arrays of type !
   !  real(c_real) and of dimension(3) (i.e. indices range from 1..3).    !
   !----------------------------------------------------------------------!
    pure subroutine calc_facet_edge (p0, v, h1, h2, n1, n2)
        implicit none

        real(c_real), dimension(3), intent(  out) :: p0, v
        real(c_real), dimension(3), intent(in   ) :: n1, n2
        real(c_real),               intent(in   ) :: h1, h2

        real(c_real) :: c1, c2, c_dp, c_norm

        c_dp = dot_3d_real(n1, n2)
        c_norm = 1 - c_dp * c_dp

        c1 = ( h1 - h2 * c_dp ) / c_norm
        c2 = ( h2 - h1 * c_dp ) / c_norm

        p0(:) = c1 * n1(:) + c2 * n2(:)
        v = cross_3d_real(n1, n2)

    end subroutine calc_facet_edge

   !----------------------------------------------------------------------!
   !                                                                      !
   !  Pure Subroutine: LINES_NEAREST_PT                                   !
   !                                                                      !
   !  Purpose: Given an a line an a point, this subroutine finds the point!
   !  one the line which minimizes the cartesian distance. It also finds  !
   !  the corresponing distance along the line corresponding to this point!
   !                                                                      !
   !  Comments: Vectors are represented as one-dimensional arrays of type !
   !  real(c_real) and of dimension(3) (i.e. indices range from 1..3).    !
   !----------------------------------------------------------------------!
    pure subroutine lines_nearest_pt (lambda_min, nearest_pt, p0, v, pt)
        implicit none

        real(c_real),               intent(  out) :: lambda_min
        real(c_real), dimension(3), intent(  out) :: nearest_pt
        real(c_real), dimension(3), intent(in   ) :: p0, v, pt

        real(c_real), dimension(3) :: c

        c(:) = p0(:) - pt(:)
        lambda_min = - dot_3d_real(v, c) / dot_3d_real(v, v)

        nearest_pt(:) = p0(:) + lambda_min*v(:)

    end subroutine lines_nearest_pt

   !----------------------------------------------------------------------!
   !                                                                      !
   !  Pure Subroutine: SWAP_REALS                                         !
   !                                                                      !
   !  Purpose: Stupid little subroutine which swaps the values of its     !
   !  inputs.                                                             !
   !                                                                      !
   !  Comments: Inputs are of type real(c_real)                           !
   !                                                                      !
   !----------------------------------------------------------------------!
    pure subroutine swap_reals(a, b)
        implicit none

        real(c_real), intent(inout) :: a, b
        real(c_real)                :: bucket

        bucket = a
        a      = b
        b      = bucket

    end subroutine swap_reals

   !----------------------------------------------------------------------!
   !                                                                      !
   !  Pure Subroutine: LAMBDA_BOUNDS                                      !
   !                                                                      !
   !  Purpose: Given a line which passes through a box in three dimensions!
   !  (it can pass through the edges). Let lambda be a real value         !
   !  representing the coordinate along the line. This subroutine finds   !
   !  teh min/max values of lambda, in order for the point described by   !
   !  lambda to be contained within the box.                              !
   !                                                                      !
   !  Comments: Vectors are represented as one-dimensional arrays of type !
   !  real(c_real) and of dimension(3) (i.e. indices range from 1..3).    !
   !----------------------------------------------------------------------!
    pure subroutine lambda_bounds(lambda_min, lambda_max, id_cell, p0, v)
        implicit none

        real(c_real),               intent(  out) :: lambda_min, lambda_max
        integer,      dimension(3), intent(in   ) :: id_cell
        real(c_real), dimension(3), intent(in   ) :: p0, v

        ! c... are the preliminary boundaries
        real(c_real) :: cx_lo, cy_lo, cz_lo, cx_hi, cy_hi, cz_hi

        ! defaults such that if skipped, min/max will not choose these values anyway
        cx_lo = -huge(cx_lo)
        cy_lo = -huge(cy_lo)
        cz_lo = -huge(cz_lo)

        cx_hi = huge(cx_hi)
        cy_hi = huge(cy_hi)
        cz_hi = huge(cz_hi)

        ! if the line runs parrallel to any of these dimensions (which is true for
        ! EB edges), then skip -> the min/max functions at the end will skip them
        ! due to the +/-huge(c...) defaults (above).

        if ( abs(v(1)) .gt. epsilon(v) ) then
            cx_lo = -( p0(1) - dble(id_cell(1)) * dx(1) ) / v(1)
            cx_hi = -( p0(1) - ( dble(id_cell(1)) + 1. ) * dx(1) ) / v(1)

            if ( v(1) .lt. 0. ) then
                call swap_reals(cx_lo, cx_hi)
            end if
        end if

        if ( abs(v(2)) .gt. epsilon(v) ) then
            cy_lo = -( p0(2) - dble(id_cell(2)) * dx(2) ) / v(2)
            cy_hi = -( p0(2) - ( dble(id_cell(2)) + 1. ) * dx(2) ) / v(2)

            if ( v(2) .lt. 0. ) then
                call swap_reals(cy_lo, cy_hi)
            end if
        end if

        if ( abs(v(3)) .gt. epsilon(v) )  then
            cz_lo = -( p0(3) - dble(id_cell(3)) * dx(3) ) / v(3)
            cz_hi = -( p0(3) - ( dble(id_cell(3)) + 1. ) * dx(3) ) / v(3)

            if ( v(3) .lt. 0. ) then
                call swap_reals(cz_lo, cz_hi)
            end if
        endif

        lambda_min = max(cx_lo, cy_lo, cz_lo)
        lambda_max = min(cx_hi, cy_hi, cz_hi)

    end subroutine lambda_bounds

   !----------------------------------------------------------------------!
   !                                                                      !
   !  Pure Function: FACETS_NEAREST_PT                                    !
   !                                                                      !
   !  Purpose: Given a collision between particle and EB surface, and     !
   !  given that a neighbour cell owns the EB surface, a collision between!
   !  the particle and the EDGE of the EB facet might occur. This         !
   !  function returns the coordinates of the closest point on the edge of!
   !  an EB facet. This function does not check of collisions.            !
   !                                                                      !
   !  Comments: Position and normal vectors are represented as            !
   !  one-dimensional arrays of type real(c_real) and of dimension(3)     !
   !  (i.e. indices range from 1..3). Cells are enumerated using arrays of!
   !  type integer and of dimension(3).                                   !
   !----------------------------------------------------------------------!
    pure function facets_nearest_pt ( ind_pt, ind_loop, r_vec, eb_normal, eb_p0)
        implicit none

        real(c_real), dimension(3)             :: facets_nearest_pt
        integer,      dimension(3), intent(in) :: ind_pt, ind_loop
        real(c_real), dimension(3), intent(in) :: r_vec, eb_normal, eb_p0

        integer,      dimension(3) :: ind_facets
        integer                    :: n_facets, i_facet, tmp_facet, ind_cell, ind_nb
        real(c_real), dimension(3) :: c_vec, c_vec_tmp, rc_vec
        real(c_real), dimension(3) :: facet_normal, facet_p0, edge_p0, edge_v
        real(c_real)               :: min_dist, min_dist_tmp, eb_h, facet_h

        ! variables keeping track of coordinates on EB edges
        ! lambda_tmp: current lambda-value being used in testing for edge collions
        ! lambda: minimum (closets to bcentre) lambda value satisfying potential collision
        real(c_real) :: f_c, lambda_tmp, lambda_max, lambda_min


        ! Enumerate the possible EB facet edges invovlved.
        n_facets = 0

        if ( .not. (ind_pt(1) .eq. ind_loop(1)) ) then
            n_facets = n_facets + 1
            ind_facets(n_facets) = 1
        end if

        if ( .not. (ind_pt(2) .eq. ind_loop(2)) ) then
            n_facets = n_facets + 1
            ind_facets(n_facets) = 2
        end if

        if ( .not. (ind_pt(3) .eq. ind_loop(3)) ) then
            n_facets = n_facets + 1
            ind_facets(n_facets) = 3
        end if

        ! scalar characterizing EB facet position
        eb_h = dot_3d_real(eb_normal, eb_p0)

        ! itterate over EB facet edges and find whichever has the closest nearest point
        min_dist = huge(min_dist)
        do i_facet = 1, n_facets
            tmp_facet = ind_facets(i_facet)

            ! determine the normal of the cell's facet (cube faces)
            facet_normal = (/ 0., 0., 0. /)
            facet_normal(tmp_facet) = 1.  ! whether facing inwards or outwards is not important here

            ind_cell = ind_loop(tmp_facet)
            ind_nb = ind_pt(tmp_facet)

            ! determine position of the cell's facet
            if (ind_cell .lt. ind_nb) then
                f_c = ( dble(ind_cell) + 1.0 ) * dx(tmp_facet)
            else ! if (ind_cell .gt. ind_nb) then
                f_c = dble(ind_cell) * dx(tmp_facet)
            end if

            facet_p0 = (/                             &
                ( dble(ind_loop(1)) + 0.5 ) * dx(1) , &
                ( dble(ind_loop(2)) + 0.5 ) * dx(2) , &
                ( dble(ind_loop(3)) + 0.5 ) * dx(3)   &
            /)
            facet_p0(tmp_facet) = f_c

            ! scalar characterizing cell facet position
            facet_h = dot_3d_real(facet_normal, facet_p0)

            ! compute EB facet edge by finding the intercept between EB surface (first plane)
            ! and the cell's facet (second plane)
            call calc_facet_edge (edge_p0, edge_v, eb_h, facet_h, eb_normal, facet_normal)
            ! this solution is a line representing the closest EB edge, now compute the point
            ! on the line which minimizes the distance to the particle
            call lines_nearest_pt (lambda_tmp, c_vec_tmp, edge_p0, edge_v, r_vec)

            ! IMPORTANT: this point might be outside the cell
            !  -> in that case, it will be one of the cell's corners
            
            ! but don't use the PT_IN_BOX function, as it can yield false positives / negatives
            !  -> but if you do want to use it, include test [1] below to avoid rounding errors
            !if (.not. pt_in_box(c_vec_tmp, ind_loop, tmp_facet)) then
            
            ! if closest point is outside cell, determine the furthest we can go along the
            ! EB edge line whilst staying within the cell.
            call lambda_bounds(lambda_min, lambda_max, ind_loop, edge_p0, edge_v)
            
            ! [1]: why this test? What if (due to rounding 
            if (lambda_tmp .lt. lambda_min) then
                lambda_tmp = lambda_min
            elseif ( lambda_tmp .gt. lambda_max) then  ! [1] (see above)
                lambda_tmp = lambda_max
            end if
            c_vec_tmp(:) = edge_p0(:) + lambda_tmp*edge_v(:)
            
            !end if

            ! determine new distance to particle
            rc_vec(:) = c_vec_tmp(:) - r_vec(:)
            min_dist_tmp = dot_3d_real(rc_vec, rc_vec)

            ! minimize distance
            if (min_dist_tmp .lt. min_dist) then
                min_dist = min_dist_tmp
                c_vec(:) = c_vec_tmp(:)
            end if
        end do

        facets_nearest_pt(:) = c_vec(:)

    end function facets_nearest_pt

end subroutine fill_levelset_eb

subroutine validate_levelset(lo,    hi,   n_pad, &
                             impf,  imlo, imhi,  &
                             valid, vlo,  vhi,   &
                             phi,   phlo, phhi)  &
           bind(C, name="validate_levelset")

    implicit none

    integer,      dimension(3), intent(in   ) :: lo, hi, imlo, imhi, vlo, vhi, phlo, phhi
    integer,                    intent(in   ) :: n_pad
    real(c_real),               intent(in   ) :: impf  ( imlo(1):imhi(1), imlo(2):imhi(2), imlo(3):imhi(3) )
    integer,                    intent(in   ) :: valid (  vlo(1):vhi(1),   vlo(2):vhi(2),   vlo(3):vhi(3)  )
    real(c_real),               intent(inout) :: phi   ( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )

    integer      :: ii, jj, kk, klo, khi, jlo, jhi, ilo, ihi
    real(c_real) :: levelset_node

    klo = lo(3)! - n_pad
    khi = hi(3)! + n_pad
    jlo = lo(2)! - n_pad
    jhi = hi(2)! + n_pad
    ilo = lo(1)! - n_pad
    ihi = hi(1)! + n_pad

    do kk = klo, khi
        do jj = jlo, jhi
            do ii = ilo, ihi
                !write(*,*) "validating: ", ii, jj, kk, "v:", valid(ii, jj, kk), "if:", impf(ii, jj, kk)
                if ( valid(ii, jj, kk)  == 0 ) then
                    levelset_node = dabs( phi(ii, jj, kk) )
                    if ( impf(ii, jj, kk) <= 0 ) then
                        phi(ii, jj, kk) = levelset_node
                    else
                        phi(ii, jj, kk) = -levelset_node
                    end if
                end if
            end do
        end do
    end do
 
end subroutine validate_levelset


subroutine update_levelset(lo,    hi,           &
                           ls_in, lslo, lshi,   &
                           valid, vlo,  vhi,    &
                           phi,   phlo, phhi,   &
                           dx,    n_refine    ) &
           bind(C, name="update_levelset")

    implicit none

    integer,      dimension(3), intent(in   ) :: lo, hi, lslo, lshi, vlo, vhi, phlo, phhi
    integer,                    intent(in   ) :: n_refine
    real(c_real),               intent(in   ) :: ls_in (lslo(1):lshi(1),lslo(2):lshi(2),lslo(3):lshi(3))
    integer,                    intent(  out) :: valid ( vlo(1):vhi(1),  vlo(2):vhi(2),  vlo(3):vhi(3) )
    real(c_real),               intent(  out) :: phi   (phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))
    real(c_real), dimension(3), intent(in   ) :: dx

    real(c_real), dimension(3) :: pos_node
    real(c_real)               :: levelset_node

    integer :: i, j, k, ii, jj, kk, klo, khi, jlo, jhi, ilo, ihi
    logical :: valid_cell

    klo = lo(3)! - n_refine
    khi = hi(3)! + n_refine
    jlo = lo(2)! - n_refine
    jhi = hi(2)! + n_refine
    ilo = lo(1)! - n_refine
    ihi = hi(1)! + n_refine

    do kk = klo, khi
        do jj = jlo, jhi
            do ii = ilo, ihi
                levelset_node = ls_in(ii, jj, kk)

                !if ( dabs(levelset_node) < 5e-5 ) then
                !    write(*,*) ii, jj, kk, "pt=", (/ii, jj, kk/)*dx(:)/n_refine, "levelset=", levelset_node
                !end if

                !write(*,*) ii, jj, kk, phi(ii, jj, kk), levelset_node, sqrt(   &
                !    dot_product ( (/ii, jj/)*dx(1:2)/n_refine - (/0.0016, 0.0016/) , &
                !                  (/ii, jj/)*dx(1:2)/n_refine - (/0.0016, 0.0016/)  ) &
                !) + levelset_node

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

end subroutine update_levelset

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

subroutine eb_as_list(lo,       hi,   c_facets,  &
                      flag,     flo,  fhi,       &
                      norm,     nlo,  nhi,       &
                      bcent,    blo,  bhi,       &
                      list_out, lsize,           &
                      dx,       n_refine       ) &
           bind(C, name="eb_as_list")

    implicit none

    integer,                        intent(in   ) :: n_refine, lsize
    integer, dimension(3),          intent(in   ) :: lo, hi, flo, fhi, nlo, nhi, blo, bhi
    real(c_real),                   intent(in   ) :: dx(3)
    integer,                        intent(in   ) :: flag  ( flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3) )
    real(c_real),                   intent(in   ) :: norm  ( nlo(1):nhi(1), nlo(2):nhi(2), nlo(3):nhi(3), 3)
    real(c_real),                   intent(in   ) :: bcent ( blo(1):bhi(1), blo(2):bhi(2), blo(3):bhi(3), 3)
    real(c_real), dimension(lsize), intent(  out) :: list_out
    integer,                        intent(inout) :: c_facets

    integer                    :: i, j, k, i_facet
    real(c_real), dimension(3) :: eb_cent

    i_facet = 6 * c_facets + 1

    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                if ( is_single_valued_cell( flag(i, j, k) ) ) then

                    c_facets = c_facets + 1

                    eb_cent(:) = ( bcent(i, j, k, :)                           &
                                   + (/ dble(i), dble(j), dble(k) /)           &
                                   + (/ 0.5d0, 0.5d0, 0.5d0 /) ) * dx(:)/n_refine

                    !write(*,*) "generating eb_cent at: ", eb_cent(:), sqrt( &
                    !                dot_product ( ( eb_cent(1:2) - (/0.0016, 0.0016/) ), &
                    !                              ( eb_cent(1:2) - (/0.0016, 0.0016/) ) )&
                    !               )

                     
                    list_out( i_facet     : i_facet + 2) = eb_cent(:)
                    list_out( i_facet + 3 : i_facet + 5) = norm(i, j, k, :)

                    i_facet = i_facet + 6

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
