  subroutine calc_wall_collisions ( particles, np, nrp, tow, fc, dtsolid, &
       flag, fglo, fghi, normal, nlo, nhi, bcent, blo, bhi, apx, axlo, axhi, apy, aylo, ayhi, &
       apz, azlo, azhi, dx) &
      bind(C, name="calc_wall_collisions")

    use amrex_fort_module, only : c_real => amrex_real
    use iso_c_binding    , only: c_int

    use discretelement, only: des_coll_model_enum
    use discretelement, only: des_etat_wall, des_etan_wall, hert_kwn, hert_kwt, hertzian
    use discretelement, only: des_crossprdct
    use discretelement, only: kn_w, kt_w, mew_w!, dtsolid

    use param        , only: small_number, zero, one
    use particle_mod,  only: particle_t
    use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, is_single_valued_cell, &
                                        get_neighbor_cells

    implicit none

    integer, intent(in) :: np, nrp

    type(particle_t), intent(in   ), target :: particles(np)
    real(c_real)  ,   intent(inout)         :: tow(np,3), fc(np,3)
    real(c_real)  ,   intent(in   )         :: dtsolid

    integer, dimension(3), intent(in) :: axlo,axhi
    integer, dimension(3), intent(in) :: aylo,ayhi
    integer, dimension(3), intent(in) :: azlo,azhi
    integer, dimension(3), intent(in) :: fglo,fghi
    integer, dimension(3), intent(in) :: blo,bhi
    integer, dimension(3), intent(in) :: nlo,nhi

    integer,      intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    real(c_real), intent(in) :: bcent(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)
    real(c_real), intent(in) :: normal(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),3)
    real(c_real), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(c_real), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(c_real), intent(in) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(c_real), intent(in) :: dx(3)

    type(particle_t), pointer :: p

    real(c_real) :: lx, ly, lz
    real(c_real) :: bcentx, bcenty, bcentz
    real(c_real) :: sqrt_overlap

    real(c_real) :: normul(3)

    integer :: ll, ii, jj, kk, i, j, k, i_pt, j_pt, k_pt

    real(c_real) :: fudge
    real(c_real) :: overlap_n
    real(c_real) :: inv_dx(3)
    integer      :: nbr(-1:1,-1:1,-1:1)
    integer      :: ilo,ihi,jlo,jhi,klo,khi

    real(c_real) :: v_rel_trans_norm

    ! local normal and tangential forces
    real(c_real) :: vrel_t(3), cur_distmod
    real(c_real) :: ft(3), fn(3), overlap_t(3)

    ! worst-case: overlap with 27 neighbours
    integer                        :: n_collisions, i_collision
    integer, dimension(3,3,3)      :: central_collision, edge_collision
    real(c_real), dimension(27)    :: distmod 
    real(c_real), dimension(3, 27) :: collision_norms, collision_points 
    real(c_real), dimension(3)     :: r_vec, eb_normal, eb_p0

    integer :: PHASELL

    real(c_real) :: tangent(3)

    real(c_real) :: xp, yp, zp, rp
    real(c_real) :: pt_x, pt_y, pt_z

    ! local values used spring constants and damping coefficients
    real(c_real) :: ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W
    real(c_real) :: fnmd
    real(c_real) :: mag_overlap_t
    real(c_real) :: distmod_temp

    ! inverse cell size: used to convert positions to cell indices
    ! dx is a vector, why does this work?
    inv_dx = 1.0d0 / dx

    ! fudge factor: a little less than 1
    fudge = one - 1.d-8

    ! itterate over particles
    do ll = 1, nrp
       ! get current particle
       p => particles(ll)

       ! particle centre position
       xp = p%pos(1)
       yp = p%pos(2)
       zp = p%pos(3)

       ! particle radia
       rp = p%radius

       ! particle position (in units of cells)
       lx = xp*inv_dx(1)
       ly = yp*inv_dx(2)
       lz = zp*inv_dx(3)

       ! cell indices for position corresponding to (lx, ly, lz)
       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

       ! ignore disconnected cells
       ! is this really necessary?
       ! -> one reason for this is the get accurate wall normals...
       call get_neighbor_cells(flag(i,j,k),nbr)

       ! 27-point stencil: immediate neighbours in 3-D cubic grid
       klo = k-1
       khi = k+1
       jlo = j-1
       jhi = j+1
       ilo = i-1
       ihi = i+1

       ! ingore stencil element that couldn't possibly be interacting with particle
       if ( (xp-i*dx(1)) .gt. rp) ilo = i
       if ( (yp-j*dx(2)) .gt. rp) jlo = j
       if ( (zp-k*dx(3)) .gt. rp) klo = k

       if ( ((i+1)*dx(1)-xp) .gt. rp) ihi = i
       if ( ((j+1)*dx(2)-yp) .gt. rp) jhi = j
       if ( ((k+1)*dx(3)-zp) .gt. rp) khi = k

       n_collisions = 0
       central_collision(:, :, :) = 0
       edge_collision(:, :, :) = 0
       collision_points(:, :) = 0.

       do kk = klo, khi
          do jj = jlo, jhi
             do ii = ilo, ihi
                ! only consider cells that contain EB's
                if ( (nbr(ii-i, jj-j, kk-k) .eq. 1) .and. ( .not. is_regular_cell(flag(ii, jj, kk))) ) then
                    ! convert bcent to global coordinate system centered at plo
                    bcentx = bcent(ii, jj, kk, 1) * dx(1) + (dble(ii) + 0.5d0) * dx(1)
                    bcenty = bcent(ii, jj, kk, 2) * dx(2) + (dble(jj) + 0.5d0) * dx(2)
                    bcentz = bcent(ii, jj, kk, 3) * dx(3) + (dble(kk) + 0.5d0) * dx(3)

                    ! Distance to boundary
                    distmod_temp = dabs( (xp - bcentx) * (-normal(ii, jj, kk, 1)) + &
                                           (yp - bcenty) * (-normal(ii, jj, kk, 2)) + &
                                           (zp - bcentz) * (-normal(ii, jj, kk, 3)) )

                    if (distmod_temp .lt. rp) then
                        ! Point where the normal from the particle to the plane intersects the plane
                        ! ->Note: slightly less than the full normal avoiding float precision errors
                        pt_x = xp + normal(ii, jj, kk, 1) * distmod_temp * fudge
                        pt_y = yp + normal(ii, jj, kk, 2) * distmod_temp * fudge
                        pt_z = zp + normal(ii, jj, kk, 3) * distmod_temp * fudge

                        ! Cell that point is in
                        i_pt = floor(pt_x * inv_dx(1))
                        j_pt = floor(pt_y * inv_dx(2))
                        k_pt = floor(pt_z * inv_dx(3))

                        ! check if point-of-collision is in current cell, that's great!
                        ! Particle is interacting with "flat" EB's
                        if ( (i_pt .eq. ii) .and. (j_pt .eq. jj) .and. (k_pt .eq. kk)) then
                            n_collisions = n_collisions + 1
                            central_collision(ii - i, jj - j, kk - k) = 1

                            ! distmod = distmod_temp
                            distmod(n_collisions) = distmod_temp
                            collision_points(:, n_collisions) = (/ pt_x, pt_y, pt_z /)

                            collision_norms(1, n_collisions) = normal(ii, jj, kk, 1)
                            collision_norms(2, n_collisions) = normal(ii, jj, kk, 2)
                            collision_norms(3, n_collisions) = normal(ii, jj, kk, 3)
                        else if (                                                                          &
                            ( abs( normal(ii, jj, kk, 1) - normal(i_pt, j_pt, k_pt, 1) ) .gt. 1.d-8 ) .or. &
                            ( abs( normal(ii, jj, kk, 2) - normal(i_pt, j_pt, k_pt, 2) ) .gt. 1.d-8 ) .or. &
                            ( abs( normal(ii, jj, kk, 3) - normal(i_pt, j_pt, k_pt, 3) ) .gt. 1.d-8 )      &
                            ) then
                            block
                                ! variables keepint track of edges/corner of neighbour EB's
                                integer, dimension(3) :: ind_loop, ind_pt
                                ! c_vec: closest point (to particle) on EB-facet edge 
                                real(c_real), dimension(3) :: c_vec
                                    
                                ! assing value to vectors (arrays)
                                ind_loop = (/ ii, jj, kk /)
                                ind_pt   = (/ i_pt, j_pt, k_pt /)
                                
                                r_vec = (/ xp, yp, zp /)
                                eb_normal(:) = normal(ii, jj, kk, :)
                                eb_p0 = (/ bcentx, bcenty, bcentz /)

                                ! p_vec is outside cell, determine with cell boundaries might be cut
                                c_vec = facets_nearest_pt ( ind_pt, ind_loop, r_vec, eb_normal, eb_p0 )
                                    
                                ! correct collision coordinates
                                pt_x = c_vec(1)
                                pt_y = c_vec(2)
                                pt_z = c_vec(3)

                                ! check if there is still overlap
                                distmod_temp = ((xp - pt_x) * ( xp - pt_x ) + & 
                                    (yp - pt_y) * ( yp - pt_y ) + &
                                    (zp - pt_z) * ( zp - pt_z ))

                                if (distmod_temp .lt. rp * rp) then
                                    if (.not. point_in_list(collision_points, c_vec, n_collisions)) then
                                        distmod_temp = sqrt(distmod_temp)
                                        
                                        n_collisions = n_collisions + 1
                            
                                        distmod(n_collisions) = distmod_temp
                                        collision_points(:, n_collisions) = (/ pt_x, pt_y, pt_z /)

                                        collision_norms(1, n_collisions) = -(xp - pt_x) / distmod_temp
                                        collision_norms(2, n_collisions) = -(yp - pt_y) / distmod_temp
                                        collision_norms(3, n_collisions) = -(zp - pt_z) / distmod_temp
                                    end if
                                end if
                            end block
                        end if
                    end if
                end if
             end do
          end do
       end do

       do i_collision = 1, n_collisions
          
          cur_distmod = distmod(i_collision)
          normul(:) = collision_norms(:, i_collision)

          ! Calculate the particle/wall overlap.
          overlap_n = rp - cur_distmod

          ! *****************************************************************************
          ! Calculate the translational relative velocity

          call cfrelvel_wall(ll, v_rel_trans_norm, vrel_t, normul, cur_distmod, particles)

          ! subroutine cfrelvel_wall (ll, vrn, vrt, norm, dist, particles )
          ! *****************************************************************************
          ! total relative velocity + rotational contribution
          ! v_rot = distmod * particles(ll) % omega
          ! vreltrans(:) =  particles(ll) % vel  + des_crossprdct(v_rot, normal)

          ! magnitude of normal component of relative velocity (scalar)
          ! v_rel_trans_norm = dot_product(vreltrans,normal)

          ! total relative translational slip velocity at the contact point
          ! equation (8) in tsuji et al. 1992
          ! vrel_t(:) =  vreltrans(:) - v_rel_trans_norm*normal(:)
          ! *****************************************************************************

          ! Calculate the spring model parameters.
          phaseLL = particles(ll) % phase

          ! Hertz vs linear spring-dashpot contact model
          if ( DES_COLL_MODEL_ENUM == HERTZIAN ) then
             sqrt_overlap = sqrt(overlap_n)
             KN_DES_W     = hert_kwn(phaseLL)*sqrt_overlap
             KT_DES_W     = hert_kwt(phaseLL)*sqrt_overlap
             sqrt_overlap = SQRT(sqrt_overlap)
             ETAN_DES_W   = DES_ETAN_WALL(phaseLL)*sqrt_overlap
             ETAT_DES_W   = DES_ETAT_WALL(phaseLL)*sqrt_overlap
          else
             KN_DES_W     = KN_W
             KT_DES_W     = KT_W
             ETAN_DES_W   = DES_ETAN_WALL(phaseLL)
             ETAT_DES_W   = DES_ETAT_WALL(phaseLL)
          end if

          ! Calculate the normal contact force
          FN(:) = -(KN_DES_W * overlap_n  + ETAN_DES_W * V_REL_TRANS_NORM) * normul(:)

          ! Calculate the tangential displacement.
          overlap_t(:) = dtsolid*vrel_t(:)
          mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t))

          ! Check for Coulombs friction law and limit the maximum value of the
          ! tangential force on a particle in contact with a wall.
          if ( mag_overlap_t > 0.0 ) then

             ! Max force before the on set of frictional slip.
             fnmd = MEW_W*sqrt(dot_product(FN,FN))

             ! Direction of tangential force.
             tangent = overlap_t/mag_overlap_t

             FT = -fnmd * tangent

          else
             FT = 0.0
          end if

          ! Add the collision force to the total forces acting on the particle.
          FC(LL,:) = FC(LL,:) + FN(:) + FT(:)
          ! print *,'Vertical dist / force ',distmod, fn(1), fc(ll,1)

          ! Add the torque force to the total torque acting on the particle.
          TOW(LL,:) = TOW(LL,:) + cur_distmod*DES_CROSSPRDCT(normul(:),FT)

       !end if ! if test on (d < radius)
       end do
    end do ! loop over particles

contains
    function point_in_list (list, pt, n_list)
        implicit none

        logical                                   :: point_in_list
        real(c_real), dimension(3,27), intent(in) :: list
        real(c_real), dimension(3),    intent(in) :: pt
        integer,                       intent(in) :: n_list

        real(c_real), dimension(3) :: current_pt
        integer                    :: i
        
        point_in_list = .false.
        do i = 1, n_list
            current_pt(:) = list(:, i)
            !write(*,*) current_pt
            if ( &
                    ( abs( pt(1) - current_pt(1) ) .lt. 1.d-8 ) .and. &
                    ( abs( pt(2) - current_pt(2) ) .lt. 1.d-8 ) .and. &
                    ( abs( pt(3) - current_pt(3) ) .lt. 1.d-8 )       &
                ) then
                point_in_list = .true.
                !write(*,*) "Point disqualified:"
                !write(*,*) pt
                exit
            end if
        end do


    end function point_in_list

    pure function dot_3d_real (v1, v2)
        implicit none
        
        real(c_real)                           :: dot_3d_real
        real(c_real), dimension(3), intent(in) :: v1, v2

        dot_3d_real = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

    end function dot_3d_real

    pure function cross_3d_real (v1, v2)
        implicit none

        real(c_real), dimension(3)             :: cross_3d_real
        real(c_real), dimension(3), intent(in) :: v1, v2

        cross_3d_real(1) = v1(2)*v2(3) - v1(3)*v2(2)
        cross_3d_real(2) = v1(3)*v2(1) - v1(1)*v2(3)
        cross_3d_real(3) = v1(1)*v2(2) - v1(2)*v2(1)

    end function cross_3d_real

    pure function pt_in_box (pt, id, id_ignore)
        implicit none

        logical                                :: pt_in_box
        real(c_real), dimension(3), intent(in) :: pt
        integer,      dimension(3), intent(in) :: id
        integer,                    intent(in) :: id_ignore

        real(c_real), dimension(3) :: box_max, box_min
        integer                    :: i

        box_min(:) = id(:) * dx(:)
        box_max(:) = (id(:) + 1.) * dx(:)

        pt_in_box = .true.

        do i = 1, 3
            if (.not. (i .eq. id_ignore) ) then
                if ( pt(i) .lt. box_min(i) ) then
                    pt_in_box = .false.
                    exit
                end if

                if ( pt(i) .gt. box_max(i) ) then
                    pt_in_box = .false.
                    exit
                end if
            end if
        end do

    end function pt_in_box

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

    pure subroutine swap_reals(a, b)
        implicit none
            
        real(c_real), intent(inout) :: a, b
        real(c_real)                :: bucket
        
        bucket = a
        a      = b
        b      = bucket

    end subroutine swap_reals

    pure subroutine lambda_bounds(lambda_min, lambda_max, id_cell, p0, v)
        implicit none 

        real(c_real),               intent(  out) :: lambda_min, lambda_max
        integer,      dimension(3), intent(in   ) :: id_cell
        real(c_real), dimension(3), intent(in   ) :: p0, v

        real(c_real) :: cx_lo, cy_lo, cz_lo, cx_hi, cy_hi, cz_hi!, c_bucket

        cx_lo = -huge(cx_lo)
        cy_lo = -huge(cy_lo)
        cz_lo = -huge(cz_lo)

        cx_hi = huge(cx_hi)
        cy_hi = huge(cy_hi)
        cz_hi = huge(cz_hi)

        if ( abs(v(1)) .gt. 1.d-8 ) then
            cx_lo = -( p0(1) - dble(id_cell(1)) * dx(1) ) / v(1)
            cx_hi = -( p0(1) - ( dble(id_cell(1)) + 1. ) * dx(1) ) / v(1)

            if ( v(1) .lt. 0. ) then
                call swap_reals(cx_lo, cx_hi)
            end if
        end if

        if ( abs(v(2)) .gt. 1.d-8 ) then
            cy_lo = -( p0(2) - dble(id_cell(2)) * dx(2) ) / v(2)
            cy_hi = -( p0(2) - ( dble(id_cell(2)) + 1. ) * dx(2) ) / v(2)

            if ( v(2) .lt. 0. ) then
                call swap_reals(cy_lo, cy_hi)
            end if
        end if

        if ( abs(v(3)) .gt. 1.d-8 )  then
            cz_lo = -( p0(3) - dble(id_cell(3)) * dx(3) ) / v(3)
            cz_hi = -( p0(3) - ( dble(id_cell(3)) + 1. ) * dx(3) ) / v(3)

            if ( v(3) .lt. 0. ) then 
                call swap_reals(cz_lo, cz_hi)
            end if
        endif

        lambda_min = max(cx_lo, cy_lo, cz_lo)
        lambda_max = min(cx_hi, cy_hi, cz_hi)
        
        ! write(*,*) cx_lo, cy_lo, cz_lo
        ! write(*,*) cx_hi, cy_hi, cz_hi
        ! write(*,*) lambda_min, lambda_max
        ! write(*,*) v

    end subroutine lambda_bounds

    pure function facets_nearest_pt ( ind_pt, ind_loop, r_vec, eb_normal, eb_p0)
        implicit none
        
        real(c_real), dimension(3)             :: facets_nearest_pt
        integer,      dimension(3), intent(in) :: ind_pt, ind_loop
        real(c_real), dimension(3), intent(in) :: r_vec, eb_normal, eb_p0
        
        integer,      dimension(3) :: ind_facets
        integer                    :: n_facets, i_facet, tmp_facet, cur_facet, ind_cell, ind_nb, i
        real(c_real), dimension(3) :: v_vec, c_vec, c_vec_tmp, rc_vec
        real(c_real), dimension(3) :: facet_normal, facet_p0, edge_p0, edge_v
        real(c_real)               :: min_dist, min_dist_tmp, eb_h, facet_h

        ! variables keeping track of coordinates on EB edges
        ! lambda_tmp: current lambda-value being used in testing for edge collions
        ! lambda: minimum (closets to bcentre) lambda value satisfying potential collision
        real(c_real) :: v_c, f_c, p_c, lambda, lambda_tmp, lambda_max, lambda_min

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

        eb_h = dot_3d_real(eb_normal, eb_p0)

        min_dist = huge(min_dist)
        do i_facet = 1, n_facets
            tmp_facet = ind_facets(i_facet)

            facet_normal = (/ 0., 0., 0. /)
            facet_normal(tmp_facet) = 1.  ! wheter facing inwards or outwards is not important here

            ind_cell = ind_loop(tmp_facet)
            ind_nb = ind_pt(tmp_facet)

            if (ind_cell .lt. ind_nb) then
                f_c = ( dble(ind_cell) + 1.0 ) * dx(tmp_facet)
            else ! if (ind_cell .gt. ind_nb) then
                f_c = dble(ind_cell) * dx(tmp_facet)
            end if
            
            facet_p0 = (/                            & 
                ( dble(ind_loop(1)) + 0.5 ) * dx(1), &
                ( dble(ind_loop(2)) + 0.5 ) * dx(2), &
                ( dble(ind_loop(3)) + 0.5 ) * dx(3)  &
            /)
            facet_p0(tmp_facet) = f_c
            facet_h = dot_3d_real(facet_normal, facet_p0)

            call calc_facet_edge (edge_p0, edge_v, eb_h, facet_h, eb_normal, facet_normal)
            call lines_nearest_pt (lambda_tmp, c_vec_tmp, edge_p0, edge_v, r_vec)
            
            if (.not. pt_in_box(c_vec_tmp, ind_loop, tmp_facet)) then
                call lambda_bounds(lambda_min, lambda_max, ind_loop, edge_p0, edge_v)
                if (lambda_tmp .lt. lambda_min) then
                    lambda_tmp = lambda_min
                else
                    lambda_tmp = lambda_max
                end if

                c_vec_tmp(:) = edge_p0(:) + lambda_tmp*edge_v(:)
            end if
            
            rc_vec(:) = c_vec_tmp(:) - r_vec(:)
            min_dist_tmp = dot_3d_real(rc_vec, rc_vec)

            if (min_dist_tmp .lt. min_dist) then
                min_dist = min_dist_tmp
                c_vec(:) = c_vec_tmp(:)
            end if
        end do

        facets_nearest_pt(:) = c_vec(:)

       end function facets_nearest_pt


   !----------------------------------------------------------------------!
   !                                                                      !
   !  Subroutine: CFRELVEL_WALL                                           !
   !                                                                      !
   !  Purpose: Calculate the normal and tangential components of the      !
   !  relative velocity between a particle and wall contact.              !
   !                                                                      !
   !  Comments: Only the magnitude of the normal component is returned    !
   !  whereas the full tangential vector is returned.                     !
   !----------------------------------------------------------------------!
   subroutine cfrelvel_wall (ll, vrn, vrt, norm, dist, particles )

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      ! function for calculating the cross prodcut
      use discretelement, only: des_crossprdct
      use particle_mod,   only: particle_t

      type(particle_t), intent(in) :: particles(:)

      ! dummy arguments:
      !---------------------------------------------------------------------//
      ! particle index.
      integer, intent(in) :: ll
      ! magnitude of the total relative translational velocity.
      real(c_real), intent(out):: vrn
      ! total relative translational velocity (vector).
      real(c_real), intent(out):: vrt(3)
      ! unit normal from particle center to closest point on wall
      real(c_real), intent(in) :: norm(3)
      ! distance between particle center and wall.
      real(c_real), intent(in) :: dist

      ! local variables
      !---------------------------------------------------------------------//
      ! additional relative translational motion due to rotation
      real(c_real) :: v_rot(3)
      ! total relative velocity at contact point
      real(c_real) :: vreltrans(3)

      ! total relative velocity + rotational contribution
      v_rot = dist * particles(ll) % omega
      vreltrans(:) =  particles(ll) % vel  + des_crossprdct(v_rot, norm)

      ! magnitude of normal component of relative velocity (scalar)
      vrn = dot_product(vreltrans,norm)

      ! total relative translational slip velocity at the contact point
      ! equation (8) in tsuji et al. 1992
      vrt(:) =  vreltrans(:) - vrn*norm(:)

   end subroutine cfrelvel_wall

  end subroutine calc_wall_collisions
