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
    real(c_real) :: sqrt_overlap

    real(c_real) :: normul(3)
    
    ! facet barycenter (bcent) in global coordinates
    real(c_real), dimension(3) :: eb_cent
    ! position on facet plane, closest to particle
    real(c_real), dimension(3) :: eb_min_pt
    ! vector-indices of loop (over neighbour cells) and collision pt
    integer, dimension(3)      :: vindex_loop, vindex_pt

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
    real(c_real), dimension(27)    :: distmod 
    real(c_real), dimension(3, 27) :: collision_norms
    real(c_real), dimension(3)     :: c_vec, delta_c_vec
    real(c_real)                   :: len2_collision_norm

    integer :: PHASELL

    real(c_real) :: tangent(3)

    real(c_real) :: xp, yp, zp, rp

    ! local values used spring constants and damping coefficients
    real(c_real) :: etan_des_W, etat_des_W, kn_des_W, kt_des_W
    real(c_real) :: fnmd
    real(c_real) :: mag_overlap_t
    real(c_real) :: distmod_temp

    ! inverse cell size: used to convert positions to cell indices
    !   -> dx is a vector, fortran is amazing!
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
       ! -> one reason for this is the get accurate wall normals...
       call get_neighbor_cells(flag(i,j,k),nbr)
       
       ! ******************************************************************** !
       !                                                                      !
       ! BEGIN Enumerating Collisons with embedded boundary                   !
       !                                                                      !
       ! ******************************************************************** !

       n_collisions = 0

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

       do kk = klo, khi
          do jj = jlo, jhi
             do ii = ilo, ihi

                ! only consider cells that contain EB's
                if ( (nbr(ii-i, jj-j, kk-k) .eq. 1) .and. ( .not. is_regular_cell(flag(ii, jj, kk))) ) then
                    ! convert bcent to global coordinate system centered at plo
                    ! bcentx = bcent(ii, jj, kk, 1) * dx(1) + (dble(ii) + 0.5d0) * dx(1)
                    eb_cent(:) = ( bcent(ii, jj, kk, :) + (/dble(ii), dble(jj), dble(kk)/) + (/.5d0, .5d0, .5d0/) )*dx(:)

                    ! Distance to closest point on EB
                    distmod_temp = dabs( dot_3d_real( p%pos(:) - eb_cent(:), -normal(ii, jj, kk, :) ) )

                    if (distmod_temp .lt. rp) then
                        ! vector-index of triple-loop
                        vindex_loop = (/ii, jj, kk/)

                        ! Point where the normal from the particle to the plane intersects the plane
                        !  -> Note: slightly less than the full normal avoiding float precision errors
                        eb_min_pt(:) = p%pos(:) + normal(ii, jj, kk, :) * distmod_temp * fudge

                        ! Cell that point is in
                        vindex_pt = floor( eb_min_pt(:) * inv_dx(:) )

                        ! Check if point-of-collision is in current cell, that's great!
                        !  -> particle is interacting with "flat" EB's
                        if ( all( vindex_loop == vindex_pt ) ) then
                            ! register "head-on" collision:
                            n_collisions = n_collisions + 1
                            distmod(n_collisions) = distmod_temp
                            collision_norms(:, n_collisions) = normal(ii, jj, kk, :)
                        
                        ! only consider EB-edges if not colliding with EB-face already 
                        !  -> EB-facet has different normal => collide with EB-edge
                        else
                            i_pt = vindex_pt(1)
                            j_pt = vindex_pt(2) 
                            k_pt = vindex_pt(3)

                            ! Don't ask for a normal in a cell that doesn't have a wall
                            if ( is_regular_cell(flag(i_pt,j_pt,k_pt)) .or. &
                                 is_covered_cell(flag(i_pt,j_pt,k_pt)) ) then
                                cycle
                            end if 

                            if ( all( abs( normal(ii, jj, kk, :) - normal(i_pt, j_pt, k_pt, :) ) < 1.d-8 ) ) then
                                ! EB-facet has same normal => skip
                                cycle
                            end if 

                            ! Particle-EB-plain intercept is outside cell
                            !  -> determine nearest point on EB boundary to particle (candidate for collision)
                            !  -> c_vec will be on the "edge" of the EB boundary
                            c_vec = facets_nearest_pt (vindex_pt, vindex_loop, p%pos, normal(ii, jj, kk, :), eb_cent)

                            ! check if there is still overlap
                            delta_c_vec(:) = p%pos(:) - c_vec(:)
                            distmod_temp = dot_3d_real(delta_c_vec, delta_c_vec)

                            if (distmod_temp .lt. rp * rp) then
                                ! The nearest point on EB edge overlaps with particle
                                !  => Collison still occurs! => Register collision with EB edge
                                distmod_temp = sqrt(distmod_temp)
                                n_collisions = n_collisions + 1
                                distmod(n_collisions) = distmod_temp

                                ! assign average normal to edge collision:
                                collision_norms(:, n_collisions) = normal(ii, jj, kk, :) + normal(i_pt, j_pt, k_pt, :)
                                len2_collision_norm = dot_3d_real( collision_norms(:, n_collisions), &
                                                                   collision_norms(:, n_collisions))
                                collision_norms(:, n_collisions) = collision_norms(:, n_collisions)  &
                                                                   / sqrt(len2_collision_norm)
                            end if
                        end if
                    end if
                end if
             end do
          end do
       end do

       ! ******************************************************************** !
       !                                                                      !
       ! APPLY Enumerated Collisons to particles                              !
       !                                                                      !
       ! ******************************************************************** !

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
          if ( des_coll_model_enum == HERTZIAN ) then
             sqrt_overlap = sqrt(overlap_n)
             kn_des_w     = hert_kwn(phaseLL)*sqrt_overlap
             kt_des_w     = hert_kwt(phaseLL)*sqrt_overlap
             sqrt_overlap = sqrt(sqrt_overlap)
             etan_des_w   = des_etan_wall(phaseLL)*sqrt_overlap
             etat_des_w   = des_etat_wall(phaseLL)*sqrt_overlap
          else
             kn_des_w     = kn_w
             kt_des_w     = kt_w
             etan_des_w   = des_etan_wall(phaseLL)
             etat_des_w   = des_etat_wall(phaseLL)
          end if

          ! Calculate the normal contact force
          fn(:) = -(kn_des_w * overlap_n  + etan_des_w * v_rel_trans_norm) * normul(:)

          ! Calculate the tangential displacement.
          overlap_t(:) = dtsolid*vrel_t(:)
          mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t))

          ! Check for Coulombs friction law and limit the maximum value of the
          ! tangential force on a particle in contact with a wall.
          if ( mag_overlap_t > 0.0 ) then

             ! Max force before the on set of frictional slip.
             fnmd = MEW_W*sqrt(dot_product(fn,fn))

             ! Direction of tangential force.
             tangent = overlap_t/mag_overlap_t

             ft = -fnmd * tangent

          else
             ft = 0.0
          end if

          ! Add the collision force to the total forces acting on the particle.
          fc(LL,:) = fc(LL,:) + fn(:) + ft(:)

          ! Add the torque force to the total torque acting on the particle.
          tow(LL,:) = tow(LL,:) + cur_distmod*des_CROSSPRDCT(normul(:),ft)

       !end if ! if test on (d < radius)
       end do
    end do ! loop over particles

contains

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
    pure function pt_in_box (pt, id, id_ignore)
        implicit none

        logical                                :: pt_in_box
        real(c_real), dimension(3), intent(in) :: pt
        integer,      dimension(3), intent(in) :: id
        integer,                    intent(in) :: id_ignore

        real(c_real), dimension(3) :: box_max, box_min
        integer                    :: i
        
        ! Determine box boundaries
        box_min(:) = id(:) * dx(:)
        box_max(:) = (id(:) + 1.) * dx(:)

        pt_in_box = .true.

        ! Check each coordinate. Skip ignored coordinate.
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
        ! due to the huge(c...) defaults (above).
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
            facet_normal(tmp_facet) = 1.  ! wheter facing inwards or outwards is not important here

            ind_cell = ind_loop(tmp_facet)
            ind_nb = ind_pt(tmp_facet)

            ! determine position of the cell's facet
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
            if (.not. pt_in_box(c_vec_tmp, ind_loop, tmp_facet)) then
                ! if closest point is outside cell, determine the furthest we can go along the 
                ! EB edge line whilst staying within the cell.
                call lambda_bounds(lambda_min, lambda_max, ind_loop, edge_p0, edge_v)
                if (lambda_tmp .lt. lambda_min) then
                    lambda_tmp = lambda_min
                else
                    lambda_tmp = lambda_max
                end if

                c_vec_tmp(:) = edge_p0(:) + lambda_tmp*edge_v(:)
            end if
            
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


   !----------------------------------------------------------------------!
   !                                                                      !
   !  Subroutine: cfrelvel_wall                                           !
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
