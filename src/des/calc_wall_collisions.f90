module wall_collisions
    use amrex_fort_module, only : c_real => amrex_real
    use iso_c_binding    , only: c_int

    implicit none
contains

  subroutine calc_wall_collisions ( particles, np, nrp, tow, fc, dtsolid, &
       flag, fglo, fghi, normal, nlo, nhi, bcent, blo, bhi, &
       dx) &
      bind(C, name="calc_wall_collisions")

    use amrex_fort_module, only : c_real => amrex_real
    use iso_c_binding    , only: c_int

    use discretelement, only: des_coll_model_enum
    use discretelement, only: des_etat_wall, des_etan_wall, hert_kwn, hert_kwt, hertzian
    use discretelement, only: des_crossprdct
    use discretelement, only: kn_w, kt_w, mew_w!, dtsolid

    use param,         only: small_number, zero, one
    use particle_mod,  only: particle_t
    use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, is_single_valued_cell, &
                                        get_neighbor_cells

    use eb_geometry, only: dot_3d_real, facets_nearest_pt
    
    implicit none

    integer, intent(in) :: np, nrp

    type(particle_t), intent(in   ), target :: particles(np)
    real(c_real)  ,   intent(inout)         :: tow(np,3), fc(np,3)
    real(c_real)  ,   intent(in   )         :: dtsolid

    integer, dimension(3), intent(in) :: fglo,fghi
    integer, dimension(3), intent(in) :: blo,bhi
    integer, dimension(3), intent(in) :: nlo,nhi

    integer,      intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    real(c_real), intent(in) :: bcent(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)
    real(c_real), intent(in) :: normal(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),3)
    real(c_real), intent(in) :: dx(3)

    type(particle_t), pointer :: p

    real(c_real) :: lx, ly, lz
    real(c_real) :: sqrt_overlap

    real(c_real) :: normul(3)
    !real(c_real) :: wca_overlap_factor, wca_strength, wca_radius, wca_inv_r, wca_offset, f_wca
    !real(c_real) :: v_normal, wca_dist

    ! facet barycenter (bcent) in global coordinates
    real(c_real), dimension(3) :: eb_cent
    ! position on facet plane, closest to particle
    real(c_real), dimension(3) :: eb_min_pt
    ! vector-indices of loop (over neighbour cells) and collision pt
    integer, dimension(3)      :: vindex_loop, vindex_pt
    logical                    :: nb_normal_valid;

    integer :: ll, ii, jj, kk, i, j, k, i_pt, j_pt, k_pt

    real(c_real) :: fudge
    real(c_real) :: overlap_n
    real(c_real) :: inv_dx(3)
    integer      :: ilo,ihi,jlo,jhi,klo,khi

    real(c_real) :: v_rel_trans_norm

    ! local normal and tangential forces
    real(c_real) :: vrel_t(3), cur_distmod, cur_plane_distmod
    real(c_real) :: ft(3), fn(3), overlap_t(3)

    ! worst-case: overlap with 27 neighbours
    integer                        :: n_collisions, i_collision
    real(c_real), dimension(27)    :: distmod, plane_distmod
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
    real(c_real) :: distmod_temp, distmod_edge

    ! inverse cell size: used to convert positions to cell indices
    !   -> dx is a vector, fortran is amazing!
    inv_dx = 1.0d0 / dx

    ! fudge factor: a little less than 1
    fudge = one - 1.d-8

    ! from WCA potential: r_min = wca_overlap_factor * wca_radius
    !wca_overlap_factor = 2.**(1./6.)

    ! iterate over particles
    do ll = 1, nrp
       ! get current particle
       p => particles(ll)

       ! particle centre position
       xp = p%pos(1)
       yp = p%pos(2)
       zp = p%pos(3)

       ! particle radia
       rp = p%radius

       ! WCA-interaction kicks in _after_ MFIX-interaction
       !  => wca_radius < particle radius
       !  => most MFIX collisions don't see WCA interaction
       !wca_radius = rp
       !wca_strength = p%mass * 10000;

       ! particle position (in units of cells)
       lx = xp*inv_dx(1)
       ly = yp*inv_dx(2)
       lz = zp*inv_dx(3)

       ! cell indices for position corresponding to (lx, ly, lz)
       i = floor(lx)
       j = floor(ly)
       k = floor(lz)

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
                if ( ( is_single_valued_cell(flag(ii, jj, kk))) ) then
!                    call bl_proffortfuncstart_int(1)
                    ! convert bcent to global coordinate system centered at plo
                    ! bcentx = bcent(ii, jj, kk, 1) * dx(1) + (dble(ii) + 0.5d0) * dx(1)
                    eb_cent(:) = ( bcent(ii, jj, kk, :) + (/dble(ii), dble(jj), dble(kk)/) + (/.5d0, .5d0, .5d0/) )*dx(:)

                    ! Distance to closest point on EB
                    ! distmod_temp = dabs( dot_3d_real( p%pos(:) - eb_cent(:), -normal(ii, jj, kk, :) ) )
                    distmod_temp =  dot_3d_real( p%pos(:) - eb_cent(:), -normal(ii, jj, kk, :) )


                    if (dabs(distmod_temp) .lt. rp) then
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
                            distmod(n_collisions) = dabs(distmod_temp)
                            !if ( dabs( distmod_temp ) .lt. 1e-8) then
                            !    write(*,*) distmod_temp
                            !end if

                            plane_distmod(n_collisions) = distmod_temp
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
                                 nb_normal_valid = .false.
                             else
                                 nb_normal_valid = .true.
                            end if

                            if (nb_normal_valid) then
                                if ( all( abs( normal(ii, jj, kk, :) - normal(i_pt, j_pt, k_pt, :) ) < 1.d-8 ) ) then
                                    ! EB-facet has same normal => skip
!                                    call bl_proffortfuncstop_int(1)
                                    cycle
                                end if

                            end if

                            ! Particle-EB-plain intercept is outside cell
                            !  -> determine nearest point on EB boundary to particle (candidate for collision)
                            !  -> c_vec will be on the "edge" of the EB boundary
                            c_vec = facets_nearest_pt (vindex_pt, vindex_loop, p%pos, normal(ii, jj, kk, :), eb_cent, dx)

                            ! check if there is still overlap
                            delta_c_vec(:) = p%pos(:) - c_vec(:)
                            distmod_edge = dot_3d_real(delta_c_vec, delta_c_vec)

                            if (distmod_edge .lt. rp * rp) then
                                ! The nearest point on EB edge overlaps with particle
                                !  => Collison still occurs! => Register collision with EB edge
                                distmod_edge = sqrt(distmod_edge)
                                n_collisions = n_collisions + 1
                                distmod(n_collisions) = distmod_edge
                                plane_distmod(n_collisions) = distmod_temp

                                ! assign average normal to edge collision:
                                if (nb_normal_valid) then
                                    collision_norms(:, n_collisions) = normal(ii, jj, kk, :) + normal(i_pt, j_pt, k_pt, :)
                                    len2_collision_norm = dot_3d_real( collision_norms(:, n_collisions), &
                                                                       collision_norms(:, n_collisions))
                                    collision_norms(:, n_collisions) = collision_norms(:, n_collisions)  &
                                                                     / sqrt(len2_collision_norm)
                                else
                                    collision_norms(:, n_collisions) = normal(ii, jj, kk, :)
                                endif
                            end if ! distmod < rp*rp

                        end if ! else of "head-on" test
                    end if ! if distmod < rp
!                    call bl_proffortfuncstop_int(1)

                end if ! if cut cell
             end do  ! ii
          end do  ! jj
       end do  ! kk

       ! ******************************************************************** !
       !                                                                      !
       ! APPLY Enumerated Collisons to particles                              !
       !                                                                      !
       ! ******************************************************************** !

       do i_collision = 1, n_collisions
!          call bl_proffortfuncstart_int(2)

          cur_distmod = distmod(i_collision)
          cur_plane_distmod = plane_distmod(i_collision)
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

          ! Add WCA force (to mittigate wall-penetration)
          !f_wca = 0;
          !wca_dist = cur_plane_distmod + rp/2 ! * wca_overlap_factor
          !
          !if ( wca_dist <= wca_radius * wca_overlap_factor ) then
          !   wca_inv_r = wca_radius / dabs(wca_dist)
          !   f_wca = 4. * wca_strength * (       &
          !                 12. * wca_inv_r**11   &
          !               - 6.  * wca_inv_r**5    )
          !
          !   !v_normal = dabs( dot_3d_real( p%vel(:) , -normul(:) ) )
          !   !p%vel(:) = p%vel(:) + normul(:) * v_normal
          !    write(*,*) "f_wca = ", f_wca, "wca_radius = ", wca_radius, "distmod = ", cur_distmod, "plane_distmod = ", cur_plane_distmod
          !    write(*,*) "normul = ", normul(:)
          !    write(*,*) "f old = ", fn
          !    write(*,*) "f new = ", -f_wca*normul(:)
          !    write(*,*) " "
          !end if
          !
          !fn(:) = fn(:) - f_wca * normul(:)


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

!          call bl_proffortfuncstop_int(2)
       !end if ! if test on (d < radius)
       end do
    end do ! loop over particles
end subroutine calc_wall_collisions

subroutine calc_wall_collisions_ls(particles, np,   nrp,     &
                                   tow,       fc,   dtsolid, &
                                   valid,     vlo,  vhi,     &
                                   phi,       phlo, phhi,    &
                                   dx,        n_refine     ) &
           bind(c, name="calc_wall_collisions_ls")

    use particle_mod,   only: particle_t
    use discretelement, only: des_coll_model_enum
    use discretelement, only: des_etat_wall, des_etan_wall, hert_kwn, hert_kwt, hertzian
    use discretelement, only: des_crossprdct
    use discretelement, only: kn_w, kt_w, mew_w

    use eb_levelset, only: interp_levelset, normal_levelset

    implicit none

    ! ** input varaibles

    type(particle_t), intent(in   ), target :: particles(np)
    integer,          intent(in   )         :: np, nrp, vlo(3), vhi(3), phlo(3), phhi(3), n_refine

    real(c_real),     intent(inout)         :: tow(np,3), fc(np,3)
    real(c_real),     intent(in   )         :: dtsolid, dx(3)
    integer,          intent(in   )         :: valid( vlo(1): vhi(1),  vlo(2): vhi(2),  vlo(3): vhi(3) )
    real(c_real),     intent(in   )         :: phi(  phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )

    ! ** declare local varaibles:
    integer :: phasell

    ! p: pointer to current particle (in particles) array
    type(particle_t), pointer  :: p

    ! inv_dx: inverse cell size
    ! pos   : current particle's position
    ! plo   : origin of paritcle coordinates
    real(c_real), dimension(3) :: inv_dx, pos, plo, ft, fn, overlap_t, vrel_t, normal, tangent

    ! xp, yp, zp: current particle's x, y, z position
    ! rp:         current particle's radius
    ! lx, ly, lz: position in units of of cell length
    ! ls_value:   value of the level-set function
    ! overlap:    particle-wall overlap
    real(c_real)               :: xp, yp, zp, rp, lx, ly, lz, ls_value, overlap_n

    ! ll:      do-loop counter itterating over particles
    ! i, j, k: indices of cell containing current particle
    integer                    :: ll, i, j, k

    ! ls_valid: indicates if particle is near wall (and level-set needs to be tested)
    logical                    :: ls_valid

    real(c_real) :: v_rel_trans_norm, sqrt_overlap, kn_des_w, kt_des_w, etan_des_w, etat_des_w
    real(c_real) :: mag_overlap_t, fnmd

    ! inverse cell size: used to convert positions to cell indices
    inv_dx(:) = 1.0d0 / dx(:)

    plo = (/ 0., 0., 0. /)

    do ll = 1, nrp
        ! get current particle
        p => particles(ll)

        ! particle position
        xp = p%pos(1)
        yp = p%pos(2)
        zp = p%pos(3)
        ! in units of cell length
        lx = xp * inv_dx(1)
        ly = yp * inv_dx(2)
        lz = zp * inv_dx(3)
        ! get particle cell index
        i = floor(lx)
        j = floor(ly)
        k = floor(lz)

        ! checking valid seems to take more time than just using the level-set
        !if (valid(i*n_refine, j*n_refine, k*n_refine) .eq. 1) then
            ! particle radius
            rp  = p%radius

            ! compute particle/wall overlap
            pos = (/ xp, yp, zp /)

            ! interpolates levelset from nodal phi to position pos
            call interp_levelset(pos, plo, n_refine, phi, phlo, phhi, dx, ls_value);
            overlap_n = rp - ls_value

            if (ls_value .lt. rp) then

                call normal_levelset(pos, plo, n_refine, phi, phlo, phhi, dx, normal)
                normal(:) = -normal(:)

                ! *****************************************************************************
                ! calculate the translational relative velocity

                call cfrelvel_wall(ll, v_rel_trans_norm, vrel_t, normal, ls_value, particles)

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

                ! calculate the spring model parameters.
                phasell = particles(ll) % phase

                ! hertz vs linear spring-dashpot contact model
                if ( des_coll_model_enum == hertzian ) then
                    sqrt_overlap = sqrt(overlap_n)
                    kn_des_w     = hert_kwn(phasell)*sqrt_overlap
                    kt_des_w     = hert_kwt(phasell)*sqrt_overlap
                    sqrt_overlap = sqrt(sqrt_overlap)
                    etan_des_w   = des_etan_wall(phasell)*sqrt_overlap
                    etat_des_w   = des_etat_wall(phasell)*sqrt_overlap
                else
                    kn_des_w     = kn_w
                    kt_des_w     = kt_w
                    etan_des_w   = des_etan_wall(phasell)
                    etat_des_w   = des_etat_wall(phasell)
                end if

                ! calculate the normal contact force
                fn(:) = -(kn_des_w * overlap_n  + etan_des_w * v_rel_trans_norm) * normal(:)

                ! calculate the tangential displacement.
                overlap_t(:) = dtsolid*vrel_t(:)
                mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t))

                ! check for coulombs friction law and limit the maximum value of the
                ! tangential force on a particle in contact with a wall.
                if ( mag_overlap_t > 0.0 ) then

                    ! max force before the on set of frictional slip.
                    fnmd = mew_w*sqrt(dot_product(fn,fn))

                    ! direction of tangential force.
                    tangent = overlap_t/mag_overlap_t

                    ft = -fnmd * tangent

                else
                    ft = 0.0
                end if

                ! add the collision force to the total forces acting on the particle.
                fc(ll,:) = fc(ll,:) + fn(:) + ft(:)

                ! add the torque force to the total torque acting on the particle.
                tow(ll,:) = tow(ll,:) + ls_value*des_crossprdct(normal(:),ft)

            end if

        !end if

    end do


end subroutine calc_wall_collisions_ls


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


end module wall_collisions
