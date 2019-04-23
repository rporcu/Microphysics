module wall_collisions
    use amrex_fort_module, only : rt => amrex_real
    use iso_c_binding    , only: c_int

    implicit none
contains

  subroutine ls_has_walls(has_wall, phi, phlo, phhi, tol) bind(c, name="ls_has_walls")

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




    subroutine calc_wall_collisions(particles, np,   nrp,     &
                                    tow,       fc,   dtsolid, &
                                    valid,     vlo,  vhi,     &
                                    phi,       phlo, phhi,    &
                                    dx,        n_refine     ) &
               bind(c, name="calc_wall_collisions")

        use particle_mod,   only: particle_t
        use discretelement, only: des_etat_wall, des_etan_wall
        use discretelement, only: des_crossprdct
        use discretelement, only: kn_w, kt_w, mew_w

        use amrex_eb_levelset_module, only: amrex_eb_interp_levelset, amrex_eb_normal_levelset

        implicit none

        ! ** input varaibles

        type(particle_t), intent(in   ), target :: particles(np)
        integer,          intent(in   )         :: np, nrp, vlo(3), vhi(3), phlo(3), phhi(3), n_refine

        real(rt),     intent(inout)         :: tow(np,3), fc(np,3)
        real(rt),     intent(in   )         :: dtsolid, dx(3)
        integer,          intent(in   )         :: valid( vlo(1): vhi(1),  vlo(2): vhi(2),  vlo(3): vhi(3) )
        real(rt),     intent(in   )         :: phi(  phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )

        ! ** declare local varaibles:
        integer :: phasell

        ! p: pointer to current particle (in particles) array
        type(particle_t), pointer  :: p

        ! inv_dx: inverse cell size
        ! pos   : current particle's position
        ! plo   : origin of paritcle coordinates
        real(rt), dimension(3) :: inv_dx, pos, plo, ft, fn, overlap_t, vrel_t, normal, tangent

        ! xp, yp, zp: current particle's x, y, z position
        ! rp:         current particle's radius
        ! lx, ly, lz: position in units of of cell length
        ! ls_value:   value of the level-set function
        ! overlap:    particle-wall overlap
        real(rt)               :: rp, ls_value, overlap_n

        !---------------------------------------------------------------------
        ! additional relative translational motion due to rotation
        real(rt) :: v_rot(3)
        ! total relative velocity at contact point
        real(rt) :: vreltrans(3)
        !---------------------------------------------------------------------

        ! ll:      do-loop counter itterating over particles
        ! i, j, k: indices of cell containing current particle
        integer                    :: ll

        ! ls_valid: indicates if particle is near wall (and level-set needs to be tested)
        !logical                    :: ls_valid

        real(rt) :: v_rel_trans_norm, kn_des_w, kt_des_w, etan_des_w, etat_des_w
        real(rt) :: mag_overlap_t, fnmd

        !---------------------------------------------------------------------

        ! inverse cell size: used to convert positions to cell indices
        inv_dx(:) = 1.0d0 / dx(:)

        plo = (/ 0., 0., 0. /)

        do ll = 1, nrp
            ! get current particle
            p => particles(ll)

            ! checking valid seems to take more time than just using the level-set
            !if (valid(i*n_refine, j*n_refine, k*n_refine) .eq. 1) then
                ! particle radius
                rp  = p%radius

                ! compute particle/wall overlap
                !pos = (/ xp, yp, zp /)
                pos = p%pos

                ! interpolates levelset from nodal phi to position pos
                call amrex_eb_interp_levelset(pos, plo, n_refine, phi, phlo, phhi, dx, ls_value);
                overlap_n = rp - ls_value

                if (ls_value .lt. rp) then

                    call amrex_eb_normal_levelset(pos, plo, n_refine, phi, phlo, phhi, dx, normal)
                    normal(:) = -normal(:)

                    ! *****************************************************************************
                    !  Calculate the normal and tangential components of the      !
                    !  relative velocity between a particle and wall contact.              !
                    !  Note: Only the magnitude of the normal component is returned    !
                    !  whereas the full tangential vector is returned.                     !
                    ! call       cfrelvel_wall (ll, v_rel_trans_norm, vrel_t, normal, ls_value, particles)
                    ! subroutine cfrelvel_wall (ll, vrn,              vrt,    norm,   dist,     particles )
                    ! *****************************************************************************
                    ! Total relative velocity + rotational contribution
                    v_rot = ls_value * particles(ll) % omega
                    vreltrans(:) =  particles(ll) % vel  + des_crossprdct(v_rot, normal)

                    ! Magnitude of normal component of relative velocity (scalar)
                    v_rel_trans_norm = dot_product(vreltrans,normal)

                    ! Total relative translational slip velocity at the contact point
                    ! equation (8) in tsuji et al. 1992
                    vrel_t(:) =  vreltrans(:) - v_rel_trans_norm*normal(:)
                    ! *****************************************************************************

                    ! calculate the spring model parameters.
                    phasell = particles(ll) % phase

                    kn_des_w     = kn_w
                    kt_des_w     = kt_w
                    etan_des_w   = des_etan_wall(phasell)
                    etat_des_w   = des_etat_wall(phasell)

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

    end subroutine calc_wall_collisions

end module wall_collisions
