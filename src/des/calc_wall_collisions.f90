
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
    use error_manager,  only: err_msg, flush_err_msg, init_err_msg

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
    real(c_real), intent(in) :: bcent  (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)
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
    real(c_real) :: vrel_t(3), distmod
    real(c_real) :: ft(3), fn(3), overlap_t(3)

    integer :: PHASELL

    real(c_real) :: tangent(3)

    real(c_real) :: xp, yp, zp, rp
    real(c_real) :: pt_x, pt_y, pt_z

    ! local values used spring constants and damping coefficients
    real(c_real) :: ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W
    real(c_real) :: fnmd
    real(c_real) :: mag_overlap_t
    real(c_real) :: distmod_temp

    integer :: debug_id

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

       ! variable storing EB-particle overlap distance
       ! is there a FORTRAN MAX_DBL?
       distmod = 1.e20 ! HACK: larger number

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

       ! itterate over stencil
       do kk = klo, khi
          do jj = jlo, jhi
             do ii = ilo, ihi

                if (nbr(ii-i,jj-j,kk-k) .eq. 1) then
                ! if (.true.) then

                   ! if (is_single_valued_cell(flag(ii,jj,kk))) then
                   ! only consider cells that contain EB's
                   if ( .not. is_regular_cell(flag(ii,jj,kk))) then

                      ! convert bcent to global coordinate system centered at plo
                      bcentx = bcent(ii, jj, kk, 1)*dx(1) + (dble(ii) + 0.5d0)*dx(1)
                      bcenty = bcent(ii, jj, kk, 2)*dx(2) + (dble(jj) + 0.5d0)*dx(2)
                      bcentz = bcent(ii, jj, kk, 3)*dx(3) + (dble(kk) + 0.5d0)*dx(3)

                      ! Distance to boundary
                      distmod_temp = dabs( (xp - bcentx) * (-normal(ii,jj,kk,1)) + &
                                           (yp - bcenty) * (-normal(ii,jj,kk,2)) + &
                                           (zp - bcentz) * (-normal(ii,jj,kk,3)) )

                      ! Only keep the closest one
                      if (distmod_temp .lt. distmod) then
                         if (distmod_temp .lt. rp) then

                            ! Point where the normal from the particle to the plane intersects the plane
                            !   Note we go slightly less than the full normal in order to avoid precision problems
                            pt_x = xp + normal(ii,jj,kk,1) * distmod_temp * fudge
                            pt_y = yp + normal(ii,jj,kk,2) * distmod_temp * fudge
                            pt_z = zp + normal(ii,jj,kk,3) * distmod_temp * fudge

                            ! Cell that point is in
                            i_pt = floor(pt_x*inv_dx(1))
                            j_pt = floor(pt_y*inv_dx(2))
                            k_pt = floor(pt_z*inv_dx(3))

                            ! check if point-of-collision is in current cell, that's great!
                            ! Particle is interacting with "flat" EB's
                            if ( (i_pt .eq. ii) .and. (j_pt .eq. jj) .and. (k_pt .eq. kk)) then
                               distmod = distmod_temp
   
                               normul(1) = normal(ii,jj,kk,1)
                               normul(2) = normal(ii,jj,kk,2)
                               normul(3) = normal(ii,jj,kk,3)
                            else ! NOTE: particle can still interact with the edge (coner) of neighbour EB's
                                ! this should be moved to a subroutine
                                block
                                    ! variables keepint track of edges/corner of neighbour EB's
                                    integer, dimension(3) :: ind_loop, ind_pt, ind_facets
                                    ! n_facets: number of possible sides of the cell containing an EB edge
                                    ! tmp_facet: current facet/edge being tested
                                    ! cur_facet: best-guess of facet most likely a particle-EB edge collision 
                                    ! index_cell: cell index (along direction being tested) in loop
                                    ! ind_nb: neighbour cell index (along direction being tested)
                                    ! i_dim: direction (1,2, or 3) being tested
                                    integer :: n_facets, cur_facet, tmp_facet, ind_cell, ind_nb, i_dim, i_facet
                                    
                                    ! variables keeping track of coordinates on EB edges
                                    ! lambda_tmp: current lambda-value being used in testing for edge collions
                                    ! lambda: minimum (closets to bcentre) lambda value satisfying potential collision
                                    real(c_real) :: v_c, f_c, p_c, lambda, lambda_tmp
                                    ! b_vec: vector representation of EB-facet centre in cell
                                    ! p_vec: vector representation of EB-plain/collision point
                                    ! v_vec: vector from b_vec to p_vec
                                    ! c_vec: closest point (to particle) on line from b_vec to p_vec still in cell
                                    real(c_real), dimension(3) :: p_vec, b_vec, v_vec, c_vec

                                    ! assing value to vectors (arrays)
                                    ind_loop = (/ ii, jj, kk /)
                                    ind_pt   = (/ i_pt, j_pt, k_pt /)

                                    ! p_vec is outside cell, determine with cell boundaries might be cut
                                    call determine_cut_facets( ind_pt, ind_loop, ind_facets, n_facets)

                                    ! there could be several candidates, itterate picking the smallest lambda
                                    lambda = 1.e20 ! HACK: giant number!
                                    if (n_facets .gt. 0) then
                                        ! assign vectorial quantities
                                        p_vec = (/ pt_x, pt_y, pt_z /)
                                        b_vec = (/ bcentx, bcenty, bcentz /)

                                        ! v_vec = p_vec - b_vec
                                        call vector_sub_3d_real (v_vec, p_vec, b_vec)
                                        
                                        ! find minimal lamba
                                        do i_facet = 1, n_facets
                                            ! current guess
                                            tmp_facet = ind_facets(i_facet)
                                            
                                            ! cell cell indices for current cell (loop index) and neighbour
                                            ! in the current direction (facet) being tested
                                            ind_cell = ind_loop(tmp_facet)
                                            ind_nb = ind_pt(tmp_facet)
                                            
                                            ! the call face being tested depend on if the neighbour cell is to the
                                            ! "left" or the "right" of the current (loop) cell
                                            if (ind_cell .lt. ind_nb) then
                                                f_c = ( dble(ind_cell) + 1.0 )*dx(tmp_facet)
                                            else ! if (ind_cell .gt. ind_nb) then
                                                f_c = dble(ind_cell)*dx(tmp_facet)
                                            end if
                                            
                                            ! determin the lambda-value corresponding to the interception between the
                                            ! line b_vec->p_vec and the cell face being tested
                                            p_c = b_vec(tmp_facet)
                                            v_c = v_vec(tmp_facet)
                                            lambda_tmp = (f_c - p_c)/v_c
                                            
                                            ! keep minimal lambda-value and correponding direction (facet) index
                                            if (lambda_tmp .lt. lambda) then
                                                lambda = lambda_tmp
                                                cur_facet = tmp_facet
                                            end if
                                        end do

                                        ! now compute the coordinate of this interception point
                                        ! this will be closest to bcentre and is guaranteed to be int the same cell
                                        ! that "owns" the EB-facet.
                                        do i_dim = 1, 3
                                            if (i_dim .eq. cur_facet) then
                                                c_vec(i_dim) = f_c
                                            else
                                                c_vec(i_dim) = lambda*v_vec(i_dim) + b_vec(i_dim)
                                            end if
                                        end do
                                        
                                        ! correct collision coordinates
                                        pt_x = c_vec(1)
                                        pt_y = c_vec(2)
                                        pt_z = c_vec(3)

                                        ! check if there is still overlap
                                        distmod_temp = ((xp - pt_x) * ( xp - pt_x ) + & 
                                            (yp - pt_y) * ( yp - pt_y ) + &
                                            (zp - pt_z) * ( zp - pt_z ))
                                        if (distmod_temp .lt. distmod) then 
                                            if (distmod_temp .lt. rp*rp) then
                                                distmod = sqrt(distmod_temp)
                                                
                                                ! this normal is wrong, but it will do for now...
                                                normul(1) = -(xp - pt_x) / distmod ! normal(ii,jj,kk,1)
                                                normul(2) = -(yp - pt_y) / distmod ! normal(ii,jj,kk,2)
                                                normul(3) = -(zp - pt_z) / distmod ! normal(ii,jj,kk,3)

                                            end if
                                        end if
                                    end if
                                end block
                            end if
   
                         end if ! if dist lt rp
                      end if ! if dist lt distmod

                   end if ! if single_valued

                end if ! if nbr = 1

             end do
          end do
       end do

       if (distmod .lt. p%radius) then

          ! Calculate the particle/wall overlap.
          overlap_n = rp - distmod

          ! *****************************************************************************
          ! Calculate the translational relative velocity

          call cfrelvel_wall(ll, v_rel_trans_norm, vrel_t, normul, distmod, particles)

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
          TOW(LL,:) = TOW(LL,:) + distmod*DES_CROSSPRDCT(normul(:),FT)

       end if ! if test on (d < radius)

    end do ! loop over particles
   
   contains
       ! subroutine assign_vector_3d_int (vec, v_i, v_j, v_k)
       !     implicit none
       !
       !     integer, dimension(3), intent(  out) :: vec
       !     integer,               intent(in   ) :: v_i, v_j, v_k
       ! 
       !     vec(1) = v_i
       !     vec(2) = v_j
       !     vec(3) = v_k
       !
       ! end subroutine assign_vector_3d_int

       ! subroutine assign_vector_3d_real (vec, v_i, v_j, v_k)
       !     implicit none
       ! 
       !     real(c_real), dimension(3), intent(  out) :: vec
       !     real(c_real),               intent(in   ) :: v_i, v_j, v_k
       ! 
       !     vec(1) = v_i
       !     vec(2) = v_j
       !     vec(3) = v_k
       ! 
       ! end subroutine assign_vector_3d_real

       pure subroutine vector_sub_3d_real (res, a, b)
           implicit none

           real(c_real), dimension(3), intent(  out) :: res
           real(c_real), dimension(3), intent(in   ) :: a
           real(c_real), dimension(3), intent(in   ) :: b

           res(1) = a(1) - b(1)
           res(2) = a(2) - b(2)
           res(3) = a(3) - b(3)

       end subroutine vector_sub_3d_real

       pure subroutine determine_cut_facets ( ind_pt, ind_loop, ind_facets, n_facets )
           implicit none

           integer, dimension(1:3), intent(in   ) :: ind_pt
           integer, dimension(1:3), intent(in   ) :: ind_loop
           integer, dimension(1:3), intent(  out) :: ind_facets
           integer,                 intent(  out) :: n_facets

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

       end subroutine determine_cut_facets


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
