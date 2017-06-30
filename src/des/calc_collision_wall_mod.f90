!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_COLLISION_WALL                                    C
!  Author: Rahul Garg                               Date: 1-Dec-2013   C
!                                                                      C
!  Purpose: subroutines for particle-wall collisions when cutcell is   C
!           used. Also contains rehack of routines for cfslide and     C
!           cffctow which might be different from the stand alone      C
!           routines. Eventually all the DEM routines will be          C
!           consolidated.                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
module calc_collision_wall

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   use discretelement, only: des_coll_model_enum
   use discretelement, only: des_etat_wall, des_etan_wall, hert_kwn, hert_kwt, hertzian
   use discretelement, only: des_crossprdct
   use discretelement, only: kn_w, kt_w, mew_w!, dtsolid
   use error_manager, only: err_msg, flush_err_msg, init_err_msg
   use param, only: small_number, zero

   use stl_functions_des, only: closestptpointtriangle
   use discretelement, only: normal_particle

   implicit none
   private
   

   public :: calc_dem_force_with_wall_stl
   
contains


   subroutine calc_dem_force_with_wall_stl ( particles, fc, tow, xlength, ylength, zlength, dtsolid )

      use bc,            only: cyclic_x, cyclic_y, cyclic_z
      use param,         only: zero, one
      use particle_mod,  only: particle_t
      

      type(particle_t), intent(in)    :: particles(:)
      real(c_real)  ,   intent(inout) :: fc(:,:), tow(:,:)
      real(c_real)  ,   intent(in   ) :: xlength, ylength, zlength, dtsolid 
     
      integer :: ll, nf

      real(c_real) ::OVERLAP_N, SQRT_OVERLAP

      real(c_real) :: V_REL_TRANS_NORM, DISTSQ, RADSQ, CLOSEST_PT(3)
      ! local normal and tangential forces
      real(c_real) :: NORMAL(3), VREL_T(3), DIST(3), DISTMOD
      real(c_real) :: FT(3), FN(3), OVERLAP_T(3)

      integer :: PHASELL

      real(c_real) :: TANGENT(3)
      real(c_real) :: FNMD
      ! local values used spring constants and damping coefficients
      real(c_real) ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W

      real(c_real) :: MAG_OVERLAP_T

      real(c_real) :: line_t
      ! flag to tell if the orthogonal projection of sphere center to
      ! extended plane detects an overlap

      real(c_real) :: MAX_DISTSQ
      integer :: MAX_NF
      real(c_real), dimension(3) :: PARTICLE_MIN, PARTICLE_MAX, POS_TMP
      !     Vertex Coordinates X ,Y and Z
      real(c_real), dimension(3,3,6) :: VERTEX
      !     Face normal vector (normalized)
      real(c_real), dimension(3,6) :: NORM_FACE

      ! Skip this routine if the system is fully periodic.
      if(cyclic_x .and. cyclic_y .and. cyclic_z) return

      ! West Face
      VERTEX(1,:,1) = (/ZERO, ZERO, ZERO/)
      VERTEX(2,:,1) = (/ZERO, 2*YLENGTH, ZERO/)
      VERTEX(3,:,1) = (/ZERO, ZERO, 2*ZLENGTH/)
      NORM_FACE(:,1) = (/ONE, ZERO, ZERO/)

      ! East Face
      VERTEX(1,:,2) = (/XLENGTH, ZERO, ZERO/)
      VERTEX(2,:,2) = (/XLENGTH, 2*YLENGTH, ZERO/)
      VERTEX(3,:,2) = (/XLENGTH, ZERO, 2*ZLENGTH/)
      NORM_FACE(:,2) = (/-ONE, ZERO, ZERO/)

      ! South Face
      VERTEX(1,:,3) = (/ZERO, ZERO, ZERO/)
      VERTEX(2,:,3) = (/2*XLENGTH, ZERO, ZERO/)
      VERTEX(3,:,3) = (/ZERO, ZERO, 2*ZLENGTH/)
      NORM_FACE(:,3) = (/ZERO, ONE, ZERO/)

      ! North Face
      VERTEX(1,:,4) = (/ZERO, YLENGTH, ZERO/)
      VERTEX(2,:,4) = (/2*XLENGTH, YLENGTH, ZERO/)
      VERTEX(3,:,4) = (/ZERO, YLENGTH, 2*ZLENGTH/)
      NORM_FACE(:,4) = (/ZERO, -ONE, ZERO/)

      ! Bottom Face
      VERTEX(1,:,5) = (/ZERO, ZERO, ZERO/)
      VERTEX(2,:,5) = (/2*XLENGTH, ZERO, ZERO/)
      VERTEX(3,:,5) = (/ZERO, 2*YLENGTH, ZERO/)
      NORM_FACE(:,5) = (/ZERO, ZERO, ONE/)

      ! Top Face
      VERTEX(1,:,6) = (/ZERO, ZERO, ZLENGTH/)
      VERTEX(2,:,6) = (/2*XLENGTH, ZERO, ZLENGTH/)
      VERTEX(3,:,6) = (/ZERO, 2*YLENGTH, ZLENGTH/)
      NORM_FACE(:,6) = (/ZERO, ZERO, -ONE/)

      do ll = 1, size( particles )


         associate ( radius => particles(ll) % radius, pos => particles(ll) % pos, &
              & vel => particles(ll) % vel, omega => particles(ll) % omega )
        
            ! skipping non-existent particles or ghost particles
            if ( .not. ( NORMAL_PARTICLE == particles(ll) % state ) ) cycle

            ! Check particle LL for wall contacts
            radsq = radius * radius

            particle_max(:) = pos + radius
            particle_min(:) = pos - radius

            do nf = 1, 6

               if ( nf == 1 ) then
                  if ( pos(1) >  radius ) cycle
               else if ( nf == 2) then
                  if ( pos(1) < ( xlength - radius ) ) cycle
               else if ( nf == 3 ) then
                  if ( pos(2) > radius ) cycle
               else if ( nf == 4 ) then
                  if ( pos(2) < ( ylength - radius ) ) cycle
               else if ( nf == 5 ) then
                  if ( pos(3) > radius ) cycle
               else if ( nf == 6 ) then
                  if ( pos(3) < ( zlength - radius ) ) cycle
               end if


               ! Checking all the facets is time consuming due to the expensive
               ! separating axis test. Remove this facet from contention based on
               ! a simple orthogonal projection test.
               
               ! Parametrize a line as p = p_0 + t normal and intersect with the
               ! triangular plane. If t>0, then point is on the non-fluid side of
               ! the plane. If the plane normal is assumed to point toward the fluid.
               
               ! -undefined, because non zero values will imply the sphere center
               ! is on the non-fluid side of the plane. Since the testing
               ! is with extended plane, this could very well happen even
               ! when the particle is well inside the domain (assuming the plane
               ! normal points toward the fluid). See the pic below. So check
               ! only when line_t is negative
               
               !                 \   Solid  /
               !                  \  Side  /
               !                   \      /
               !                    \    /
               ! Wall 1, fluid side  \  /  Wall 2, fluid side
               !                      \/
               !                        o particle
               !
               ! line_t will be positive for wall 1 (incorrectly indicating center
               ! is outside the domain) and line_t will be negative for wall 2.
               !
               ! Therefore, only stick with this test when line_t is negative and let
               ! the separating axis test take care of the other cases.
               
               ! Since this is for checking static config, line's direction is the
               ! same as plane's normal. For moving particles, the line's normal will
               ! be along the point joining new and old positions.
               
               line_t = DOT_product( ( VERTEX(1,:,nf) - pos ), NORM_FACE(:,nf) )

               ! k - rad >= tol_orth, where k = -line_t, then orthogonal
               ! projection is false. Substituting for k
               ! => line_t + rad <= -tol_orth
               ! choosing tol_orth = 0.01% of des_radius = 0.0001*des_radius

               ! Orthogonal projection will detect false positives even
               ! when the particle does not overlap the triangle.
               ! However, if the orthogonal projection shows no overlap, then
               ! that is a big fat negative and overlaps are not possible.
               if ( line_t  <= ( -1.0001d0 * radius ) ) cycle

               POS_TMP = pos
               
               call ClosestPtPointTriangle( POS_TMP,  VERTEX(:,:,nf), CLOSEST_PT(:) )

               DIST(:) = CLOSEST_PT(:) - pos
               DISTSQ = DOT_PRODUCT(DIST, DIST)

               IF(DISTSQ .GE. RADSQ - SMALL_NUMBER) CYCLE

               MAX_DISTSQ = DISTSQ
               MAX_NF = NF

               ! Assign the collision normal based on the facet with the
               ! largest overlap.
               NORMAL(:) = DIST(:)/sqrt(DISTSQ)

               ! Facet's normal is correct normal only when the intersection is with
               ! the face. When the intersection is with edge or vertex, then the
               ! normal is based on closest pt and sphere center. The definition above
               ! of the normal is generic enough to account for differences between
               ! vertex, edge, and facet.

               ! Calculate the particle/wall overlap.
               DISTMOD = sqrt(MAX_DISTSQ)
               OVERLAP_N = radius - DISTMOD

               ! Calculate the translational relative velocity
               call cfrelvel_wall(ll, v_rel_trans_norm, vrel_t, normal, distmod, &
                    & particles )

               ! Calculate the spring model parameters.
               phaseLL = particles(ll) % phase

               ! Hertz vs linear spring-dashpot contact model
               if ( DES_COLL_MODEL_ENUM == HERTZIAN ) then
                  sqrt_overlap = sqrt(OVERLAP_N)
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
               FN(:) = -(KN_DES_W * OVERLAP_N * NORMAL(:) + &
                    ETAN_DES_W * V_REL_TRANS_NORM * NORMAL(:))

               ! Calculate the tangential displacement.
               overlap_t(:) = dtsolid*vrel_t(:)
               mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t))

               ! Check for Coulombs friction law and limit the maximum value of the
               ! tangential force on a particle in contact with a wall.
               if ( MAG_OVERLAP_T > 0.0 ) then
                  ! Max force before the on set of frictional slip.
                  FNMD = MEW_W*sqrt(DOT_PRODUCT(FN,FN))
                  ! Direction of tangential force.
                  TANGENT = OVERLAP_T/MAG_OVERLAP_T
                  FT = -FNMD * TANGENT
               else
                  FT = 0.0
               end if

               ! Add the collision force to the total forces acting on the particle.
               FC(LL,:) = FC(LL,:) + FN(:) + FT(:)

               ! Add the torque force to the total torque acting on the particle.
               TOW(LL,:) = TOW(LL,:) + DISTMOD*DES_CROSSPRDCT(NORMAL,FT)
               
            end do

         end associate

      end do

   end subroutine calc_dem_force_with_wall_stl



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
      ! unit normal from particle center to closest point on stl (wall)
      real(c_real), intent(in) :: norm(3)
      ! distance between particle center and stl (wall).
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

   
end module CALC_COLLISION_WALL
