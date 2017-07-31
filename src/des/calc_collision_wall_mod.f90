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

      real(c_real) ::overlap_n, SQRT_OVERLAP

      real(c_real) :: v_rel_trans_norm, distsq, radsq, closest_pt(3)

      ! local normal and tangential forces
      real(c_real) :: normal(3), vrel_t(3), dist(3), distmod
      real(c_real) :: ft(3), fn(3), overlap_t(3)

      integer :: PHASELL

      real(c_real) :: tangent(3)
      real(c_real) :: fnmd
      ! local values used spring constants and damping coefficients
      real(c_real) ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W

      real(c_real) :: MAG_OVERLAP_T

      real(c_real) :: line_t
      ! flag to tell if the orthogonal projection of sphere center to
      ! extended plane detects an overlap

      real(c_real) :: MAX_distsq
      integer :: MAX_NF
      real(c_real), dimension(3) :: PARTICLE_MIN, PARTICLE_MAX, POS_TMP
      !     Vertex Coordinates X ,Y and Z
      real(c_real), dimension(3,3,6) :: vertex
      !     Face normal vector (normalized)
      real(c_real), dimension(3,6) :: norm_face

      ! additional relative translational motion due to rotation
      real(c_real) :: v_rot(3)
      ! total relative velocity at contact point
      real(c_real) :: vreltrans(3)

      ! Skip this routine if the system is fully periodic.
      if (cyclic_x .and. cyclic_y .and. cyclic_z) return

      ! West Face
      vertex(1,:,1) = (/ZERO, ZERO, ZERO/)
      vertex(2,:,1) = (/ZERO, 2*YLENGTH, ZERO/)
      vertex(3,:,1) = (/ZERO, ZERO, 2*ZLENGTH/)
      norm_face(:,1) = (/ONE, ZERO, ZERO/)

      ! East Face
      vertex(1,:,2) = (/XLENGTH, ZERO, ZERO/)
      vertex(2,:,2) = (/XLENGTH, 2*YLENGTH, ZERO/)
      vertex(3,:,2) = (/XLENGTH, ZERO, 2*ZLENGTH/)
      norm_face(:,2) = (/-ONE, ZERO, ZERO/)

      ! South Face
      vertex(1,:,3) = (/ZERO, ZERO, ZERO/)
      vertex(2,:,3) = (/2*XLENGTH, ZERO, ZERO/)
      vertex(3,:,3) = (/ZERO, ZERO, 2*ZLENGTH/)
      norm_face(:,3) = (/ZERO, ONE, ZERO/)

      ! North Face
      vertex(1,:,4) = (/ZERO, YLENGTH, ZERO/)
      vertex(2,:,4) = (/2*XLENGTH, YLENGTH, ZERO/)
      vertex(3,:,4) = (/ZERO, YLENGTH, 2*ZLENGTH/)
      norm_face(:,4) = (/ZERO, -ONE, ZERO/)

      ! Bottom Face
      vertex(1,:,5) = (/ZERO, ZERO, ZERO/)
      vertex(2,:,5) = (/2*XLENGTH, ZERO, ZERO/)
      vertex(3,:,5) = (/ZERO, 2*YLENGTH, ZERO/)
      norm_face(:,5) = (/ZERO, ZERO, ONE/)

      ! Top Face
      vertex(1,:,6) = (/ZERO, ZERO, ZLENGTH/)
      vertex(2,:,6) = (/2*XLENGTH, ZERO, ZLENGTH/)
      vertex(3,:,6) = (/ZERO, 2*YLENGTH, ZLENGTH/)
      norm_face(:,6) = (/ZERO, ZERO, -ONE/)

      do ll = 1, size( particles )


         associate ( radius => particles(ll) % radius, pos => particles(ll) % pos, &
              & vel => particles(ll) % vel, omega => particles(ll) % omega )
        
         ! Check particle LL for wall contacts
         radsq = radius * radius

         do nf = 1, 6

               if ( nf == 1 .and. .not.cyclic_x) then
                  if ( pos(1) >  radius ) cycle
               else if ( nf == 2.and. .not.cyclic_x) then
                  if ( pos(1) < ( xlength - radius ) ) cycle
               else if ( nf == 3 .and. .not.cyclic_y) then
                  if ( pos(2) > radius ) cycle
               else if ( nf == 4 .and. .not.cyclic_y) then
                  if ( pos(2) < ( ylength - radius ) ) cycle
               else if ( nf == 5 .and. .not.cyclic_z) then
                  if ( pos(3) > radius ) cycle
               else if ( nf == 6 .and. .not.cyclic_z) then
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
               
               line_t = DOT_product( ( vertex(1,:,nf) - pos ), norm_face(:,nf) )

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
               
               call ClosestPtPointTriangle( POS_TMP,  vertex(:,:,nf), closest_pt(:) )

               dist(:) = closest_pt(:) - pos
               distsq = dot_product(dist, dist)

               if (distsq .ge. radsq - small_number) cycle

               max_distsq = distsq
               max_nf = nf

               ! Assign the collision normal based on the facet with the
               ! largest overlap.
               normal(:) = dist(:)/sqrt(distsq)

               ! Facet's normal is correct normal only when the intersection is with
               ! the face. When the intersection is with edge or vertex, then the
               ! normal is based on closest pt and sphere center. The definition above
               ! of the normal is generic enough to account for differences between
               ! vertex, edge, and facet.

               ! Calculate the particle/wall overlap.
               distmod = sqrt(MAX_distsq)
               overlap_n = radius - distmod

               ! *****************************************************************************
               ! Calculate the translational relative velocity
               call cfrelvel_wall(ll, v_rel_trans_norm, vrel_t, normal, distmod, particles)
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
               FN(:) = -(KN_DES_W * overlap_n * normal(:) + &
                    ETAN_DES_W * V_REL_TRANS_NORM * normal(:))

               ! Calculate the tangential displacement.
               overlap_t(:) = dtsolid*vrel_t(:)
               mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t))

               ! Check for Coulombs friction law and limit the maximum value of the
               ! tangential force on a particle in contact with a wall.
               if ( MAG_OVERLAP_T > 0.0 ) then
                  ! Max force before the on set of frictional slip.
                  FNMD = MEW_W*sqrt(dot_product(FN,FN))
                  ! Direction of tangential force.
                  tangent = OVERLAP_T/MAG_OVERLAP_T
                  FT = -FNMD * tangent
               else
                  FT = 0.0
               end if

               ! Add the collision force to the total forces acting on the particle.
               FC(LL,:) = FC(LL,:) + FN(:) + FT(:)

               ! Add the torque force to the total torque acting on the particle.
               TOW(LL,:) = TOW(LL,:) + distmod*DES_CROSSPRDCT(normal,FT)
               
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
