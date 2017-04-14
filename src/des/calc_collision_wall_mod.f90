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
      MODULE CALC_COLLISION_WALL

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use discretelement, only: des_coll_model_enum
      use discretelement, only: des_etat_wall, des_etan_wall, hert_kwn, hert_kwt, hertzian
      use discretelement, only: des_crossprdct
      use discretelement, only: kn_w, kt_w, mew_w, dtsolid
      use error_manager, only: err_msg, flush_err_msg, init_err_msg
      use param, only: small_number, zero

      use stl_functions_des, only: closestptpointtriangle
      use discretelement, only: normal_particle

      PRIVATE
      PUBLIC :: CALC_DEM_FORCE_WITH_WALL_STL

      CONTAINS

!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL(particle_phase, &
            particle_state,  des_radius, des_pos_new, &
            des_vel_new, omega_new, fc, tow, xlength, ylength, zlength)

      use bc, only: cyclic_x, cyclic_y, cyclic_z
      use param, only: zero, one

      IMPLICIT NONE

      real(c_real)  , intent(in   ) :: des_radius(:)
      real(c_real)  , intent(in   ) :: des_pos_new(:,:), des_vel_new(:,:), omega_new(:,:)
      real(c_real)  , intent(inout) :: fc(:,:), tow(:,:)
      real(c_real)  , intent(in   ) :: xlength, ylength, zlength
      integer(c_int), intent(in   ) :: particle_state(:)
      integer(c_int), intent(in   ) :: particle_phase(:)

      integer :: LL
      integer :: NF
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
      real(c_real), DIMENSION(3) :: PARTICLE_MIN, PARTICLE_MAX, POS_TMP
!     Vertex Coordinates X ,Y and Z
      real(c_real), DIMENSION(3,3,6) :: VERTEX
!     Face normal vector (normalized)
      real(c_real), DIMENSION(3,6) :: NORM_FACE

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

      DO LL = 1, size(des_radius)

! skipping non-existent particles or ghost particles
! make sure the particle is not classified as a new 'entering' particle
! or is already marked as a potential exiting particle
         IF(.NOT.NORMAL_PARTICLE==PARTICLE_STATE(LL)) CYCLE

! Check particle LL for wall contacts
         RADSQ = DES_RADIUS(LL)*DES_RADIUS(LL)

         particle_max(:) = des_pos_new(LL,:) + des_radius(LL)
         particle_min(:) = des_pos_new(LL,:) - des_radius(LL)

         DO NF = 1, 6

            if(nf==1    ) then
               if(des_pos_new(ll,1) > des_radius(LL)) cycle
            elseif(nf==2) then
               if(des_pos_new(ll,1) < xlength - des_radius(LL)) cycle
            elseif(nf==3) then
               if(des_pos_new(ll,2) > des_radius(LL)) cycle
            elseif(nf==4) then
               if(des_pos_new(ll,2) < ylength - des_radius(LL)) cycle
            elseif(nf==5) then
               if(des_pos_new(ll,3) > des_radius(LL)) cycle
            elseif(nf==6) then
               if(des_pos_new(ll,3) < zlength - des_radius(LL)) cycle
            endif


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

            line_t = DOT_PRODUCT(VERTEX(1,:,NF) - des_pos_new(LL,:),&
               NORM_FACE(:,NF))

! k - rad >= tol_orth, where k = -line_t, then orthogonal
! projection is false. Substituting for k
! => line_t + rad <= -tol_orth
! choosing tol_orth = 0.01% of des_radius = 0.0001*des_radius

! Orthogonal projection will detect false positives even
! when the particle does not overlap the triangle.
! However, if the orthogonal projection shows no overlap, then
! that is a big fat negative and overlaps are not possible.
            if((line_t.le.-1.0001d0*des_radius(LL))) CYCLE

            POS_TMP = DES_POS_NEW(LL,:)
            CALL ClosestPtPointTriangle(POS_TMP, &
               VERTEX(:,:,NF), CLOSEST_PT(:))

            DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(LL,:)
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
            DISTMOD = SQRT(MAX_DISTSQ)
            OVERLAP_N = DES_RADIUS(LL) - DISTMOD

! Calculate the translational relative velocity
            CALL CFRELVEL_WALL(LL, V_REL_TRANS_NORM,VREL_T,           &
               NORMAL, DISTMOD, DES_VEL_NEW, OMEGA_NEW)

! Calculate the spring model parameters.
            phaseLL = particle_phase(LL)

! Hertz vs linear spring-dashpot contact model
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               sqrt_overlap = SQRT(OVERLAP_N)
               KN_DES_W = hert_kwn(phaseLL)*sqrt_overlap
               KT_DES_W = hert_kwt(phaseLL)*sqrt_overlap
               sqrt_overlap = SQRT(sqrt_overlap)
               ETAN_DES_W = DES_ETAN_WALL(phaseLL)*sqrt_overlap
               ETAT_DES_W = DES_ETAT_WALL(phaseLL)*sqrt_overlap
            ELSE
               KN_DES_W = KN_W
               KT_DES_W = KT_W
               ETAN_DES_W = DES_ETAN_WALL(phaseLL)
               ETAT_DES_W = DES_ETAT_WALL(phaseLL)
            ENDIF

! Calculate the normal contact force
            FN(:) = -(KN_DES_W * OVERLAP_N * NORMAL(:) + &
               ETAN_DES_W * V_REL_TRANS_NORM * NORMAL(:))

! Calculate the tangential displacement.
            OVERLAP_T(:) = DTSOLID*VREL_T(:)
            MAG_OVERLAP_T = sqrt(DOT_PRODUCT(OVERLAP_T, OVERLAP_T))

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall.
            IF(MAG_OVERLAP_T > 0.0) THEN
! Max force before the on set of frictional slip.
               FNMD = MEW_W*sqrt(DOT_PRODUCT(FN,FN))
! Direction of tangential force.
               TANGENT = OVERLAP_T/MAG_OVERLAP_T
               FT = -FNMD * TANGENT
            ELSE
               FT = 0.0
            ENDIF

! Add the collision force to the total forces acting on the particle.
            FC(LL,:) = FC(LL,:) + FN(:) + FT(:)

! Add the torque force to the toal torque acting on the particle.
            TOW(LL,:) = TOW(LL,:) + DISTMOD*DES_CROSSPRDCT(NORMAL,FT)

         ENDDO

      ENDDO

      RETURN

      END SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL


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
      SUBROUTINE CFRELVEL_WALL(LL, VRN, VRT, NORM, DIST, DES_VEL_NEW, OMEGA_NEW)

! Function for calculating the cross prodcut
      use discretelement, only: DES_CROSSPRDCT

      IMPLICIT NONE

! Particle translational velocity
      real(c_real), DIMENSION(:,:), INTENT(IN) :: DES_VEL_NEW
! Particle rotational velocity
      real(c_real), DIMENSION(:,:), INTENT(IN) :: OMEGA_NEW

! Dummy arguments:
!---------------------------------------------------------------------//
! Particle index.
      integer, INTENT(IN) :: LL
! Magnitude of the total relative translational velocity.
      real(c_real), INTENT(OUT):: VRN
! Total relative translational velocity (vector).
      real(c_real), INTENT(OUT):: VRT(3)
! Unit normal from particle center to closest point on stl (wall)
      real(c_real), INTENT(IN) :: NORM(3)
! Distance between particle center and stl (wall).
      real(c_real), INTENT(IN) :: DIST

! Local variables
!---------------------------------------------------------------------//
! Additional relative translational motion due to rotation
      real(c_real) :: V_ROT(3)
! Total relative velocity at contact point
      real(c_real) :: VRELTRANS(3)

! Total relative velocity + rotational contribution
      V_ROT = DIST*OMEGA_NEW(LL,:)
      VRELTRANS(:) =  DES_VEL_NEW(LL,:) + DES_CROSSPRDCT(V_ROT, NORM)

! magnitude of normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! total relative translational slip velocity at the contact point
! Equation (8) in Tsuji et al. 1992
      VRT(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN
      END SUBROUTINE CFRELVEL_WALL

      end module CALC_COLLISION_WALL
