module CFRELVEL_MODULE

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   implicit none
   private

   public CFRELVEL
   public cfrelvel_aos

   
   
contains
   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
   !
   !  Subroutine: CFRELVEL
   !  Purpose: Calculate the normal and tangential components of the
   !           relative velocity between contacting particles
   !
   !  Author: Jay Boyalakuntla                           Date: 12-Jun-04
   !  Reviewer: Rahul Garg                               Date: 01-Aug-07
   !
   !  Comments: Relative (translational) velocity required for eqn 6
   !  from the following paper:
   !    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
   !    simulation of plug glow of cohesionless particles in a
   !    horizontal pipe", Powder technology, 71, 239-250, 1992
   !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   SUBROUTINE CFRELVEL(L, II, VRN, VSLIP, NORM, DIST_LI, DES_VEL_NEW, DES_RADIUS, OMEGA_NEW)

      !-----------------------------------------------
      ! Modules
      !-----------------------------------------------
      use discretelement, only: DES_CROSSPRDCT
      IMPLICIT NONE

      real(c_real), DIMENSION(:), INTENT(IN) :: des_radius
      real(c_real), DIMENSION(:,:), INTENT(IN) :: des_vel_new, omega_new

      !-----------------------------------------------
      ! Dummy arguments
      !-----------------------------------------------
      ! indices of particle-particle contact pair
      integer, INTENT(IN) :: L, II
      ! distance between particle centers
      real(c_real), INTENT(IN) :: DIST_LI
      ! unit normal vector along the line of contact pointing from
      ! particle L to particle II
      real(c_real), INTENT(IN) :: NORM(3)
      ! slip velocity at point of contact
      real(c_real), INTENT(OUT) :: VSLIP(3)
      ! normal component of relative contact velocity (scalar)
      real(c_real), INTENT(OUT) :: VRN
      !-----------------------------------------------
      ! Local variables
      !-----------------------------------------------
      ! translational relative velocity
      real(c_real) :: VRELTRANS(3)
      ! rotational velocity at point of contact
      real(c_real) :: V_ROT(3), OMEGA_SUM(3)
      ! distance from the contact point to the particle centers
      real(c_real) :: DIST_CL, DIST_CI
      !-----------------------------------------------

      ! translational relative velocity
      VRELTRANS(:) = (DES_VEL_NEW(L,:) - DES_VEL_NEW(II,:))

      ! calculate the distance from the particle center to the contact point,
      ! which is taken as the radical line
      ! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
      DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
           (2.d0*DIST_LI)
      DIST_CI = DIST_LI - DIST_CL

      OMEGA_SUM(:) = OMEGA_NEW(L,:)*DIST_CL + &
           OMEGA_NEW(II,:)*DIST_CI

      ! calculate the rotational relative velocity
      V_ROT = DES_CROSSPRDCT(OMEGA_SUM, NORM)

      ! total relative velocity
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

      ! normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

      ! slip velocity of the contact point
      ! Equation (8) in Tsuji et al. 1992
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN
   END SUBROUTINE CFRELVEL



   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
   !
   !  Subroutine: CFRELVEL
   !  Purpose: Calculate the normal and tangential components of the
   !           relative velocity between contacting particles
   !
   !  Author: Jay Boyalakuntla                           Date: 12-Jun-04
   !  Reviewer: Rahul Garg                               Date: 01-Aug-07
   !
   !  Comments: Relative (translational) velocity required for eqn 6
   !  from the following paper:
   !    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
   !    simulation of plug glow of cohesionless particles in a
   !    horizontal pipe", Powder technology, 71, 239-250, 1992
   !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   subroutine cfrelvel_aos(L, II, VRN, VSLIP, NORM, DIST_LI, particles)

      !-----------------------------------------------
      ! Modules
      !-----------------------------------------------
      use discretelement, only: DES_CROSSPRDCT
      use particle_mod,   only: particle_t

      type(particle_t),  intent(in)  :: particles(:) 

      !-----------------------------------------------
      ! Dummy arguments
      !-----------------------------------------------
      ! indices of particle-particle contact pair
      integer, INTENT(IN) :: L, II
      ! distance between particle centers
      real(c_real), INTENT(IN) :: DIST_LI
      ! unit normal vector along the line of contact pointing from
      ! particle L to particle II
      real(c_real), INTENT(IN) :: NORM(3)
      ! slip velocity at point of contact
      real(c_real), INTENT(OUT) :: VSLIP(3)
      ! normal component of relative contact velocity (scalar)
      real(c_real), INTENT(OUT) :: VRN
      !-----------------------------------------------
      ! Local variables
      !-----------------------------------------------
      ! translational relative velocity
      real(c_real) :: VRELTRANS(3)
      ! rotational velocity at point of contact
      real(c_real) :: V_ROT(3), OMEGA_SUM(3)
      ! distance from the contact point to the particle centers
      real(c_real) :: DIST_CL, DIST_CI
      !-----------------------------------------------

      ! translational relative velocity
      VRELTRANS(:) = ( particles(l) % vel  -  particles(ii) % vel )

      ! calculate the distance from the particle center to the contact point,
      ! which is taken as the radical line
      ! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
      DIST_CL = ( DIST_LI**2 + particles(l) % radius**2 - particles(ii) % radius**2 )/&
           (2.d0*DIST_LI)
      DIST_CI = DIST_LI - DIST_CL

      OMEGA_SUM(:) = particles(l) % omega * DIST_CL + &
           particles(ii) % omega * DIST_CI

      ! calculate the rotational relative velocity
      V_ROT = DES_CROSSPRDCT(OMEGA_SUM, NORM)

      ! total relative velocity
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

      ! normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

      ! slip velocity of the contact point
      ! Equation (8) in Tsuji et al. 1992
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)

     
   end subroutine CFRELVEL_AOS


   
END MODULE CFRELVEL_MODULE
