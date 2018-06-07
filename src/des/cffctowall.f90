MODULE CFFCTOWALL_MODULE

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFFCTOWALL
!  Purpose: Calculate the total force and torque on a particle
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer: Rahul Garg                               Date: 02-Aug-07!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE CFFCTOWALL(L, NORM, DIST_LI, FN, FT, DES_RADIUS, FC, TOW)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use discretelement, only: des_crossprdct
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! particle index
      integer, INTENT(IN) :: L
! distance between particle center and wall
      real(rt), INTENT(IN) :: DIST_LI
! unit normal vector along the line of contact pointing from
! particle L to wall
      real(rt), INTENT(IN) :: NORM(3)
! normal and tangential force
      real(rt), INTENT(IN) :: FN(3), FT(3)

      real(rt), DIMENSION(:), INTENT(IN) :: DES_RADIUS
      real(rt), DIMENSION(:,:), INTENT(INOUT) :: FC, TOW
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variable for calculating torque on particle
      real(rt) :: CROSSP(3)
! distance from the contact point to the particle center
      real(rt) DIST_CL
!------------------------------------------------

! total contact force
      FC(L,:) = FC(L,:) + FN(:) + FT(:)

! calculate the distance from the particle center to the wall
      DIST_CL = DIST_LI - DES_RADIUS(L)

! total torque
      CROSSP = DES_CROSSPRDCT(NORM, FT)
      TOW(:,L) = TOW(:,L) + DIST_CL*CROSSP(:)

      RETURN
      END SUBROUTINE CFFCTOWALL

END MODULE CFFCTOWALL_MODULE
