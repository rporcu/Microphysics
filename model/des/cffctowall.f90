MODULE CFFCTOWALL_MODULE

   use amrex_fort_module, only : c_real => amrex_real
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
      USE discretelement, only: des_crossprdct
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! particle index
      INTEGER, INTENT(IN) :: L
! distance between particle center and wall
      real(c_real), INTENT(IN) :: DIST_LI
! unit normal vector along the line of contact pointing from
! particle L to wall
      real(c_real), INTENT(IN) :: NORM(3)
! normal and tangential force
      real(c_real), INTENT(IN) :: FN(3), FT(3)

      real(c_real), DIMENSION(:), INTENT(IN) :: DES_RADIUS
      real(c_real), DIMENSION(:,:), INTENT(INOUT) :: FC, TOW
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variable for calculating torque on particle
      real(c_real) :: CROSSP(3)
! distance from the contact point to the particle center
      real(c_real) DIST_CL
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
