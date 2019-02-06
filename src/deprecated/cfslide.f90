MODULE CFSLIDE_MODULE

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int
   use param, only: zero, equal

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFSLIDE(V_TANG, PARTICLE_SLIDE, MU)
!  Purpose:  Check for Coulombs friction law - calculate sliding
!            friction
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer:                                          Date:
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFSLIDE(V_TANG, PARTICLE_SLIDE, MU, FT_tmp, FN_tmp)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! tangent to the plane of contact
      real(rt), INTENT(IN) :: V_TANG(3)
! logic set to T when a sliding contact occurs
      logical, INTENT(OUT) :: PARTICLE_SLIDE
! Coefficient of friction
      real(rt), INTENT(IN) :: MU
! normal force
      real(rt), DIMENSION(3), INTENT(IN) :: FN_tmp
! tangential force
      real(rt), DIMENSION(3), INTENT(INOUT) :: FT_tmp
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! squared magnitude of tangential and normal forces
      real(rt) FTMD, FNMD
!-----------------------------------------------

      FTMD = dot_product(FT_tmp(:),FT_tmp(:))
      FNMD = dot_product(FN_tmp(:),FN_tmp(:))

      IF (FTMD.GT.(MU*MU*FNMD)) THEN
! tangential force based on sliding friction
         PARTICLE_SLIDE = .TRUE.
         IF(EQUAL(MINVAL(ABS(V_TANG)), ZERO)) THEN
            FT_tmp(:) =  MU * FT_tmp(:) * SQRT(FNMD/FTMD)
         ELSE
            FT_tmp(:) = -MU * V_TANG(:) * &
               SQRT(FNMD/dot_product(V_TANG,V_TANG))
         ENDIF


      ENDIF

      RETURN
      END SUBROUTINE CFSLIDE
   END MODULE CFSLIDE_MODULE
