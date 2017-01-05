MODULE CALC_D_MOD

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_d_n                                                !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!           pressure correction -- North                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_D(D, AXIS, A_M, ep_g, f_gds, flag, dx, dy, dz)

! Global Variables:
!---------------------------------------------------------------------//
! Flag: Coupled DEM simulation
      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: DES_ONEWAY_COUPLED
! Pressure scale factor
      use scales, only: P_SCALE

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and size of solids phase arrays.

      use compar, only: istart2,iend2,jstart2,jend2,kstart2,kend2
      use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3

      use functions, only: AVG
      use param1, only: ZERO, SMALL_NUMBER
      use functions, only: ieast, jnorth, ktop

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Pressure correction
      real(c_real), INTENT(OUT) :: d(:,:,:)
! "X", "Y", or "Z"
      CHARACTER, INTENT(IN) :: axis
      real(c_real), INTENT(IN):: A_M&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
      real(c_real), INTENT(IN):: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
      real(c_real), INTENT(IN   ) :: dx, dy, dz

! Local variables:
!---------------------------------------------------------------------//
! Usual Indices
      INTEGER :: I,J,K
      real(c_real) :: axy, axz, ayz, vol

      ! Temp variable 
      real(c_real) :: AM0
      real(c_real) :: EPGA
      logical :: COUPLED
!......................................................................!

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz
      vol = dx*dy*dz

      COUPLED = (DES_CONTINUUM_COUPLED .AND. .NOT.DES_ONEWAY_COUPLED)

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

         AM0 = -A_M(I,J,K,0)

         if (axis.eq.'X') then
            if(flag(i,j,k,2) >= 2000 .and. &
               flag(i,j,k,2) <= 2011) then
               EPGA = AYZ*AVG(EP_G(I,J,K),EP_G(ieast(i,j,k),j,k))
               IF(COUPLED) AM0 = AM0 + 0.5d0*VOL* &
                  (F_GDS(i,j,k) + F_GDS(ieast(i,j,k),j,k))
            ELSE
               EPGA = ZERO
            ENDIF

         else if (axis.eq.'Y') then
            if(flag(i,j,k,3) >= 2000 .and. &
               flag(i,j,k,3) <= 2011) then
               EPGA = AXZ*AVG(EP_G(I,J,K),EP_G(i,jnorth(i,j,k),k))
               IF(COUPLED) AM0 = AM0 + 0.5d0*VOL* &
                  (F_GDS(i,j,k) + F_GDS(i,jnorth(i,j,k),k))
            ELSE
               EPGA = ZERO
            ENDIF

         else if (axis.eq.'Z') then
            if(flag(i,j,k,4) >= 2000 .and. &
               flag(i,j,k,4) <= 2011) then
               EPGA = AXY*AVG(EP_G(I,J,K),EP_G(i,j,ktop(i,j,k)))
               IF(COUPLED) AM0 = AM0 + 0.5d0*VOL* &
                  (F_GDS(I,J,K) + F_GDS(i,j,ktop(i,j,k)))
            ELSE
               EPGA = ZERO
            ENDIF
         endif

         IF(abs(AM0) > SMALL_NUMBER) THEN
            D(I,J,K) = P_SCALE*EPGA/AM0

         ELSE
            D(I,J,K) = ZERO
         ENDIF

      ENDDO
      ENDDO
      ENDDO

   END SUBROUTINE CALC_D

   END MODULE CALC_D_MOD
