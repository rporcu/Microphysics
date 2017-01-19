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
      SUBROUTINE CALC_D(slo, shi, lo, hi, D, AXIS, A_m, alo, ahi, ep_g, f_gds, flag, dx, dy, dz)

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

      use functions, only: AVG
      use param1, only: ZERO, SMALL_NUMBER
      use functions, only: ieast, jnorth, ktop

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)

      ! Pressure correction
      real(c_real), intent(  out) :: d&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! "X", "Y", or "Z"
      CHARACTER, intent(in   ) :: axis

      real(c_real), intent(in   ):: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3), -3:3)
      real(c_real), intent(in   ):: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer     , intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dx, dy, dz

! Local variables:
!---------------------------------------------------------------------//
! Usual Indices
      INTEGER :: i,j,k
      INTEGER :: ie,je,ke
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

      ie = hi(1)
      je = hi(2)
      ke = hi(3)

      if (axis.eq.'x') then
          ie = hi(1)+1
      else if (axis.eq.'y') then
          je = hi(2)+1
      else if (axis.eq.'z') then
          ke = hi(3)+1
      endif

      do K = lo(3),ke
        do J = lo(2),je
          do I = lo(1),ie

         AM0 = -A_m(I,J,K,0)

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
