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
      SUBROUTINE CALC_D(slo, shi, lo, hi, D, axis, A_m, ep_g, f_gds, flag, dx, dy, dz)

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

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! Pressure correction
      real(c_real), intent(  out) :: d&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! "X", "Y", or "Z"
      character, intent(in   ) :: axis

      real(c_real), intent(in   ):: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3), -3:3)
      real(c_real), intent(in   ):: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer     , intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dx, dy, dz

!---------------------------------------------------------------------//
! Local variables:
!---------------------------------------------------------------------//
      integer      :: i,j,k,is,js,ks
      real(c_real) :: axy, axz, ayz, vol

      ! Temp variable 
      real(c_real) :: Am0
      real(c_real) :: EPGA
      logical :: COUPLED
!---------------------------------------------------------------------//

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz
      vol = dx*dy*dz

      COUPLED = (DES_CONTINUUM_COUPLED .AND. .NOT.DES_ONEWAY_COUPLED)

      is = lo(1)
      js = lo(2)
      ks = lo(3)

      if (axis.eq.'x') then
          is = lo(1)-1
      else if (axis.eq.'y') then
          js = lo(2)-1
      else if (axis.eq.'z') then
          ks = lo(3)-1
      endif
 
      DO K = ks, hi(3)
        DO J = js, hi(2)
          DO I = is, hi(1)

         Am0 = -A_m(I,J,K,0)

         if (axis.eq.'X') then
            if(flag(i,j,k,2) >= 2000 .and. &
               flag(i,j,k,2) <= 2011) then
               EPGA = AYZ*AVG(EP_G(I,J,K),EP_G(i+1,j,k))
               IF(COUPLED) Am0 = Am0 + 0.5d0*VOL* &
                  (F_GDS(i,j,k) + F_GDS(i+1,j,k))
            ELSE
               EPGA = ZERO
            ENDIF

         else if (axis.eq.'Y') then
            if(flag(i,j,k,3) >= 2000 .and. &
               flag(i,j,k,3) <= 2011) then
               EPGA = AXZ*AVG(EP_G(I,J,K),EP_G(i,j+1,k))
               IF(COUPLED) Am0 = Am0 + 0.5d0*VOL* &
                  (F_GDS(i,j,k) + F_GDS(i,j+1,k))
            ELSE
               EPGA = ZERO
            ENDIF

         else if (axis.eq.'Z') then
            if(flag(i,j,k,4) >= 2000 .and. &
               flag(i,j,k,4) <= 2011) then
               EPGA = AXY*AVG(EP_G(I,J,K),EP_G(i,j,k+1))
               IF(COUPLED) Am0 = Am0 + 0.5d0*VOL* &
                  (F_GDS(I,J,K) + F_GDS(i,j,k+1))
            ELSE
               EPGA = ZERO
            ENDIF
         endif

         IF(abs(Am0) > SMALL_NUMBER) THEN
            D(I,J,K) = P_SCALE*EPGA/Am0

         ELSE
            D(I,J,K) = ZERO
         ENDIF

      ENDDO
      ENDDO
      ENDDO

   END SUBROUTINE CALC_D

   END MODULE CALC_D_MOD
