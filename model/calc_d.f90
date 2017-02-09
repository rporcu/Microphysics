module calc_d_mod

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   use functions, only: AVG
   use param1, only: ZERO, SMALL_NUMBER

   ! Flag: Coupled DEM simulation
   use discretelement, only: DES_CONTINUUM_COUPLED
   use discretelement, only: DES_ONEWAY_COUPLED

   ! Pressure scale factor
   use scales, only: P_SCALE

   implicit none

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_d_*                                                !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!           pressure correction
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_d_e(slo, shi, ulo, uhi, lo, hi, d_e, A_m, &
                       ep_g, f_gds, flag, dx, dy, dz)

      integer, intent(in   ) :: slo(3),shi(3),ulo(3),uhi(3),lo(3),hi(3)

      ! Pressure correction
      real(c_real), intent(  out) :: d_e&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))

      real(c_real), intent(in   ):: A_m&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3), -3:3)

      real(c_real), intent(in   ):: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer     , intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dx, dy, dz

      integer      :: i,j,k
      real(c_real) :: ayz, vol
      real(c_real) :: Am0, epga

      logical      :: coupled

      COUPLED = (DES_CONTINUUM_COUPLED .AND. .NOT.DES_ONEWAY_COUPLED)

      ayz = dy*dz
      vol = dx*dy*dz

      DO K = lo(3), hi(3)
        DO J = lo(2), hi(2)
          DO I = lo(1)-1, hi(1)

            Am0 = -A_m(I,J,K,0)

            IF(abs(Am0) > SMALL_NUMBER) THEN

               if(flag(i,j,k,2) >= 2000 .and. &
                  flag(i,j,k,2) <= 2011) then
                  epga = ayz*AVG(EP_G(I,J,K),EP_G(i+1,j,k))
                  IF(COUPLED) Am0 = Am0 + 0.5d0*VOL* &
                     (F_GDS(i,j,k) + F_GDS(i+1,j,k))
               ELSE
                  epga = ZERO
               ENDIF

               d_e(I,J,K) = P_SCALE*epga/Am0

            ELSE
               d_e(I,J,K) = ZERO
            ENDIF

          ENDDO
        ENDDO
      ENDDO

   end subroutine calc_d_e

   subroutine calc_d_n(slo, shi, vlo, vhi, lo, hi, d_n, A_m, ep_g, f_gds, flag, dx, dy, dz)

      integer     , intent(in   ) :: slo(3),shi(3),vlo(3),vhi(3),lo(3),hi(3)

      ! Pressure correction
      real(c_real), intent(  out) :: d_n&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

      real(c_real), intent(in   ):: A_m&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3), -3:3)

      real(c_real), intent(in   ):: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer     , intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dx, dy, dz

      integer      :: i,j,k
      real(c_real) :: axy, axz, ayz, vol
      real(c_real) :: Am0, epga
      logical      :: coupled

      COUPLED = (DES_CONTINUUM_COUPLED .AND. .NOT.DES_ONEWAY_COUPLED)

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz
      vol = dx*dy*dz
 
      DO K = lo(3), hi(3)
        DO J = lo(2)-1, hi(2)
          DO I = lo(1), hi(1)

         Am0 = -A_m(I,J,K,0)

         IF(abs(Am0) > SMALL_NUMBER) THEN

            if(flag(i,j,k,3) >= 2000 .and. &
               flag(i,j,k,3) <= 2011) then
               epga = AXZ*AVG(EP_G(I,J,K),EP_G(i,j+1,k))
               IF(COUPLED) Am0 = Am0 + 0.5d0*VOL* &
                  (F_GDS(i,j,k) + F_GDS(i,j+1,k))
            ELSE
               epga = ZERO
            ENDIF

            d_n(I,J,K) = P_SCALE*epga/Am0

         ELSE
            d_n(I,J,K) = ZERO
         ENDIF

          ENDDO
        ENDDO
      ENDDO

   end subroutine calc_d_n

   subroutine calc_d_t(slo, shi, wlo, whi, lo, hi, d_t, A_m, ep_g, f_gds, flag, dx, dy, dz)

      integer     , intent(in   ) :: slo(3),shi(3),wlo(3),whi(3),lo(3),hi(3)

      ! Pressure correction
      real(c_real), intent(  out) :: d_t&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ):: A_m&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3), -3:3)

      real(c_real), intent(in   ):: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer     , intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dx, dy, dz

      integer      :: i,j,k
      real(c_real) :: axy, vol
      real(c_real) :: Am0, epga
      logical      :: coupled

      COUPLED = (DES_CONTINUUM_COUPLED .AND. .NOT.DES_ONEWAY_COUPLED)

      axy = dx*dy
      vol = dx*dy*dz
 
      DO K = lo(3)-1, hi(3)
        DO J = lo(2), hi(2)
          DO I = lo(1), hi(1)

         Am0 = -A_m(I,J,K,0)

         IF(abs(Am0) > SMALL_NUMBER) THEN

            if(flag(i,j,k,4) >= 2000 .and. &
               flag(i,j,k,4) <= 2011) then
               epga = axy*AVG(EP_G(I,J,K),EP_G(i,j,k+1))
               IF(COUPLED) Am0 = Am0 + 0.5d0*VOL* &
                  (F_GDS(I,J,K) + F_GDS(i,j,k+1))
            ELSE
               epga = ZERO
            ENDIF

            d_t(I,J,K) = P_SCALE*epga/Am0

         ELSE
            d_t(I,J,K) = ZERO
         ENDIF

          ENDDO
        ENDDO
      ENDDO

   end subroutine calc_d_t

   end module calc_d_mod
