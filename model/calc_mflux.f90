MODULE CALC_MFLUX_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MFLUX                                              C
!  Purpose: Calculate the convection fluxes. Master routine.           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_MFLUX (slo, shi, lo, hi, u, v, w, rop_e, rop_n, rop_t, &
         flux_e, flux_n, flux_t, flag, dx, dy, dz) bind(C, name="calc_mflux")

      USE compar, only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE functions, only: iminus, jminus, kminus

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      real(c_real), intent(in   ) :: u&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: rop_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in) :: dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Indices
      integer      :: i,j,k
      real(c_real) :: ayz, axz, axy
!---------------------------------------------------------------------//

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

         IF (1.eq.flag(i,j,k,1)) THEN

            ! East face (i+1/2, j, k)
            Flux_E(i,j,k) = ROP_E(i,j,k)*AYZ*U(i,j,k)

            ! West face (i-1/2, j, k)
            IF (.NOT.1.eq.flag(iminus(i,j,k),j,k,1)) then
               Flux_E(iminus(i,j,k),j,k) = ROP_E(iminus(i,j,k),j,k)*AYZ*U(iminus(i,j,k),j,k)
            ENDIF

            ! North face (i, j+1/2, k)
            Flux_N(i,j,k) = ROP_N(i,j,k)*AXZ*V(i,j,k)
            ! South face (i, j-1/2, k)
            IF (.NOT.1.eq.flag(i,jminus(i,j,k),k,1)) then
              Flux_N(i,jminus(i,j,k),k) = ROP_N(i,jminus(i,j,k),k)*AXZ*V(i,jminus(i,j,k),k)
            ENDIF


            ! Top face (i, j, k+1/2)
            Flux_T(i,j,k) = ROP_T(i,j,k)*AXY*W(i,j,k)

            ! Bottom face (i, j, k-1/2)
            IF (.NOT.1.eq.flag(i,j,kminus(i,j,k),1)) then
               Flux_T(i,j,kminus(i,j,k)) = ROP_T(i,j,kminus(i,j,k))*AXY*W(i,j,kminus(i,j,k))
            ENDIF

         ENDIF
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE CALC_MFLUX
END MODULE CALC_MFLUX_MODULE
