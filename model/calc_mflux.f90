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
      SUBROUTINE CALC_MFLUX (u, v, w, rop_e, rop_n, rop_t, &
         flux_e, flux_n, flux_t, flag) bind(C, name="calc_mflux")

      USE compar, only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE functions, only: iminus, jminus, kminus
      USE geometry, only: ayz, axz, axy

      implicit none

      real(c_real), intent(in   ) :: u&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: v&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: w&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: rop_e&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: rop_n&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: rop_t&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(inout) :: flux_e&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(inout) :: flux_n&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(inout) :: flux_t&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      integer(c_int), intent(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)

! Local variables
!---------------------------------------------------------------------//
! Indices
      integer :: i,j,k
!---------------------------------------------------------------------//

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
