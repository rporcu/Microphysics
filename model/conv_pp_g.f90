module conv_pp_g_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_Pp_g                                               C
!  Purpose: Determine convection terms for Pressure correction         C
!           equation.                                                  C
!                                                                      C
!  Notes: The off-diagonal coefficients calculated here must be        C
!         positive. The center coefficient and the source vector are   C
!         negative.                                                    C
!         Multiplication with factors d_e, d_n, and d_t are carried    C
!         out in source_pp_g.  Constant pressure boundaries are        C
!         handled by holding the Pp_g at the boundaries zero.  For     C
!         specified mass flow boundaries (part of) a's are calculated  C
!         here since b is calculated from a's in source_pp_g.  After   C
!         calculating b, a's are multiplied by d and at the flow       C
!         boundaries are set to zero.                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: to use face densities calculated in CONV_ROP               C
!  Author: M. Syamlal                                 Date: 1-JUN-05   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_PP_G(A_M, rop_ge, rop_gn, rop_gt, flag)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar   , only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use matrix   , only: e, w, s, n, t, b
      USE geometry , only: axy, ayz, axz
      USE functions, only: iplus, iminus, jminus, jplus, kminus, kplus

      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)

      real(c_real), INTENT(IN   ) :: rop_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: rop_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(IN   ) :: rop_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      integer ::  i,j,k
! local value of A_m
      real(c_real) :: am
!-----------------------------------------------

! Calculate convection fluxes through each of the faces
      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2
         IF (flag(i,j,k,1)==1) THEN

! East face (i+1/2, j, k)
            AM = ROP_GE(I,J,K)*AYZ
            A_M(I,J,K,E) = AM
            A_M(iplus(i,j,k),j,k,W) = AM

! North face (i, j+1/2, k)
            AM = ROP_GN(I,J,K)*AXZ
            A_M(I,J,K,N) = AM
            A_M(I,jplus(i,j,k),K,S) = AM

! Top face (i, j, k+1/2)
            AM = ROP_GT(I,J,K)*AXY
            A_M(I,J,K,T) = AM
            A_M(I,J,kplus(i,j,k),B) = AM

! West face (i-1/2, j, k)
            IF (flag(iminus(i,j,k),j,k,1) /= 1) THEN
               AM = ROP_GE(iminus(i,j,k),j,k)*AYZ
               A_M(I,J,K,W) = AM
            ENDIF

! South face (i, j-1/2, k)
            IF (flag(i,jminus(i,j,k),k,1) /=1 ) THEN
               AM = ROP_GN(i,jminus(i,j,k),k)*AXZ
               A_M(I,J,K,S) = AM
            ENDIF

! Bottom face (i, j, k-1/2)
            IF (flag(i,j,kminus(i,j,k),1)/=1) THEN
               AM = ROP_GT(i,j,kminus(i,j,k))*AXY
               A_M(I,J,K,B) = AM
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE CONV_PP_G

end module conv_pp_g_module
