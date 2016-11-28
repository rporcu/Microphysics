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

      SUBROUTINE CONV_PP_G(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar, only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE param
      USE param1
      USE fldvar
      USE run
      use matrix, only: e, w, s, n, t, b
      USE physprop
      USE geometry
      USE functions
      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IPJK, IJPK, IJKP
      INTEGER :: IMJK, IJMK, IJKM
      INTEGER ::  I, J, K
! local value of A_m
      DOUBLE PRECISION :: am
!-----------------------------------------------


! Calculate convection fluxes through each of the faces
      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2
         IJK = FUNIJK(i,j,k)
         IF (fluid_at(i,j,k)) THEN
            IPJK = FUNIJK(iplus(i,j,k),j,k)
            IJPK = FUNIJK(i,jplus(i,j,k),k)
            IJKP = FUNIJK(i,j,kplus(i,j,k))

! East face (i+1/2, j, k)
            AM = ROP_GE(I,J,K)*AYZ
            A_M(I,J,K,E) = AM
            A_M(iplus(i,j,k),j,k,W) = AM

! North face (i, j+1/2, k)
            AM = ROP_GN(I,J,K)*AXZ
            A_M(I,J,K,N) = AM
            A_M(I,jplus(i,j,k),K,S) = AM

! Top face (i, j, k+1/2)
            IF (DO_K) THEN
               AM = ROP_GT(I,J,K)*AXY
               A_M(I,J,K,T) = AM
               A_M(I,J,kplus(i,j,k),B) = AM
            ENDIF

! West face (i-1/2, j, k)
            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IF (.NOT.fluid_at(iminus(i,j,k),j,k)) THEN
               AM = ROP_GE(iminus(i,j,k),j,k)*AYZ
               A_M(I,J,K,W) = AM
            ENDIF

! South face (i, j-1/2, k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IF (.NOT.fluid_at(i,jminus(i,j,k),k)) THEN
               AM = ROP_GN(i,jminus(i,j,k),k)*AXZ
               A_M(I,J,K,S) = AM
            ENDIF

! Bottom face (i, j, k-1/2)
            IF (DO_K) THEN
               IJKM = FUNIJK(i,j,kminus(i,j,k))
               IF (.NOT.fluid_at(i,j,kminus(i,j,k))) THEN
                  AM = ROP_GT(i,j,kminus(i,j,k))*AXY
                  A_M(I,J,K,B) = AM
               ENDIF
            ENDIF
         ENDIF   ! end if (fluid_at(i,j,k))
      ENDDO
      ENDDO
      ENDDO


      RETURN
      END SUBROUTINE CONV_PP_G
