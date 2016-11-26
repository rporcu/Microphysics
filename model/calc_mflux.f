!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MFLUX                                              C
!  Purpose: Calculate the convection fluxes. Master routine.           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_MFLUX()

! Modules
!---------------------------------------------------------------------//
      USE fldvar, only: u_g, v_g, w_g
      USE fldvar, only: rop_ge, rop_gn, rop_gt
      USE fldvar, only: flux_ge, flux_gn, flux_gt

      IMPLICIT NONE
!---------------------------------------------------------------------//
         CALL CALC_MFLUX0 (U_g, V_g, W_g, ROP_gE, ROP_gN, ROP_gT, &
                           Flux_gE, Flux_gN, Flux_gT)

      RETURN
      END SUBROUTINE CALC_MFLUX

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine:: CALC_MFLUX0                                            C
!  Purpose: Calculate the convection fluxes.                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_MFLUX0(U, V, W, ROP_E, ROP_N, ROP_T,&
                             Flux_E, Flux_N, Flux_T)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE functions, only: funijk, fluid_at
      USE functions, only: iminus, jminus, kminus
      USE geometry, only: do_k
      USE geometry, only: ayz, axz, axy
      USE param, only: dimension_3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: W(DIMENSION_3)
! Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION, INTENT(IN) :: ROP_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: ROP_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: ROP_T(DIMENSION_3)
! Convective mass fluxes
      DOUBLE PRECISION, INTENT(OUT) :: Flux_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: Flux_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: Flux_T(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, IMJK, IJMK, IJKM, I, J, K
!---------------------------------------------------------------------//

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

         IJK = FUNIJK(i,j,k)
         IF (fluid_at(i,j,k)) THEN

            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)

! East face (i+1/2, j, k)
            Flux_E(IJK) = ROP_E(IJK)*AYZ*U(IJK)
! West face (i-1/2, j, k)
            IF (.NOT.fluid_at(iminus(i,j,k),j,k)) then
               Flux_E(IMJK) = ROP_E(IMJK)*AYZ*U(IMJK)
            ENDIF

! North face (i, j+1/2, k)
            Flux_N(IJK) = ROP_N(IJK)*AXZ*V(IJK)
! South face (i, j-1/2, k)
            IF (.NOT.fluid_at(i,jminus(i,j,k),k)) then
              Flux_N(IJMK) = ROP_N(IJMK)*AXZ*V(IJMK)
            ENDIF

            IF (DO_K) THEN
               IJKM = FUNIJK(i,j,kminus(i,j,k))
! Top face (i, j, k+1/2)
               Flux_T(IJK) = ROP_T(IJK)*AXY*W(IJK)
! Bottom face (i, j, k-1/2)
               IF (.NOT.fluid_at(i,j,kminus(i,j,k))) then
                 Flux_T(IJKM) = ROP_T(IJKM)*AXY*W(IJKM)
               ENDIF
            ENDIF   ! end if do_k
         ENDIF   ! end if fluid_at
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CALC_MFLUX0
