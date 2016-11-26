!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP                                                C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. Master routine.                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_ROP()

! Modules
!---------------------------------------------------------------------//
      USE fldvar, only: rop_g, u_g, v_g, w_g
      USE fldvar, only: rop_ge, rop_gn, rop_gt
      USE run, only: discretize
      IMPLICIT NONE

!---------------------------------------------------------------------//

      IF (DISCRETIZE(1) == 0) THEN               ! 0 & 1 => first order upwinding
         CALL CONV_ROP0 (ROP_g, U_g, V_g, W_g, &
                         ROP_gE, ROP_gN, ROP_gT)
      ELSE
         CALL CONV_ROP1 (DISCRETIZE(1), ROP_g, U_g, V_g, W_g, &
                         ROP_gE, ROP_gN, ROP_gT)
      ENDIF

      RETURN
      END SUBROUTINE CONV_ROP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP                                                C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. FOU routine.                                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_ROP0(ROP, U, V, W, ROP_E, ROP_N, ROP_T)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3
      USE functions, only: funijk
      USE functions, only: fluid_at
      USE functions, only: ieast, jnorth, ktop
      USE functions, only: iwest, jsouth, kbot
      USE functions, only: iminus, jminus, kminus
      USE geometry, only: do_k
      USE param, only: dimension_3
      USE param1, only: zero
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! macroscopic density (rho_prime)
      DOUBLE PRECISION, INTENT(IN) :: ROP(DIMENSION_3)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: W(DIMENSION_3)
! Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_T(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, I, J, K
      INTEGER :: IJKE, IJKN, IJKT
      INTEGER :: IJKW, IJKS, IJKB
      INTEGER :: IMJK, IJMK, IJKM
!---------------------------------------------------------------------//


      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)

         IF (fluid_at(i,j,k)) THEN
            IJKE = FUNIJK(ieast(i,j,k),j,k)
            IJKN = FUNIJK(i,jnorth(i,j,k),k)
            IJKT = FUNIJK(i,j,ktop(i,j,k))

            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)

! East face (i+1/2, j, k)
            IF (U(IJK) >= ZERO) THEN
               ROP_E(IJK) = ROP(IJK)
            ELSE
               ROP_E(IJK) = ROP(IJKE)
            ENDIF
! West face (i-1/2, j, k)
            IF (.NOT.fluid_at(iminus(i,j,k),j,k)) THEN
               IJKW = FUNIJK(iwest(i,j,k),j,k)
               IF (U(IMJK) >= ZERO) THEN
                  ROP_E(IMJK) = ROP(IJKW)
               ELSE
                  ROP_E(IMJK) = ROP(IJK)
               ENDIF
            ENDIF


! North face (i, j+1/2, k)
            IF (V(IJK) >= ZERO) THEN
               ROP_N(IJK) = ROP(IJK)
            ELSE
               ROP_N(IJK) = ROP(IJKN)
            ENDIF
! South face (i, j-1/2, k)
            IF (.NOT.fluid_at(i,jminus(i,j,k),k)) THEN
               IJKS = FUNIJK(i,jsouth(i,j,k),k)
               IF (V(IJMK) >= ZERO) THEN
                 ROP_N(IJMK) = ROP(IJKS)
               ELSE
                 ROP_N(IJMK) = ROP(IJK)
               ENDIF
            ENDIF


            IF (DO_K) THEN
               IJKM = FUNIJK(i,j,kminus(i,j,k))
! Top face (i, j, k+1/2)
               IF (W(IJK) >= ZERO) THEN
                  ROP_T(IJK) = ROP(IJK)
               ELSE
                  ROP_T(IJK) = ROP(IJKT)
               ENDIF
! Bottom face (i, j, k-1/2)
               IF (.NOT.fluid_at(i,j,kminus(i,j,k))) THEN
                  IJKB = FUNIJK(i,j,kbot(i,j,k))
                  IF (W(IJKM) >= ZERO) THEN
                     ROP_T(IJKM) = ROP(IJKB)
                  ELSE
                     ROP_T(IJKM) = ROP(IJK)
                  ENDIF
               ENDIF
            ENDIF   ! end if do_k

         ENDIF   ! end if fluid_at
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CONV_ROP0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP1                                               C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. HR routine.  Here interpolate the face value of  C
!  density.                                                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_ROP1(DISC, ROP, U, V, W, ROP_E, ROP_N, ROP_T)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3
      USE functions, only: funijk
      USE functions, only: fluid_at
      USE functions, only: ieast, jnorth, ktop
      USE functions, only: iwest, jsouth, kbot
      USE functions, only: iminus, jminus, kminus
      USE geometry, only: do_k
      USE param, only: dimension_3
      USE param1, only: one
      USE xsi, only: calc_xsi
      USE xsi_array, only: xsi_e, xsi_n, xsi_t
      USE xsi_array, only: lock_xsi_array, unlock_xsi_array
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Discretization scheme
      INTEGER, INTENT(IN) :: DISC
! macroscopic density (rho_prime)
      DOUBLE PRECISION, INTENT(IN) :: ROP(DIMENSION_3)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: W(DIMENSION_3)
! Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_T(DIMENSION_3)
!
! Local variables
!---------------------------------------------------------------------//
      INTEGER :: I,J,K,IJK, IJKE, IJKN, IJKT
      INTEGER :: IJKW, IJKS, IJKB, IMJK, IJMK, IJKM
      Integer :: incr

!---------------------------------------------------------------------//

      call lock_xsi_array

! Calculate factors
      incr=0
      CALL CALC_XSI (DISC, ROP, U, V, W, XSI_E, XSI_N, XSI_T, incr)


      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)

         IF (fluid_at(i,j,k)) THEN
            IJKE = FUNIJK(ieast(i,j,k),j,k)
            IJKN = FUNIJK(i,jnorth(i,j,k),k)
            IJKT = FUNIJK(i,j,ktop(i,j,k))

            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)

! East face (i+1/2, j, k)
            ROP_E(IJK) = ((ONE-XSI_E(IJK))*ROP(IJK)+&
                         XSI_E(IJK)*ROP(IJKE))
! West face (i-1/2, j, k)
            IF (.NOT.fluid_at(iminus(i,j,k),j,k)) THEN
               IJKW = FUNIJK(iwest(i,j,k),j,k)
               ROP_E(IMJK) = ((ONE - XSI_E(IMJK))*ROP(IJKW)+&
                             XSI_E(IMJK)*ROP(IJK))
            ENDIF


! North face (i, j+1/2, k)
            ROP_N(IJK) = ((ONE-XSI_N(IJK))*ROP(IJK)+&
                         XSI_N(IJK)*ROP(IJKN))
! South face (i, j-1/2, k)
            IF (.NOT.fluid_at(i,jminus(i,j,k),k)) THEN
               IJKS = FUNIJK(i,jsouth(i,j,k),k)
               ROP_N(IJMK) = ((ONE - XSI_N(IJMK))*ROP(IJKS)+&
                             XSI_N(IJMK)*ROP(IJK))
            ENDIF


            IF (DO_K) THEN
               IJKM = FUNIJK(i,j,kminus(i,j,k))

! Top face (i, j, k+1/2)
               ROP_T(IJK) = ((ONE - XSI_T(IJK))*ROP(IJK)+&
                            XSI_T(IJK)*ROP(IJKT))
! Bottom face (i, j, k-1/2)
               IF (.NOT.fluid_at(i,j,kminus(i,j,k))) THEN
                  IJKB = FUNIJK(i,j,kbot(i,j,k))
                  ROP_T(IJKM) = ((ONE - XSI_T(IJKM))*ROP(IJKB)+&
                                XSI_T(IJKM)*ROP(IJK))
               ENDIF
            ENDIF   ! end if do_k

         ENDIF   ! end if fluid_at
      ENDDO
      ENDDO
      ENDDO

      call unlock_xsi_array

      RETURN
      END SUBROUTINE CONV_ROP1
