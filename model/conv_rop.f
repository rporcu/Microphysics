MODULE CONV_ROP_MODULE
   CONTAINS
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
      SUBROUTINE CONV_ROP(u_g, v_g, w_g, rop_g, rop_ge, rop_gn, rop_gt)

! Modules
!---------------------------------------------------------------------//
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE run, only: discretize

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_ge&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_gn&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_gt&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

!---------------------------------------------------------------------//

      IF (DISCRETIZE(1) == 0) THEN               ! 0 & 1 => first order upwinding
         CALL CONV_ROP0 (ROP_g, U_g, V_g, W_g, ROP_gE, ROP_gN, ROP_gT)
      ELSE
         CALL CONV_ROP1 (DISCRETIZE(1), &
                         rop_g, u_g, v_g, w_g, &
                         rop_ge, rop_gn, rop_gt)
      ENDIF

      RETURN
   CONTAINS

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
      USE functions, only: fluid_at
      USE functions, only: ieast, jnorth, ktop
      USE functions, only: iwest, jsouth, kbot
      USE functions, only: iminus, jminus, kminus
      USE param1, only: zero
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! macroscopic density (rho_prime)
      DOUBLE PRECISION, INTENT(IN) :: ROP(istart3:iend3,jstart3:jend3,kstart3:kend3)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: V(istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: W(istart3:iend3,jstart3:jend3,kstart3:kend3)
! Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_E(istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_N(istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_T(istart3:iend3,jstart3:jend3,kstart3:kend3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K
!---------------------------------------------------------------------//

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3

         IF (fluid_at(i,j,k)) THEN

! East face (i+1/2, j, k)
            IF (U(i,j,k) >= ZERO) THEN
               ROP_E(i,j,k) = ROP(i,j,k)
            ELSE
               ROP_E(i,j,k) = ROP(ieast(i,j,k),j,k)
            ENDIF
! West face (i-1/2, j, k)
            IF (.NOT.fluid_at(iminus(i,j,k),j,k)) THEN
               IF (U(iminus(i,j,k),j,k) >= ZERO) THEN
                  ROP_E(iminus(i,j,k),j,k) = ROP(iwest(i,j,k),j,k)
               ELSE
                  ROP_E(iminus(i,j,k),j,k) = ROP(i,j,k)
               ENDIF
            ENDIF


! North face (i, j+1/2, k)
            IF (V(i,j,k) >= ZERO) THEN
               ROP_N(i,j,k) = ROP(i,j,k)
            ELSE
               ROP_N(i,j,k) = ROP(i,jnorth(i,j,k),k)
            ENDIF
! South face (i, j-1/2, k)
            IF (.NOT.fluid_at(i,jminus(i,j,k),k)) THEN
               IF (V(i,jminus(i,j,k),k) >= ZERO) THEN
                 ROP_N(i,jminus(i,j,k),k) = ROP(i,jsouth(i,j,k),k)
               ELSE
                 ROP_N(i,jminus(i,j,k),k) = ROP(i,j,k)
               ENDIF
            ENDIF


! Top face (i, j, k+1/2)
            IF (W(i,j,k) >= ZERO) THEN
               ROP_T(i,j,k) = ROP(i,j,k)
            ELSE
               ROP_T(i,j,k) = ROP(i,j,ktop(i,j,k))
            ENDIF
! Bottom face (i, j, k-1/2)
            IF (.NOT.fluid_at(i,j,kminus(i,j,k))) THEN
               IF (W(i,j,kminus(i,j,k)) >= ZERO) THEN
                  ROP_T(i,j,kminus(i,j,k)) = ROP(i,j,kbot(i,j,k))
               ELSE
                  ROP_T(i,j,kminus(i,j,k)) = ROP(i,j,k)
               ENDIF
            ENDIF

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
      SUBROUTINE CONV_ROP1(DISC, &
                           rop, u, v, w, &
                           rop_e, rop_n, rop_t)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3
      USE functions, only: fluid_at
      USE functions, only: ieast, jnorth, ktop
      USE functions, only: iwest, jsouth, kbot
      USE functions, only: iminus, jminus, kminus
      USE param1, only: one
      USE xsi, only: calc_xsi
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Discretization scheme
      INTEGER, INTENT(IN) :: DISC
! macroscopic density (rho_prime)
      DOUBLE PRECISION, INTENT(in) :: rop(istart3:iend3,jstart3:jend3,kstart3:kend3)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: u(istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: v(istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: w(istart3:iend3,jstart3:jend3,kstart3:kend3)
! Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION, INTENT(OUT) :: rop_e(istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT) :: rop_n(istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT) :: rop_t(istart3:iend3,jstart3:jend3,kstart3:kend3)
!
! Local variables
!---------------------------------------------------------------------//
      INTEGER :: I,J,K
      Integer :: incr
      DOUBLE PRECISION, allocatable :: xsi_e(:,:,:), xsi_n(:,:,:), xsi_t(:,:,:)

!---------------------------------------------------------------------//

      allocate(xsi_e(istart3:iend3, jstart3:jend3, kstart3:kend3) )
      allocate(xsi_n(istart3:iend3, jstart3:jend3, kstart3:kend3) )
      allocate(xsi_t(istart3:iend3, jstart3:jend3, kstart3:kend3) )


! Calculate factors
      incr=0
      CALL CALC_XSI (DISC, ROP, U, V, W, XSI_E, XSI_N, XSI_T, incr)

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3

         IF (fluid_at(i,j,k)) THEN

! East face (i+1/2, j, k)
            ROP_E(i,j,k) = ((ONE-XSI_E(i,j,k))*ROP(i,j,k) + &
               XSI_E(i,j,k) *ROP(ieast(i,j,k),j,k) )
! West face (i-1/2, j, k)
            IF (.NOT.fluid_at(iminus(i,j,k),j,k)) THEN
               ROP_E(iminus(i,j,k),j,k) = &
                  ((ONE - XSI_E(iminus(i,j,k),j,k))*ROP(iwest(i,j,k),j,k) + &
                  XSI_E(iminus(i,j,k),j,k) *ROP(i,j,k) )
            ENDIF


! North face (i, j+1/2, k)
            ROP_N(i,j,k) = ((ONE-XSI_N(i,j,k))*ROP(i,j,k)+&
               XSI_N(i,j,k) *ROP(i,jnorth(i,j,k),k))
! South face (i, j-1/2, k)
            IF (.NOT.fluid_at(i,jminus(i,j,k),k)) THEN
               ROP_N(i,jminus(i,j,k),k) = &
                  ((ONE - XSI_N(i,jminus(i,j,k),k))*ROP(i,jsouth(i,j,k),k) + &
                  XSI_N(i,jminus(i,j,k),k) *ROP(i,j,k) )
            ENDIF

! Top face (i, j, k+1/2)
            ROP_T(i,j,k) = ((ONE - XSI_T(i,j,k))*ROP(i,j,k) + &
               XSI_T(i,j,k) *ROP(i,j,ktop(i,j,k)) )
! Bottom face (i, j, k-1/2)
            IF (.NOT.fluid_at(i,j,kminus(i,j,k))) THEN
               ROP_T(i,j,kminus(i,j,k)) = &
                  ((ONE - XSI_T(i,j,kminus(i,j,k)))*ROP(i,j,kbot(i,j,k)) + &
                  XSI_T(i,j,kminus(i,j,k)) *ROP(i,j,k) )
            ENDIF

         ENDIF   ! end if fluid_at
      ENDDO
      ENDDO
      ENDDO

      deallocate( xsi_e, xsi_n, xsi_t)

      RETURN
      END SUBROUTINE CONV_ROP1

   END SUBROUTINE CONV_ROP

END MODULE CONV_ROP_MODULE
