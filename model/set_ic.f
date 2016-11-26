!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IC                                                  !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: This module sets all the initial conditions.               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE SET_IC

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE constant
      USE physprop
      USE ic
      USE fldvar
      USE scales
      USE compar
      USE run
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I, J, K, IJK
! Local index for initial condition
      INTEGER :: L
! Temporary variable for storing IC_EP_g
      DOUBLE PRECISION :: EPGX
! Temporary variable for storing IC_P_g
      DOUBLE PRECISION :: PGX
! Temporary variable for storing IC_U_g
      DOUBLE PRECISION :: UGX
! Temporary variable for storing IC_V_g
      DOUBLE PRECISION :: VGX
! Temporary variable for storing IC_W_g
      DOUBLE PRECISION :: WGX
!-----------------------------------------------

!  Set the initial conditions.
      DO L = 1, DIMENSION_IC
         IF (IC_DEFINED(L)) THEN

            EPGX = IC_EP_G(L)
            PGX = IC_P_G(L)
            UGX = IC_U_G(L)
            VGX = IC_V_G(L)
            WGX = IC_W_G(L)

            DO K = IC_K_B(L), IC_K_T(L)
            DO J = IC_J_S(L), IC_J_N(L)
            DO I = IC_I_W(L), IC_I_E(L)
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I,J,K)

               IF (.NOT.WALL_AT(i,j,k)) THEN
                  IF (EPGX /= UNDEFINED) EP_G(IJK) = EPGX

                  IF (IC_TYPE(L) == 'PATCH') THEN
                      IF (PGX /= UNDEFINED) P_G(IJK) = SCALE_PRESSURE(PGX)
                  ELSE
                     P_G(IJK) = merge(SCALE_PRESSURE(PGX), UNDEFINED,           &
                        PGX /= UNDEFINED)
                  ENDIF

                  IF (UGX /= UNDEFINED) U_G(IJK) = UGX
                  IF (VGX /= UNDEFINED) V_G(IJK) = VGX
                  IF (WGX /= UNDEFINED) W_G(IJK) = WGX


! end of modifications for GHD theory
               ENDIF     ! Fluid at
            ENDDO   ! over i
            ENDDO   ! over j
            ENDDO   ! over k
         ENDIF   ! if (ic_defined)
      ENDDO   ! over dimension_ic

      RETURN
      END SUBROUTINE SET_IC
