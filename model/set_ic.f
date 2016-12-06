!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IC                                                  !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: This module sets all the initial conditions.               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE SET_IC(ep_g, p_g, u_g, v_g, w_g, flag)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar, only: istart3,jstart3,iend3,jend3,kstart3,kend3
      USE ic, only: dimension_ic
      USE ic, only: ic_defined, ic_ep_g, ic_p_g, ic_u_g, ic_v_g, ic_w_g, ic_type, ic_i_w, ic_j_s, ic_k_b, ic_i_e, ic_j_n, ic_k_t
      USE scales, only: scale_pressure
      USE param1, only: undefined

      IMPLICIT NONE

      double precision, intent(inout) :: ep_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) ::  p_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) ::  u_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) ::  v_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) ::  w_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      integer, intent(in) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,0:4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I, J, K
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

               IF (flag(i,j,k,1)<100) THEN
                  IF (EPGX /= UNDEFINED) EP_G(I,J,K) = EPGX

                  IF (IC_TYPE(L) == 'PATCH') THEN
                      IF (PGX /= UNDEFINED) P_G(I,J,K) = SCALE_PRESSURE(PGX)
                  ELSE
                     P_G(I,J,K) = merge(SCALE_PRESSURE(PGX), UNDEFINED,           &
                        PGX /= UNDEFINED)
                  ENDIF

                  IF (UGX /= UNDEFINED) U_G(I,J,K) = UGX
                  IF (VGX /= UNDEFINED) V_G(I,J,K) = VGX
                  IF (WGX /= UNDEFINED) W_G(I,J,K) = WGX


! end of modifications for GHD theory
               ENDIF     ! Fluid at
            ENDDO   ! over i
            ENDDO   ! over j
            ENDDO   ! over k
         ENDIF   ! if (ic_defined)
      ENDDO   ! over dimension_ic

      RETURN
      END SUBROUTINE SET_IC
