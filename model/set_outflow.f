!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_OUTFLOW                                             C
!  Purpose: Set specified outflow bc for pressure outflow,             C
!  mass outflow, outflow and now also pressure inflow bc               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!                                                                      C
!  Comments:                                                           C
!  If the outflow boundary is on the W, S or B side of the domain and  C
!  the component of velocity through the plane is defined then this    C
!  routine will NOT modify it (i.e., its value from the momentum       C
!  solver is maintained).                                              C
!                                                                      C
!  In general it would seem this routine does little more than what    C
!  is already done within the momentum routines in terms of velocity.  C
!  The normal component is either 1) untouched if the outflow is on    C
!  the W, S, B side of the domain or 2) is set to value of the         C
!  adjacent fluid cell if it is on the E, N, T side (similar to the    C
!  respective momentum bc routines). The primary addition here is      C
!  that the tangential components of a bc cell are set to that of      C
!  the adjacent fluid cell. Note the tangential components are not     C
!  explicitly handled in the momentum _BC_ routines; instead their     C
!  values are based on solution of the momentum equation which is      C
!  replaced here                                                       C
!                                                                      C
!  Several routines are called which perform the following tasks:      C
!  set_outflow_misc - several derived quantities are set in the        C
!      boundary                                                        C
!  set_outflow_ep - the void/volume fraction and bulk densities are    C
!      set in the boundary                                             C
!  set_outflow_fluxes - convective fluxes are set in the boundary      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_OUTFLOW(BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e

      use fldvar, only: rop_g
      use fldvar, only: u_g, v_g, w_g

      use functions, only: im_of, ip_of, jm_of, jp_of, km_of, kp_of
      use functions, only: fluid_at
      use functions, only: is_on_mype_plus2layers
      use functions, only: funijk
      use compar, only: dead_cell_at

      use param, only: dimension_m
      use param1, only: undefined, zero
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Boundary condition number
      INTEGER, INTENT(IN) :: BCV

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K
! index for boundary cell
      INTEGER :: IJK
! index for a fluid cell adjacent to the boundary cell
      INTEGER :: FIJK
! local value for normal component of gas and solids velocity defined
! such that
      DOUBLE PRECISION :: RVEL_G, RVEL_S(DIMENSION_M)
!---------------------------------------------------------------------//

! Loop over the range of boundary cells
      DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO J = BC_J_S(BCV), BC_J_N(BCV)
            DO I = BC_I_W(BCV), BC_I_E(BCV)
! Check if current i,j,k resides on this PE
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF(DEAD_CELL_AT(I,J,K)) CYCLE
               IJK = FUNIJK(I,J,K)

! Fluid cell at West
! --------------------------------------------------------------------//
               IF (FLUID_AT(IM_OF(IJK))) THEN
                  FIJK = IM_OF(IJK)
                  RVEL_G = U_G(FIJK)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)

! Set the boundary cell value of the normal component of velocity
! according to the value in the adjacent fluid cell. Note the value
! of the boundary velocity is a scaled version of the value of the
! adjacent fluid velocity based on the concentration ratio of the fluid
! cell to the boundary cell.
! - For the gas phase, this ratio is most likely 1 except for
!   compressible cases with a PO/PI boundary where P_g of the boundary
!   is set and may differ from the value of the adjacent fluid cell.
! - For the solids phase this seems unnecessary..? Differences may arise
                  IF (ROP_G(IJK) > ZERO) THEN
                    U_G(IJK) = ROP_G(FIJK)*U_G(FIJK)/ROP_G(IJK)
                  ELSE
                     U_G(IJK) = ZERO
                  ENDIF

! the tangential components are not explicitly handled in the boundary
! condition routines of the corresponding momentum equation
                  V_G(IJK) = V_G(FIJK)
                  W_G(IJK) = W_G(FIJK)

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! end if (fluid_at(im_of(ijk)))


! Fluid cell at East
! --------------------------------------------------------------------//
               IF (FLUID_AT(IP_OF(IJK))) THEN
                  FIJK = IP_OF(IJK)
! define normal component such that it is positive when exiting the
! domain
                  RVEL_G = -U_G(IJK)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)

! provide an initial value for the velocity component through the domain
! otherwise its present value (from solution of the corresponding
! momentum eqn) is kept. values for the velocity components in the off
! directions are modified (needed for PO or O boundaries but not MO or
! PI as velocities should be fully specified by this point)
                  IF (U_G(IJK) == UNDEFINED) THEN
                     IF (ROP_G(IJK) > ZERO) THEN
                        U_G(IJK) = ROP_G(FIJK)*U_G(FIJK)/ROP_G(IJK)
                     ELSE
                        U_G(IJK) = ZERO
                     ENDIF
                  ENDIF
                  V_G(IJK) = V_G(FIJK)
                  W_G(IJK) = W_G(FIJK)

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! end if (fluid_at(ip_of(ijk)))


! Fluid cell at South
! --------------------------------------------------------------------//
               IF (FLUID_AT(JM_OF(IJK))) THEN
                  FIJK = JM_OF(IJK)
                  RVEL_G = V_G(FIJK)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  IF (ROP_G(IJK) > ZERO) THEN
                     V_G(IJK) = ROP_G(FIJK)*V_G(FIJK)/ROP_G(IJK)
                  ELSE
                     V_G(IJK) = ZERO
                  ENDIF
                  U_G(IJK) = U_G(FIJK)
                  W_G(IJK) = W_G(FIJK)

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! end if (fluid_at(jm_of(ijk)))


! Fluid cell at North
! --------------------------------------------------------------------//
               IF (FLUID_AT(JP_OF(IJK))) THEN
                  FIJK = JP_OF(IJK)
                  RVEL_G = -V_G(IJK)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  IF (V_G(IJK) == UNDEFINED) THEN
                     IF (ROP_G(IJK) > ZERO) THEN
                        V_G(IJK) = ROP_G(FIJK)*V_G(FIJK)/ROP_G(IJK)
                     ELSE
                        V_G(IJK) = ZERO
                     ENDIF
                  ENDIF
                  U_G(IJK) = U_G(FIJK)
                  W_G(IJK) = W_G(FIJK)

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! if (fluid_at(jp_of(ijk)))


! Fluid cell at Bottom
! --------------------------------------------------------------------//
               IF (FLUID_AT(KM_OF(IJK))) THEN
                  FIJK = KM_OF(IJK)
                  RVEL_G = W_G(FIJK)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  IF (ROP_G(IJK) > ZERO) THEN
                     W_G(IJK) = ROP_G(FIJK)*W_G(FIJK)/ROP_G(IJK)
                  ELSE
                     W_G(IJK) = ZERO
                  ENDIF
                  U_G(IJK) = U_G(FIJK)
                  V_G(IJK) = V_G(FIJK)

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! if (fluid_at(km_of(ijk)))


! Fluid cell at Top
! --------------------------------------------------------------------//
               IF (FLUID_AT(KP_OF(IJK))) THEN
                  FIJK = KP_OF(IJK)
                  RVEL_G = -W_G(IJK)

                  CALL SET_OUTFLOW_MISC(BCV, IJK, FIJK)
                  CALL SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)

                  IF (W_G(IJK) == UNDEFINED) THEN
                     IF (ROP_G(IJK) > ZERO) THEN
                        W_G(IJK) = ROP_G(FIJK)*W_G(FIJK)/ROP_G(IJK)
                     ELSE
                        W_G(IJK) = ZERO
                     ENDIF
                  ENDIF
                  U_G(IJK) = U_G(FIJK)
                  V_G(IJK) = V_G(FIJK)

                  CALL SET_OUTFLOW_FLUXES(IJK, FIJK)
               ENDIF   ! if (fluid_at(kp_of(ijk)))

            ENDDO   ! end do (i=i1,i2)
         ENDDO   ! end do (j=j1,j2)
      ENDDO   ! end do (k=k1,k2)

      RETURN
      END SUBROUTINE SET_OUTFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_OUTFLOW_MISC                                        C
!  Purpose: Set the value of certain variables in the specified        C
!  outflow boundary cell that would not otherwise be set according     C
!  to their value in the adjacent fluid cell.                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_OUTFLOW_MISC(BCV, IJK, FIJK)

! Global variables
!---------------------------------------------------------------------//
      use bc, only: bc_type
      use fldvar, only: p_g, ro_g
      use fldvar, only: ro_g0, mw_avg
      use eos, only: EOSG

! Global parameters
!---------------------------------------------------------------------//
      use param1, only: undefined
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Boundary condition number
      INTEGER, INTENT(IN) :: BCV
! ijk index for boundary cell
      INTEGER, INTENT(IN) :: IJK
! ijk index for adjacent fluid cell
      INTEGER, INTENT(IN) :: FIJK
!---------------------------------------------------------------------//

      IF (BC_TYPE(BCV) /= 'P_OUTFLOW' .AND. &
          BC_TYPE(BCV) /= 'P_INFLOW') P_G(IJK) = P_G(FIJK)

      IF (RO_G0 == UNDEFINED) RO_G(IJK) = &
         EOSG(MW_AVG,P_G(IJK),295.15d0)


      RETURN
      END SUBROUTINE SET_OUTFLOW_MISC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_OUTFLOW_EP                                          C
!  Purpose: Set the volume fraction/bulk density (i.e., EP_g           C
!  in the specified outflow boundary cell that would not               C
!  otherwise be set according to their value in the adjacent fluid     C
!  cell.                                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_OUTFLOW_EP(BCV, IJK, FIJK, RVEL_G, RVEL_S)

! Global variables
!---------------------------------------------------------------------//
      use bc, only: bc_ep_g
      use fldvar, only: rop_g, ro_g, ep_g
      use discretelement, only: discrete_element
      use discretelement, only: des_rop_s
      use physprop, only: mmax, ro_s0

! Global parameters
!---------------------------------------------------------------------//
      use param, only: dimension_m
      use param1, only: undefined, zero, one
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Boundary condition number
      INTEGER, INTENT(IN) :: BCV
! ijk index for boundary cell
      INTEGER, INTENT(IN) :: IJK
! ijk index for adjacent fluid cell
      INTEGER, INTENT(IN) :: FIJK
! the gas or solids velocity in the fluid cell adjacent to the boundary
! cell dot with the outward normal of that bc plane; defines the gas or
! solids velocity component normal to the bc plane as positive when it
! is flowing into the bc cell from the fluid cell. so, for example, for
! an outflow on the eastern boundary this is the u component of velocity
! while for an outflow on the western boundary this is the -u component,
! etc.
      DOUBLE PRECISION, INTENT(IN) :: RVEL_G
      DOUBLE PRECISION, INTENT(IN), DIMENSION(DIMENSION_M) :: RVEL_S

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: M
! solids volume fraction
      DOUBLE PRECISION :: EPs
! sum of solids phases volume fractions
      DOUBLE PRECISION :: SUM_EPs
! sum of solids phases bulk densities
      DOUBLE PRECISION :: SUM_ROPS
!---------------------------------------------------------------------//

! initializing summation quantities
      SUM_ROPS = ZERO
      SUM_EPS = ZERO

! this section must be skipped until after the initial setup of the
! discrete element portion of the simulation (set_bc1 is called once
! before the initial setup).
      IF (DISCRETE_ELEMENT .AND. ALLOCATED(DES_ROP_S)) THEN
         DO M = 1, MMAX
! unlike in the two fluid model, in the discrete element model it is
! possible to actually calculate the bulk density in a flow boundary
! cell. Currently, however, such calculations are not strictly enforced.
! therefore use the bulk density of the adjacent fluid cell
            DES_ROP_S(IJK,M) = DES_ROP_S(FIJK,M)
            SUM_ROPS = SUM_ROPS + DES_ROP_S(IJK,M)
            EPS = DES_ROP_S(IJK,M)/RO_S0(M)
            SUM_EPS = SUM_EPS + EPS
         ENDDO
      ENDIF

! if bc_ep_g undefined, set ep_g accordingly (based on flow condition
! or based on bc_rop_s). if bc_ep_g is defined its set value will be
! maintained (from set_bc0).
      IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE - SUM_EPS

! now that ep_g in the boundary cell is known, define the bulk density
! of the gas phase in the boundary cell
      ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)

      RETURN
      END SUBROUTINE SET_OUTFLOW_EP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Update convective fluxes....                               C
!  Set the value of the convective fluxes in the specified boundary    C
!  cell according to their value in the adjacent fluid cell.           C
!                                                                      C
!  Comment/concern:                                                    C
!  Should these be assigned in the same method as the velocity? Note   C
!  if bc_plane is W, S, B then the normal component of velocity may be C
!  assigned a zero value as opposed to value of its neighboring fluid  C
!  cell. This routine would seem to introduce some inconsistency       C
!  between velocity and flux at boundary.                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_OUTFLOW_FLUXES(IJK,FIJK)

! Global variables
!---------------------------------------------------------------------//
      use mflux, only: flux_ge, flux_gn, flux_gt

! Dummy arguments
!---------------------------------------------------------------------//
! ijk index for boundary cell
      INTEGER, INTENT(IN) :: IJK
! ijk index for adjacent fluid cell
      INTEGER, INTENT(IN) :: FIJK
!---------------------------------------------------------------------//

      Flux_gE(IJK) = Flux_gE(FIJK)
      Flux_gN(IJK) = Flux_gN(FIJK)
      Flux_gT(IJK) = Flux_gT(FIJK)

      RETURN
      END SUBROUTINE SET_OUTFLOW_FLUXES
