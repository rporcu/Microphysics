module set_outflow_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_OUTFLOW                                             C
!  Purpose: Set specified outflow bc for pressure outflow,             C
!  mass outflow, outflow and now also pressure inflow bc               C
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
      SUBROUTINE SET_OUTFLOW(BCV,slo,shi,p_g,ep_g,ro_g,rop_g,u_g,v_g,w_g, &
                              flux_ge, flux_gn, flux_gt, flag)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e

      use functions, only: iminus,iplus,jminus,jplus,kminus,kplus

      use param1, only: is_undefined, zero
      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

! Boundary condition number
      integer, intent(in   ) :: bcv
      integer, intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_ge&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_gn&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_gt&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K
      INTEGER :: ip,im,jp,jm,kp,km
! local value for normal component of gas and solids velocity defined
! such that
      real(c_real) :: rvel_g
!---------------------------------------------------------------------//

! Loop over the range of boundary cells
      DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO J = BC_J_S(BCV), BC_J_N(BCV)
            DO I = BC_I_W(BCV), BC_I_E(BCV)

               im = iminus(i,j,k)
               ip = iplus(i,j,k)
               jm = jminus(i,j,k)
               jp = jplus(i,j,k)
               km = kminus(i,j,k)
               kp = kplus(i,j,k)

! Fluid cell at West
! --------------------------------------------------------------------//
               IF (flag(im,j,k,1)==1) THEN
                  RVEL_G = U_G(im,j,k)

                  CALL SET_OUTFLOW_MISC(BCV,slo,shi,I,J,K,im,j,k,p_g,ro_g)
                  CALL SET_OUTFLOW_EP(BCV,slo,shi,ro_g,rop_g,ep_g,I,J,K,im,j,k,RVEL_G)

! Set the boundary cell value of the normal component of velocity
! according to the value in the adjacent fluid cell. Note the value
! of the boundary velocity is a scaled version of the value of the
! adjacent fluid velocity based on the concentration ratio of the fluid
! cell to the boundary cell.
! - For the gas phase, this ratio is most likely 1 except for
!   compressible cases with a PO/PI boundary where P_g of the boundary
!   is set and may differ from the value of the adjacent fluid cell.
! - For the solids phase this seems unnecessary..? Differences may arise
                  IF (ROP_G(I,J,K) > ZERO) THEN
                     U_G(I,J,K) = ROP_G(im,j,k)*U_G(im,j,k)/ROP_G(I,J,K)
                  ELSE
                     U_G(I,J,K) = ZERO
                  ENDIF

! the tangential components are not explicitly handled in the boundary
! condition routines of the corresponding momentum equation
                  V_G(I,J,K) = V_G(im,j,k)
                  W_G(I,J,K) = W_G(im,j,k)

                  Flux_gE(I,J,K) = Flux_gE(im,j,k)
                  Flux_gN(I,J,K) = Flux_gN(im,j,k)
                  Flux_gT(I,J,K) = Flux_gT(im,j,k)
               ENDIF


! Fluid cell at East
! --------------------------------------------------------------------//
               IF (flag(ip,j,k,1)==1) THEN
! define normal component such that it is positive when exiting the
! domain
                  RVEL_G = -U_G(I,J,K)

                  CALL SET_OUTFLOW_MISC(BCV,slo,shi,I,J,K,ip,j,k,p_g,ro_g)
                  CALL SET_OUTFLOW_EP(BCV,slo,shi,ro_g,rop_g,ep_g,I,J,K, ip,j,k, RVEL_G)

! provide an initial value for the velocity component through the domain
! otherwise its present value (from solution of the corresponding
! momentum eqn) is kept. values for the velocity components in the off
! directions are modified (needed for PO or O boundaries but not MO or
! PI as velocities should be fully specified by this point)

                  IF (IS_UNDEFINED(U_G(I,J,K))) THEN
                     IF (ROP_G(I,J,K) > ZERO) THEN
                        U_G(I,J,K) = ROP_G(ip,j,k)*U_G(ip,j,k)/ROP_G(I,J,K)
                     ELSE
                        U_G(I,J,K) = ZERO
                     ENDIF
                  ENDIF
                  V_G(I,J,K) = V_G(ip,j,k)
                  W_G(I,J,K) = W_G(ip,j,k)

                  Flux_gE(I,J,K) = Flux_gE(ip,j,k)
                  Flux_gN(I,J,K) = Flux_gN(ip,j,k)
                  Flux_gT(I,J,K) = Flux_gT(ip,j,k)
               ENDIF


! Fluid cell at South
! --------------------------------------------------------------------//
               IF (flag(i,jm,k,1)==1) THEN
                  RVEL_G = V_G(i,jm,k)

                  CALL SET_OUTFLOW_MISC(BCV,slo,shi,I,J,K,i,jm,k,p_g,ro_g)
                  CALL SET_OUTFLOW_EP(BCV,slo,shi,ro_g,rop_g,ep_g,I,J,K, i,jm,k, RVEL_G)

                  IF (ROP_G(I,J,K) > ZERO) THEN
                     V_G(I,J,K) = ROP_G(i,jm,k)*V_G(i,jm,k)/ROP_G(I,J,K)
                  ELSE
                     V_G(I,J,K) = ZERO
                  ENDIF
                  U_G(I,J,K) = U_G(i,jm,k)
                  W_G(I,J,K) = W_G(i,jm,k)

                  Flux_gE(I,J,K) = Flux_gE(i,jm,k)
                  Flux_gN(I,J,K) = Flux_gN(i,jm,k)
                  Flux_gT(I,J,K) = Flux_gT(i,jm,k)
               ENDIF


! Fluid cell at North
! --------------------------------------------------------------------//
               IF (flag(i,jp,k,1)==1) THEN
                  RVEL_G = -V_G(I,J,K)

                  CALL SET_OUTFLOW_MISC(BCV,slo,shi,I,J,K,i,jp,k,p_g,ro_g)
                  CALL SET_OUTFLOW_EP(BCV,slo,shi,ro_g,rop_g,ep_g,I,J,K, i,jp,k, RVEL_G)

                  IF (IS_UNDEFINED(V_G(I,J,K))) THEN
                     IF (ROP_G(I,J,K) > ZERO) THEN
                        V_G(I,J,K) = ROP_G(i,jp,k)*V_G(i,jp,k)/ROP_G(I,J,K)
                     ELSE
                        V_G(I,J,K) = ZERO
                     ENDIF
                  ENDIF
                  U_G(I,J,K) = U_G(i,jp,k)
                  W_G(I,J,K) = W_G(i,jp,k)

                  Flux_gE(I,J,K) = Flux_gE(i,jp,k)
                  Flux_gN(I,J,K) = Flux_gN(i,jp,k)
                  Flux_gT(I,J,K) = Flux_gT(i,jp,k)
               ENDIF


! Fluid cell at Bottom
! --------------------------------------------------------------------//
               IF (flag(i,j,km,1)==1) THEN
                  RVEL_G = W_G(i,j,km)

                  CALL SET_OUTFLOW_MISC(BCV,slo,shi,I,J,K,i,j,km,p_g,ro_g)
                  CALL SET_OUTFLOW_EP(BCV,slo,shi,ro_g,rop_g,ep_g,I,J,K,i,j,km,RVEL_G)

                  IF (ROP_G(I,J,K) > ZERO) THEN
                     W_G(I,J,K) = ROP_G(i,j,km)*W_G(i,j,km)/ROP_G(I,J,K)
                  ELSE
                     W_G(I,J,K) = ZERO
                  ENDIF

                  U_G(I,J,K) = U_G(i,j,km)
                  V_G(I,J,K) = V_G(i,j,km)

                  Flux_gE(I,J,K) = Flux_gE(i,j,km)
                  Flux_gN(I,J,K) = Flux_gN(i,j,km)
                  Flux_gT(I,J,K) = Flux_gT(i,j,km)
               ENDIF


! Fluid cell at Top
! --------------------------------------------------------------------//
               IF (flag(i,j,kp,1)==1) THEN
                  RVEL_G = -W_G(I,J,K)

                  CALL SET_OUTFLOW_MISC(BCV,slo,shi,I,J,K,i,j,kp,p_g,ro_g)
                  CALL SET_OUTFLOW_EP(BCV,slo,shi,ro_g,rop_g,ep_g,I,J,K, i,j,kp,RVEL_G)

                  IF (IS_UNDEFINED(W_G(I,J,K))) THEN
                     IF (ROP_G(I,J,K) > ZERO) THEN
                        W_G(I,J,K) = ROP_G(i,j,kp)*W_G(i,j,kp)/ROP_G(I,J,K)
                     ELSE
                        W_G(I,J,K) = ZERO
                     ENDIF
                  ENDIF

                  U_G(I,J,K) = U_G(i,j,kp)
                  V_G(I,J,K) = V_G(i,j,kp)

                  Flux_gE(I,J,K) = Flux_gE(i,j,kp)
                  Flux_gN(I,J,K) = Flux_gN(i,j,kp)
                  Flux_gT(I,J,K) = Flux_gT(i,j,kp)
               ENDIF

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
      SUBROUTINE SET_OUTFLOW_MISC(BCV,slo,shi,I,J,K,FI,FJ,FK,p_g,ro_g)

! Global variables
!---------------------------------------------------------------------//
      use bc       , only: bc_type
      use fld_const, only: ro_g0, mw_avg
      use eos      , only: EOSG

! Global parameters
!---------------------------------------------------------------------//
      use param1, only: is_undefined
      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      ! Boundary condition number
      INTEGER, INTENT(IN) :: BCV

      ! i,j,k index for boundary cell
      INTEGER, INTENT(IN) :: I,J,K

       ! i,j,k index for adjacent fluid cell
      INTEGER, INTENT(IN) :: FI,FJ,FK

      real(c_real), INTENT(INOUT) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(INOUT) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
!---------------------------------------------------------------------//

      IF (BC_TYPE(BCV) /= 'P_OUTFLOW' .AND. &
          BC_TYPE(BCV) /= 'P_INFLOW') P_G(I,J,K) = P_G(FI,FJ,FK)

      IF (IS_UNDEFINED(RO_G0)) RO_G(I,J,K) = &
         EOSG(mw_avg,P_G(I,J,K),295.15d0)


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
      SUBROUTINE SET_OUTFLOW_EP(BCV,slo,shi,ro_g,rop_g,ep_g,I,J,K,FI,FJ,FK,RVEL_G)

! Global variables
!---------------------------------------------------------------------//
      use bc, only: bc_ep_g

! Global parameters
!---------------------------------------------------------------------//
      use param1, only: is_undefined, zero, one

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      real(c_real), INTENT(IN   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(INOUT) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(INOUT) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Dummy arguments
!---------------------------------------------------------------------//
! Boundary condition number
      INTEGER, INTENT(IN) :: BCV
! i,j,k index for boundary cell
      INTEGER, INTENT(IN) :: I,J,K
! i,j,k index for adjacent fluid cell
      INTEGER, INTENT(IN) :: FI,FJ,FK
! the gas or solids velocity in the fluid cell adjacent to the boundary
! cell dot with the outward normal of that bc plane; defines the gas or
! solids velocity component normal to the bc plane as positive when it
! is flowing into the bc cell from the fluid cell. so, for example, for
! an outflow on the eastern boundary this is the u component of velocity
! while for an outflow on the western boundary this is the -u component,
! etc.
      real(c_real), INTENT(IN) :: RVEL_G

! Local variables
!---------------------------------------------------------------------//
! sum of solids phases volume fractions
      real(c_real) :: SUM_EPs
! sum of solids phases bulk densities
      real(c_real) :: SUM_ROPS
!---------------------------------------------------------------------//

! initializing summation quantities
      SUM_ROPS = ZERO
      SUM_EPS = ZERO

! if bc_ep_g undefined, set ep_g accordingly (based on flow condition
! or based on bc_rop_s). if bc_ep_g is defined its set value will be
! maintained (from set_bc0).
      IF (IS_UNDEFINED(BC_EP_G(BCV))) EP_G(I,J,K) = ONE

! now that ep_g in the boundary cell is known, define the bulk density
! of the gas phase in the boundary cell
      ROP_G(i,j,k) = RO_G(I,J,K)*EP_G(I,J,K)

      RETURN
      END SUBROUTINE SET_OUTFLOW_EP

end module set_outflow_module
