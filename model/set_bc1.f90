module set_bc1_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

  contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_BC1                                                 C
!  Purpose: Set transient flow boundary conditions                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
     SUBROUTINE SET_BC1(time, dt, slo, shi, p_g, ep_g, ro_g, rop_g, &
                        u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                        flux_ge, flux_gn, flux_gt, flag, dx, dy, dz) &
       bind(C, name="set_bc1")

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_defined, bc_type
      USE param , only : dimension_bc

      use set_outflow_module, only: set_outflow

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)
      real(c_real),   intent(in   ) :: dt, time, dx, dy, dz

      integer(c_int), intent(in   ) :: flag&
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
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):shi(3))
      real(c_real), intent(inout) :: flux_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) :: flux_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) :: flux_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):shi(3))

! Local variables
!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER :: L


! Set the boundary conditions
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

            SELECT CASE(TRIM(BC_TYPE(L)))
            CASE ('P_OUTFLOW')
               CALL set_outflow(L,slo,shi,p_g,ep_g,ro_g,rop_g, &
                                u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                                flux_ge,flux_gn,flux_gt,flag)
               CALL SET_BC1_REPORT_OUTFLOW(L, time, dt, slo, shi, &
                                           u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                                           rop_g,ep_g,dx,dy,dz)
            CASE ('MASS_OUTFLOW')
               CALL set_outflow(L,slo,shi,p_g,ep_g,ro_g,rop_g, &
                                u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                                flux_ge,flux_gn,flux_gt,flag)
               CALL SET_BC1_ADJUST_OUTFLOW(L, time, dt, slo, shi, &
                                           u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                                           rop_g,ep_g,dx,dy,dz)
            CASE ('MASS_INFLOW')
            CASE ('P_INFLOW')
               CALL set_outflow(L,slo,shi,p_g,ep_g,ro_g,rop_g, &
                                u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                                flux_ge,flux_gn,flux_gt,flag)
            CASE ('OUTFLOW')
               CALL set_outflow(L,slo,shi,p_g,ep_g,ro_g,rop_g, &
                                u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                                flux_ge,flux_gn,flux_gt,flag)
               CALL SET_BC1_REPORT_OUTFLOW(L,time, dt, slo, shi, &
                                u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                                rop_g,ep_g,dx,dy,dz)
            END SELECT
         ENDIF   ! end if (bc_defined(l))
      ENDDO    ! end do loop (l=1,dimension_bc)

      RETURN
      END SUBROUTINE SET_BC1



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc1_report_outflow                                  C
!  Purpose: print out outflow conditions                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1_REPORT_OUTFLOW(BCV, time, dt, slo, shi,&
                                        u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                                        rop_g, ep_g, dx, dy, dz)

      use bc, only: bc_dt_0, bc_time
      use bc, only: bc_mout_g
      use bc, only: bc_out_n
      use bc, only: bc_vout_g
      use calc_outflow_module, only: calc_outflow
      use funits, only: dmp_log, unit_log
      use param1, only: is_undefined, zero
      use run, only: tstop

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)

      real(c_real), intent(in) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in) :: dt, time, dx, dy, dz

! Dummy arguments
!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER, intent(IN) :: BCV

      IF (IS_UNDEFINED(BC_DT_0(BCV))) RETURN

      CALL CALC_OUTFLOW(BCV,slo,shi, &
                        u_g,ulo,uhi,v_g,vlo,vhi,w_g,wlo,whi,&
                        rop_g,ep_g,dx,dy,dz)

! Calculate and accumulate the actual mass and volume outflow
      IF (TIME + 0.1d0*DT>=BC_TIME(BCV) .OR. &
          TIME+0.1d0*DT>=TSTOP) THEN
         BC_TIME(BCV) = TIME + BC_DT_0(BCV)

! Average and print out the flow rates
         BC_MOUT_G(BCV) = ABS(BC_MOUT_G(BCV))/BC_OUT_N(BCV)
         BC_VOUT_G(BCV) = ABS(BC_VOUT_G(BCV))/BC_OUT_N(BCV)
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) BCV, TIME
         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) BC_MOUT_G(BCV), &
            BC_VOUT_G(BCV)
         BC_MOUT_G(BCV) = ZERO
         BC_VOUT_G(BCV) = ZERO
         BC_OUT_N(BCV) = 0
      ENDIF

 1000 FORMAT(/,1X,'Average outflow rates at BC No. ',I2,'  At Time = ',G12.5)
 1100 FORMAT(3X,'Gas : Mass flow = ',G12.5,'     Volumetric flow = ',G12.5)

      RETURN
      END SUBROUTINE SET_BC1_REPORT_OUTFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc1_adjust_outflow                                  C
!  Purpose: Adjust velocities to get specified mass or volumetric      C
!  flow rate based on average outflow rate                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1_ADJUST_OUTFLOW(BCV, time, dt, slo, shi, &
                                        u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                                        rop_g, ep_g, dx, dy, dz)

      use bc, only: bc_dt_0, bc_time
      use bc, only: bc_i_w, bc_i_e
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_massflow_g
      use bc, only: bc_mout_g
      use bc, only: bc_out_n
      use bc, only: bc_plane
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: bc_volflow_g
      use bc, only: bc_vout_g
      use calc_outflow_module, only: calc_outflow
      use functions, only: iminus, jminus, kminus
      use funits, only: dmp_log, unit_log
      use param1, only: is_defined, zero, small_number
      use run, only: tstop

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3), vlo(3), vhi(3), wlo(3), whi(3)

! index for boundary condition
      integer, intent(in) :: bcv
      real(c_real), intent(in) :: dt, time, dx, dy, dz

      real(c_real), intent(in   ) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))


      real(c_real), intent(inout) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K
!---------------------------------------------------------------------//

      CALL CALC_OUTFLOW(BCV,slo,shi, &
                        u_g,ulo,uhi,v_g,vlo,vhi,w_g,wlo,whi,&
                        rop_g,ep_g,dx,dy,dz)

! Calculate and accumulate the actual mass and volume outflow
      IF (TIME + 0.1d0*DT>=BC_TIME(BCV) .OR. &
          TIME+0.1d0*DT>=TSTOP) THEN
         BC_TIME(BCV) = TIME + BC_DT_0(BCV)

! Average and print out the flow rates
         BC_MOUT_G(BCV) = ABS(BC_MOUT_G(BCV))/BC_OUT_N(BCV)
         BC_VOUT_G(BCV) = ABS(BC_VOUT_G(BCV))/BC_OUT_N(BCV)
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) BCV, TIME
         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) BC_MOUT_G(BCV), &
            BC_VOUT_G(BCV)
         BC_OUT_N(BCV) = 0

! Now that we know the mass and volume outflow update the bc velocities
! (gas phase)
         IF (IS_DEFINED(BC_MASSFLOW_G(BCV))) THEN
            IF (BC_MOUT_G(BCV) > SMALL_NUMBER) THEN
               SELECT CASE (TRIM(BC_PLANE(BCV)))
               CASE ('W', 'E')
                  BC_U_G(BCV) = BC_U_G(BCV)*BC_MASSFLOW_G(BCV)/&
                     BC_MOUT_G(BCV)
               CASE ('S', 'N')
                  BC_V_G(BCV) = BC_V_G(BCV)*BC_MASSFLOW_G(BCV)/&
                     BC_MOUT_G(BCV)
               CASE ('B', 'T')
                  BC_W_G(BCV) = BC_W_G(BCV)*BC_MASSFLOW_G(BCV)/&
                     BC_MOUT_G(BCV)
               END SELECT
            ENDIF
         ELSEIF (IS_DEFINED(BC_VOLFLOW_G(BCV))) THEN
            IF (BC_VOUT_G(BCV) > SMALL_NUMBER) THEN
               SELECT CASE (TRIM(BC_PLANE(BCV)))
               CASE ('W', 'E')
                  BC_U_G(BCV) = BC_U_G(BCV)*BC_VOLFLOW_G(BCV)/&
                     BC_VOUT_G(BCV)
               CASE ('S', 'N')
                  BC_V_G(BCV) = BC_V_G(BCV)*BC_VOLFLOW_G(BCV)/&
                     BC_VOUT_G(BCV)
               CASE ('B', 'T')
                  BC_W_G(BCV) = BC_W_G(BCV)*BC_VOLFLOW_G(BCV)/&
                     BC_VOUT_G(BCV)
               END SELECT
            ENDIF
         ENDIF
! zero out counter for new cycle
         BC_MOUT_G(BCV) = zero
         BC_VOUT_G(BCV) = zero

! Apply updated boundary velocities - Define the field variables at the
! boundaries according to user specifications with modifications from
! the above calculations.
! If the boundary plane is W, S, or B (i.e., the fluid cell is on the
! west, south or bottom of the boundary cell) then define the velocity
! of the adjacent fluid cell according to the boundary velocity rather
! than the velocity of the boundary cell.
! Why not set the velocity in the boundary cell itself?  Based on the
! momentum bc routine it should not really matter for MO.
         DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO J = BC_J_S(BCV), BC_J_N(BCV)
         DO I = BC_I_W(BCV), BC_I_E(BCV)
            SELECT CASE (TRIM(BC_PLANE(BCV)))
            CASE ('W'); U_G(iminus(i,j,k),j,k) = BC_U_G(BCV)
            CASE ('E'); U_G(I,J,K) = BC_U_G(BCV)
            CASE ('S'); V_G(i,jminus(i,j,k),k) = BC_V_G(BCV)
            CASE ('N'); V_G(I,J,K) = BC_V_G(BCV)
            CASE ('B'); W_G(i,j,kminus(i,j,k)) = BC_W_G(BCV)
            CASE ('T'); W_G(I,J,K) = BC_W_G(BCV)
            END SELECT

         ENDDO
         ENDDO
         ENDDO
      ENDIF   ! if time to update outflow condition

 1000 FORMAT(/,1X,'Average outflow rates at BC No. ',I2,'  At Time = ',G12.5)
 1100 FORMAT(3X,'Gas : Mass flow = ',G12.5,'     Volumetric flow = ',G12.5)

      RETURN
      END SUBROUTINE SET_BC1_ADJUST_OUTFLOW

end module set_bc1_module
