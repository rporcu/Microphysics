!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_GS_DES1                                            !
!  Author: J.Musser                                   Date: 21-NOV-14  !
!                                                                      !
!  Purpose: This routine is called from the DISCRETE side to calculate !
!  the gas-based forces acting on each particle using interpolated     !
!  values for gas velocity, gas pressure, and gas volume fraction.     !
!                                                                      !
!  Notes:                                                              !
!                                                                      !
!   F_gp is obtained from des_drag_gp subroutine is given as:          !
!    F_GP = beta*VOL_P/EP_s where VOL_P is the particle volume.        !
!                                                                      !
!  The drag force on each particle is equal to:                        !
!    D_FORCE = beta*VOL_P/EP_s*(Ug - Us) = F_GP *(Ug - Us)             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_GS_DES1

! Flag: The fluid and discrete solids are explicitly coupled.
      use discretelement, only: DES_EXPLICITLY_COUPLED
! Gas phase volume fraction
      use fldvar, only: EP_G
! Gas phase velocities
      use fldvar, only: U_G, V_G, W_G
! Size of particle arrays on this processor.
      use discretelement, only: MAX_PIP
! Flag to use interpolation
      use particle_filter, only: DES_INTERP_ON
! Interpolation cells and weights
      use particle_filter, only: FILTER_CELL, FILTER_WEIGHT
! IJK of fluid cell containing particles center
      use discretelement, only: PIJK
! Drag force on each particle
      use discretelement, only: F_GP
! Particle velocity
      use discretelement, only: DES_VEL_NEW
! Total forces acting on particle
      use discretelement, only: FC
! Gas pressure force by fluid cell
      use discretelement, only: P_FORCE
! Particle volume.
      use discretelement, only: PVOL
! Particle drag force
      use discretelement, only: DRAG_FC
! Flag for 3D simulations.
      use geometry, only: DO_K
! Function to deterine if a cell contains fluid.
      use functions, only: FLUID_AT

      use functions, only: is_normal

! Global Parameters:
!---------------------------------------------------------------------//
! Size of local fluid arrays
      use param, only: DIMENSION_3
! Double precision values.
      use param1, only: ZERO

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK, LC
! Interpolation weight
      DOUBLE PRECISION :: WEIGHT
! Interpolated gas phase quanties.
      DOUBLE PRECISION :: lEPg, VELFP(3), lPF(3)
! Drag force acting on each particle.
      DOUBLE PRECISION :: D_FORCE(3)
! Loop bound for filter
      INTEGER :: LP_BND
! Cell-center gas velocities.
      DOUBLE PRECISION, ALLOCATABLE :: UGC(:)
      DOUBLE PRECISION, ALLOCATABLE :: VGC(:)
      DOUBLE PRECISION, ALLOCATABLE :: WGC(:)
!......................................................................!

      allocate( UGC(DIMENSION_3) )
      allocate( VGC(DIMENSION_3) )
      allocate( WGC(DIMENSION_3) )

! Loop bounds for interpolation.
      LP_BND = merge(27,9,DO_K)

! Calculate the cell center gas velocities.
      CALL CALC_CELL_CENTER_GAS_VEL(U_G, V_G, W_G, UGC, VGC, WGC)

! Calculate the gas phase forces acting on each particle.

!$omp parallel default(none) private(np,lepg,velfp,ijk,weight,lpf,d_force)    &
!$omp          shared(max_pip,des_interp_on,lp_bnd,filter_cell,filter_weight, &
!$omp          ep_g,pijk,des_vel_new,f_gp,ugc,vgc,wgc,p_force,          &
!$omp          des_explicitly_coupled,drag_fc,fc,pvol)
!$omp do
      DO NP=1,MAX_PIP
         IF(.NOT.IS_NORMAL(NP)) CYCLE
! Avoid drag calculations in cells without fluid (cut-cell)
         IF(.NOT.FLUID_AT(PIJK(NP,4))) CYCLE

         lEPG = ZERO
         VELFP = ZERO
         lPF = ZERO

! Calculate the gas volume fraction, velocity, and pressure force at
! the particle's position.
         IF(DES_INTERP_ON) THEN
            DO LC=1,LP_BND
               IJK = FILTER_CELL(LC,NP)
               WEIGHT = FILTER_WEIGHT(LC,NP)
! Gas phase volume fraction.
               lEPG = lEPG + EP_G(IJK)*WEIGHT
! Gas phase velocity.
               VELFP(1) = VELFP(1) + UGC(IJK)*WEIGHT
               VELFP(2) = VELFP(2) + VGC(IJK)*WEIGHT
               VELFP(3) = VELFP(3) + WGC(IJK)*WEIGHT
! Gas pressure force.
               lPF = lPF + P_FORCE(:,IJK)*WEIGHT
            ENDDO
         ELSE
            IJK = PIJK(NP,4)
            lEPG = EP_G(IJK)
            VELFP(1) = UGC(IJK)
            VELFP(2) = VGC(IJK)
            VELFP(3) = WGC(IJK)
            lPF = P_FORCE(:,IJK)
         ENDIF

! For explicit coupling, use the drag coefficient calculated for the
! gas phase drag calculations.
         IF(DES_EXPLICITLY_COUPLED) THEN

            DRAG_FC(NP,:) = F_GP(NP)*(VELFP - DES_VEL_NEW(NP,:))

         ELSE

! Calculate the drag coefficient.
            CALL DES_DRAG_GP(NP, DES_VEL_NEW(NP,:), VELFP, lEPg)

! Calculate the gas-solids drag force on the particle
            D_FORCE = F_GP(NP)*(VELFP - DES_VEL_NEW(NP,:))

! Update the contact forces (FC) on the particle to include gas
! pressure and gas-solids drag
            FC(NP,:) = FC(NP,:) + D_FORCE(:)

            FC(NP,:) = FC(NP,:) + lPF*PVOL(NP)

         ENDIF

      ENDDO
!$omp end parallel

      deallocate( UGC, VGC, WGC )

      RETURN
      END SUBROUTINE DRAG_GS_DES1

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_GS_GAS1                                            !
!  Author: J.Musser                                   Date: 21-NOV-14  !
!                                                                      !
!                                                                      !
!  Purpose: This routine is called from the CONTINUUM. It calculates   !
!  the scalar cell center drag force acting on the fluid using         !
!  interpolated values for the gas velocity and volume fraction. The   !
!  The resulting sources are interpolated back to the fluid grid.      !
!                                                                      !
!  NOTE: The loop over particles includes ghost particles so that MPI  !
!  communications are needed to distribute overlapping force between   !
!  neighboring grid cells. This is possible because only cells "owned" !
!  by the current process will have non-zero weights.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_GS_GAS1

! Flag: The fluid and discrete solids are explicitly coupled.
      use discretelement, only: DES_EXPLICITLY_COUPLED
! Gas phase volume fraction
      use fldvar, only: EP_G
! Gas phase velocities
      use fldvar, only: U_G, V_G, W_G
      use fldvar, only: U_GO, V_GO, W_GO
! Size of particle array on this process.
      use discretelement, only: MAX_PIP
! Flag to use interpolation
      use particle_filter, only: DES_INTERP_ON
! Interpolation cells and weights
      use particle_filter, only: FILTER_CELL, FILTER_WEIGHT
! IJK of fluid cell containing particles center
      use discretelement, only: PIJK
! Drag force on each particle
      use discretelement, only: F_GP
! Particle velocity
      use discretelement, only: DES_VEL_NEW
! Contribution to gas momentum equation due to drag
      use discretelement, only: DRAG_BM
! Scalar cell center total drag force
      use discretelement, only: F_GDS
! Volume of scalar cell.
      use geometry, only: VOL
! Flag for 3D simulatoins.
      use geometry, only: DO_K
! MPI wrapper for halo exchange.
      use sendrecv, only: SEND_RECV

      use functions, only: IS_NONEXISTENT, IS_ENTERING, &
         IS_ENTERING_GHOST, IS_EXITING, IS_EXITING_GHOST

! Global Parameters:
!---------------------------------------------------------------------//
! Size of local fluid arrays
      use param, only: DIMENSION_3
! Double precision values.
      use param1, only: ZERO, ONE

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK, LC
! Interpolation weight
      DOUBLE PRECISION :: WEIGHT
! Interpolated gas phase quanties.
      DOUBLE PRECISION :: lEPg, VELFP(3)
! Loop bound for filter
      INTEGER :: LP_BND
! Drag force (intermediate calculation)
      DOUBLE PRECISION :: lFORCE
! Drag sources for fluid (intermediate calculation)
      DOUBLE PRECISION :: lDRAG_BM(3)
! Cell-center gas velocities.
      DOUBLE PRECISION, ALLOCATABLE :: UGC(:)
      DOUBLE PRECISION, ALLOCATABLE :: VGC(:)
      DOUBLE PRECISION, ALLOCATABLE :: WGC(:)
!......................................................................!

      allocate( UGC(DIMENSION_3) )
      allocate( VGC(DIMENSION_3) )
      allocate( WGC(DIMENSION_3) )

! Initialize fluid cell values.
      F_GDS = ZERO
      DRAG_BM = ZERO

! Loop bounds for interpolation.
      LP_BND = merge(27,9,DO_K)

! Calculate the cell center gas velocities.
      IF(DES_EXPLICITLY_COUPLED) THEN
         CALL CALC_CELL_CENTER_GAS_VEL(U_GO, V_GO, W_GO, UGC, VGC, WGC)
      ELSE
         CALL CALC_CELL_CENTER_GAS_VEL(U_G, V_G, W_G, UGC, VGC, WGC)
      ENDIF

! Calculate the gas phase forces acting on each particle.

!$omp parallel default(none) private(np,lepg,velfp,ijk,weight,ldrag_bm,lforce) &
!$omp          shared(max_pip,des_interp_on,lp_bnd,filter_cell,filter_weight,  &
!$omp          ep_g,pijk,des_vel_new,f_gp,vol,drag_bm,f_gds,ugc,vgc,wgc)
!$omp do
      DO NP=1,MAX_PIP
         IF(IS_NONEXISTENT(NP)) CYCLE

! The drag force is not calculated on entering or exiting particles
! as their velocities are fixed and may exist in 'non fluid' cells.
        IF(IS_ENTERING(NP) .OR. IS_EXITING(NP) .OR. IS_ENTERING_GHOST(NP) .OR. IS_EXITING_GHOST(NP)) CYCLE

         lEPG = ZERO
         VELFP = ZERO

! Calculate the gas volume fraction, velocity, and at the
! particle's position.
         IF(DES_INTERP_ON) THEN
            DO LC=1,LP_BND
               IJK = FILTER_CELL(LC,NP)
               WEIGHT = FILTER_WEIGHT(LC,NP)
! Gas phase volume fraction.
               lEPG = lEPG + EP_G(IJK)*WEIGHT
! Gas phase velocity.
               VELFP(1) = VELFP(1) + UGC(IJK)*WEIGHT
               VELFP(2) = VELFP(2) + VGC(IJK)*WEIGHT
               VELFP(3) = VELFP(3) + WGC(IJK)*WEIGHT
            ENDDO
         ELSE
            IJK = PIJK(NP,4)
            lEPG = EP_G(IJK)
            VELFP(1) = UGC(IJK)
            VELFP(2) = VGC(IJK)
            VELFP(3) = WGC(IJK)
         ENDIF

! This avoids FP exceptions for some ghost particles.
         IF(lEPg == ZERO) lEPG = EP_g(PIJK(NP,4))

! Calculate drag coefficient
         CALL DES_DRAG_GP(NP, DES_VEL_NEW(NP,:), VELFP, lEPg)

         lFORCE = F_GP(NP)

         lDRAG_BM = lFORCE*DES_VEL_NEW(NP,:)

         IF(DES_INTERP_ON) THEN
            DO LC=1,LP_BND
               IJK = FILTER_CELL(LC,NP)
               WEIGHT = FILTER_WEIGHT(LC,NP)/VOL(IJK)

               !$omp atomic
               DRAG_BM(IJK,1) = DRAG_BM(IJK,1) + lDRAG_BM(1)*WEIGHT
               !$omp atomic
               DRAG_BM(IJK,2) = DRAG_BM(IJK,2) + lDRAG_BM(2)*WEIGHT
               !$omp atomic
               DRAG_BM(IJK,3) = DRAG_BM(IJK,3) + lDRAG_BM(3)*WEIGHT
               !$omp atomic
               F_GDS(IJK) = F_GDS(IJK) + lFORCE*WEIGHT
            ENDDO
         ELSE
            IJK = PIJK(NP,4)
            WEIGHT = ONE/VOL(IJK)

            !$omp atomic
            DRAG_BM(IJK,1) = DRAG_BM(IJK,1) + lDRAG_BM(1)*WEIGHT
            !$omp atomic
            DRAG_BM(IJK,2) = DRAG_BM(IJK,2) + lDRAG_BM(2)*WEIGHT
            !$omp atomic
            DRAG_BM(IJK,3) = DRAG_BM(IJK,3) + lDRAG_BM(3)*WEIGHT

            !$omp atomic
            F_GDS(IJK) = F_GDS(IJK) + lFORCE*WEIGHT
         ENDIF

      ENDDO
!$omp end parallel

! Update the drag force and sources in ghost layers.
      CALL SEND_RECV(F_GDS, 2)
      CALL SEND_RECV(DRAG_BM, 2)

      deallocate( UGC, VGC, WGC )

      RETURN
      END SUBROUTINE DRAG_GS_GAS1

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_CELL_CENTER_GAS_VEL                                !
!  Author: J.Musser                                   Date: 07-NOV-14  !
!                                                                      !
!  Purpose: Calculate the scalar cell center gas velocity. This code   !
!  is common to the DEM and GAS calls for non-interpolated drag        !
!  routines.                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_CELL_CENTER_GAS_VEL(lUg, lVg, lWg, Uc, Vc, Wc)

! Global Variables:
!---------------------------------------------------------------------//
! Functions to average momentum to scalar cell center.
      use fun_avg, only: AVG_X_E, AVG_Y_N, AVG_Z_T
! Flags and correction factors for cut momentum cells.
      use cutcell, only: CUT_U_TREATMENT_AT, THETA_UE, THETA_UE_BAR
      use cutcell, only: CUT_V_TREATMENT_AT, THETA_VN, THETA_VN_BAR
      use cutcell, only: CUT_W_TREATMENT_AT, THETA_WT, THETA_WT_BAR
! Functions to lookup adjacent cells by index.
      use functions, only: IM_OF, JM_OF, KM_OF
      use indices, only: I_OF
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
! Flag for 3D simulatoins.
      use geometry, only: DO_K
! Function to deterine if a cell contains fluid.
      use functions, only: FLUID_AT
      use functions, only: funijk

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision parameters.
      use param, only: DIMENSION_3
      use param1, only: ZERO

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: lUg(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: lVg(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: lWg(DIMENSION_3)

      DOUBLE PRECISION, INTENT(OUT) :: Uc(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: Vc(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: Wc(DIMENSION_3)

! Local variables:
!---------------------------------------------------------------------//
! Indices of adjacent cells
      INTEGER :: i,j,k, IJK, IMJK, IJMK, IJKM

! Calculate the cell center gas velocity components.
        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
         IF(FLUID_AT(IJK)) THEN
            IMJK = IM_OF(IJK)
            IF(CUT_U_TREATMENT_AT(IMJK)) THEN
               Uc(IJK) = (THETA_UE_BAR(IMJK)*lUG(IMJK) +              &
                  THETA_UE(IMJK)*lUg(IJK))
            ELSE
               Uc(IJK) = AVG_X_E(lUG(IMJK),lUG(IJK),I_OF(IJK))
            ENDIF

            IJMK = JM_OF(IJK)
            IF(CUT_V_TREATMENT_AT(IJMK)) THEN
               Vc(IJK) = (THETA_VN_BAR(IJMK)*lVG(IJMK) +              &
                  THETA_VN(IJMK)*lVg(IJK))
            ELSE
               Vc(IJK) = AVG_Y_N(lVg(IJMK),lVg(IJK),0)
            ENDIF

            IF(DO_K) THEN
               IJKM = KM_OF(IJK)
               IF(CUT_W_TREATMENT_AT(IJKM)) THEN
                  Wc(IJK) = (THETA_WT_BAR(IJKM)*lWg(IJKM) +           &
                     THETA_WT(IJKM)* lWg(IJK))
               ELSE
                  Wc(IJK) = AVG_Z_T(lWg(IJKM),lWg(IJK),0)
               ENDIF
            ELSE
               Wc(IJK) = ZERO
            ENDIF
         ELSE
            Uc(IJK) = ZERO
            Vc(IJK) = ZERO
            Wc(IJK) = ZERO
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CALC_CELL_CENTER_GAS_VEL
