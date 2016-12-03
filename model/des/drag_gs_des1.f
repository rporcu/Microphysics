module drag_gs_des1_module

! Global Variables:
!---------------------------------------------------------------------//
! Functions to average momentum to scalar cell center.
      use fun_avg, only: AVG
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
! Function to deterine if a cell contains fluid.
      use functions, only: fluid_at
      use functions, only: funijk, iminus, jminus, kminus

! Size of particle array on this process.
      use discretelement, only: MAX_PIP
! IJK of fluid cell containing particles center
      use discretelement, only: PIJK
! Particle velocity
      use discretelement, only: DES_VEL_NEW
! Volume of scalar cell.
      use geometry, only: VOL

      use des_drag_gp_module, only: des_drag_gp

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO, ONE

! Total forces acting on particle
      use discretelement, only: FC
! Particle volume.
      use discretelement, only: PVOL
      use discretelement

  contains

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
      SUBROUTINE DRAG_GS_DES1(ep_g, u_g, v_g, w_g, ro_g, mu_g, gradPg)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: gradPg&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)

! Local variables
!---------------------------------------------------------------------//
! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK
! Interpolated gas phase quanties.
      DOUBLE PRECISION :: lEPg, VELFP(3)
! Drag force acting on each particle.
      DOUBLE PRECISION :: D_FORCE(3), f_gp
! Loop bound for filter
      integer :: i,j,k
!......................................................................!

! Calculate the gas phase forces acting on each particle.
      DO NP=1,MAX_PIP
         IF(.NOT.NORMAL_PARTICLE==PARTICLE_STATE(NP)) CYCLE
! Avoid drag calculations in cells without fluid (cut-cell)
         i = pijk(np,1)
         j = pijk(np,2)
         k = pijk(np,3)
         if (.NOT.fluid_at(i,j,k)) CYCLE

         lEPG = ZERO
         VELFP = ZERO

! Calculate the gas volume fraction, velocity, and pressure force at
! the particle's position.
         IJK = PIJK(NP,4)
         lEPG = EP_G(I,J,K)
         VELFP(1) = 0.5d0*(u_g(iminus(i,j,k),j,k) + u_g(I,J,K))
         VELFP(2) = 0.5d0*(v_g(i,jminus(i,j,k),k) + v_g(I,J,K))
         VELFP(3) = 0.5d0*(w_g(i,j,kminus(i,j,k)) + w_g(I,J,K))

! For explicit coupling, use the drag coefficient calculated for the
! gas phase drag calculations.

! Calculate the drag coefficient.
         CALL DES_DRAG_GP(NP, DES_VEL_NEW(NP,:), VELFP, lEPg, ro_g, mu_g, f_gp)

! Calculate the gas-solids drag force on the particle
         D_FORCE = F_GP*(VELFP - DES_VEL_NEW(NP,:))

! Update the contact forces (FC) on the particle to include gas
! pressure and gas-solids drag
         FC(NP,:) = FC(NP,:) + D_FORCE(:) + PVOL(NP)*gradPg(I,J,K,:)

      ENDDO

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
      SUBROUTINE DRAG_GS_GAS1(ep_g, u_g, v_g, w_g, ro_g, mu_g, &
         f_gds, drag_am, drag_bm)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT  ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT  ) :: drag_am&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT  ) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)

! Local variables
!---------------------------------------------------------------------//
! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK, I, J, K
! Interpolation weight
      DOUBLE PRECISION :: WEIGHT
! Interpolated gas phase quanties.
      DOUBLE PRECISION :: lEPg, VELFP(3)
! Drag force (intermediate calculation)
      DOUBLE PRECISION :: f_gp
! Drag sources for fluid (intermediate calculation)
      DOUBLE PRECISION :: lDRAG_BM(3)
!......................................................................!

! Initialize fluid cell values.
      F_GDS = ZERO
      DRAG_BM = ZERO

! Calculate the gas phase forces acting on each particle.

      DO NP=1,MAX_PIP
         IF(NONEXISTENT==PARTICLE_STATE(NP)) CYCLE

! The drag force is not calculated on entering or exiting particles
! as their velocities are fixed and may exist in 'non fluid' cells.
         IF(ENTERING_PARTICLE==PARTICLE_STATE(NP) .OR. &
            EXITING_PARTICLE==PARTICLE_STATE(NP) .OR. &
            ENTERING_GHOST==PARTICLE_STATE(NP) .OR. &
            EXITING_GHOST==PARTICLE_STATE(NP)) CYCLE

         lEPG = ZERO
         VELFP = ZERO
! This avoids FP exceptions for some ghost particles.
         I = PIJK(NP,1)
         J = PIJK(NP,2)
         K = PIJK(NP,3)

! Calculate the gas volume fraction, velocity, and at the
! particle's position.
         IJK = PIJK(NP,4)
         lEPG = EP_G(I,J,K)
         VELFP(1) = 0.5d0*(u_g(iminus(i,j,k),j,k) + u_g(I,J,K))
         VELFP(2) = 0.5d0*(v_g(i,jminus(i,j,k),k) + v_g(I,J,K))
         VELFP(3) = 0.5d0*(w_g(i,j,kminus(i,j,k)) + w_g(I,J,K))

         IF(lEPg == ZERO) lEPG = EP_g(I,J,K)

! Calculate drag coefficient
         CALL DES_DRAG_GP(NP, DES_VEL_NEW(NP,:), VELFP, lEPg, ro_g, mu_g, f_gp)

         lDRAG_BM = f_gp*DES_VEL_NEW(NP,:)

         IJK = PIJK(NP,4)
         WEIGHT = ONE/VOL

         DRAG_BM(I,J,K,1) = DRAG_BM(I,J,K,1) + lDRAG_BM(1)*WEIGHT
         DRAG_BM(I,J,K,2) = DRAG_BM(I,J,K,2) + lDRAG_BM(2)*WEIGHT
         DRAG_BM(I,J,K,3) = DRAG_BM(I,J,K,3) + lDRAG_BM(3)*WEIGHT

         F_GDS(i,j,k) = F_GDS(i,j,k) + f_gp*WEIGHT

      ENDDO

! Update the drag force and sources in ghost layers.
      ! CALL SEND_RECV(F_GDS, 2)
      ! CALL SEND_RECV(DRAG_BM, 2)

      RETURN
      END SUBROUTINE DRAG_GS_GAS1
end module drag_gs_des1_module
