module drag_gs_des1_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

! Global Variables:
!---------------------------------------------------------------------//
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
! Function to deterine if a cell contains fluid.
      use functions, only: iminus, jminus, kminus

! Volume of scalar cell.
      use geometry, only: VOL

      use des_drag_gp_module, only: des_drag_gp

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO, ONE

      use discretelement, only: entering_particle, entering_ghost, exiting_particle
      use discretelement, only: nonexistent, exiting_ghost, normal_particle

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
     SUBROUTINE DRAG_GS_DES(ep_g, u_g, v_g, w_g, ro_g, mu_g, gradPg, &
        flag, particle_state, pvol, des_pos_new, des_vel_new, fc, &
        des_radius, particle_phase, dx, dy, dz)


      IMPLICIT NONE

      real(c_real), intent(in   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: gradpg&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
      integer         , intent(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)

      real(c_real), intent(in   ) :: pvol(:)
      real(c_real), intent(in   ) :: des_radius(:)
      real(c_real), intent(in   ) :: des_pos_new(:,:)
      real(c_real), intent(in   ) :: des_vel_new(:,:)
      real(c_real), intent(inout) :: fc(:,:)

      integer     , intent(in   ) :: particle_state(:)
      integer     , intent(in   ) :: particle_phase(:)

      real(c_real), intent(in   ) :: dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP
! Interpolated gas phase quanties.
      real(c_real) :: lEPg, VELFP(3)
! Drag force acting on each particle.
      real(c_real) :: D_FORCE(3), f_gp
! Loop bound for filter
      integer :: i,j,k
! One over cell volume
      real(c_real) :: odx, ody, odz
!......................................................................!

      odx = 1.0d0/dx
      ody = 1.0d0/dy
      odz = 1.0d0/dz

! Calculate the gas phase forces acting on each particle.
      DO NP=1,size(pvol)

         IF(.NOT.NORMAL_PARTICLE==PARTICLE_STATE(NP)) CYCLE

         i = floor(des_pos_new(np,1)*odx) + 1
         j = floor(des_pos_new(np,2)*ody) + 1
         k = floor(des_pos_new(np,3)*odz) + 1

! Avoid drag calculations in cells without fluid (cut-cell)
         if (flag(i,j,k,1)/=1) CYCLE

! Calculate the gas volume fraction, velocity, and pressure force at
! the particle's position.
         lEPG = EP_G(I,J,K)
         VELFP(1) = 0.5d0*(u_g(iminus(i,j,k),j,k) + u_g(I,J,K))
         VELFP(2) = 0.5d0*(v_g(i,jminus(i,j,k),k) + v_g(I,J,K))
         VELFP(3) = 0.5d0*(w_g(i,j,kminus(i,j,k)) + w_g(I,J,K))

! For explicit coupling, use the drag coefficient calculated for the
! gas phase drag calculations.

! Calculate the drag coefficient.
         call des_drag_gp(np, des_vel_new(np,:), velfp, lepg, ro_g, mu_g, &
            f_gp, i,j,k, des_radius,  pvol, particle_phase)

! Calculate the gas-solids drag force on the particle
         D_FORCE = F_GP*(VELFP - DES_VEL_NEW(NP,:))

! Update the contact forces (FC) on the particle to include gas
! pressure and gas-solids drag
         FC(NP,:) = FC(NP,:) + D_FORCE(:) + PVOL(NP)*gradPg(I,J,K,:)

      ENDDO

      RETURN
      END SUBROUTINE DRAG_GS_DES

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
      SUBROUTINE DRAG_GS_GAS(ep_g, u_g, v_g, w_g, ro_g, mu_g, &
         f_gds, drag_bm, particle_phase, particle_state, &
         pvol, des_pos_new, des_vel_new, des_radius, dx, dy, dz)

      IMPLICIT NONE

      real(c_real), intent(in   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(out  ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(out  ) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)

      real(c_real), intent(in   ) :: pvol(:)
      real(c_real), intent(in   ) :: des_radius(:)
      real(c_real), intent(in   ) :: des_pos_new(:,:)
      real(c_real), intent(in   ) :: des_vel_new(:,:)

      integer     , intent(in   ) :: particle_state(:)
      integer     , intent(in   ) :: particle_phase(:)

      real(c_real), intent(in   ) :: dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, I, J, K
! Interpolation weight
      real(c_real) :: WEIGHT
! Interpolated gas phase quanties.
      real(c_real) :: lEPg, VELFP(3)
! Drag force (intermediate calculation)
      real(c_real) :: f_gp
! Drag sources for fluid (intermediate calculation)
      real(c_real) :: lDRAG_BM(3)
! One over cell volume
      real(c_real) :: odx, ody, odz
!......................................................................!

! Initialize fluid cell values.
      F_GDS = ZERO
      DRAG_BM = ZERO

      odx = 1.0d0/dx
      ody = 1.0d0/dy
      odz = 1.0d0/dz

! Calculate the gas phase forces acting on each particle.

      DO NP=1,size(pvol)

         IF(NONEXISTENT==PARTICLE_STATE(NP)) CYCLE

! The drag force is not calculated on entering or exiting particles
! as their velocities are fixed and may exist in 'non fluid' cells.
         if(entering_particle==particle_state(np) .or. &
            exiting_particle==particle_state(np) .or. &
            entering_ghost==particle_state(np) .or. &
            exiting_ghost==particle_state(np)) cycle

! This avoids FP exceptions for some ghost particles.
         i = floor(des_pos_new(np,1)*odx) + 1
         j = floor(des_pos_new(np,2)*ody) + 1
         k = floor(des_pos_new(np,3)*odz) + 1

! Calculate the gas volume fraction, velocity, and at the
! particle's position.
         lEPG = EP_G(I,J,K)
         VELFP(1) = 0.5d0*(u_g(iminus(i,j,k),j,k) + u_g(I,J,K))
         VELFP(2) = 0.5d0*(v_g(i,jminus(i,j,k),k) + v_g(I,J,K))
         VELFP(3) = 0.5d0*(w_g(i,j,kminus(i,j,k)) + w_g(I,J,K))

         IF(lEPg < EPSILON(lEPg)) lEPG = EP_g(I,J,K)

! Calculate drag coefficient
         CALL DES_DRAG_GP(NP, DES_VEL_NEW(NP,:), VELFP, lEPg, ro_g, mu_g,&
            f_gp, i,j,k, des_radius,  pvol, particle_phase)

         lDRAG_BM = f_gp*DES_VEL_NEW(NP,:)

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
      END SUBROUTINE DRAG_GS_GAS
end module drag_gs_des1_module
