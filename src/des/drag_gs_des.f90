module drag_gs_des1_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

! Global Variables:
!---------------------------------------------------------------------//
! Function to deterine if a cell contains fluid.

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
     subroutine drag_gs_des(slo, shi,ulo, uhi, vlo, vhi, wlo, whi,&
        max_pip, ep_g, u_g, v_g, w_g, ro_g, mu_g, gradpg, &
        particle_state, pvol, des_pos_new, des_vel_new,  &
        fc, des_radius, particle_phase, dx, dy, dz)

        IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)
      integer(c_int), intent(in   ) :: vlo(3),vhi(3)
      integer(c_int), intent(in   ) :: wlo(3),whi(3)
      integer(c_int), intent(in   ) :: max_pip

      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: gradpg&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      real(c_real), intent(in   ) :: pvol(max_pip)
      real(c_real), intent(in   ) :: des_radius(max_pip)
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)
      real(c_real), intent(inout) :: fc(:,:)

      integer     , intent(in   ) :: particle_state(max_pip)
      integer     , intent(in   ) :: particle_phase(max_pip)

      real(c_real), intent(in   ) :: dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Loop counters: Particle, fluid cell, neighbor cells
      integer :: np
! Interpolated gas phase quanties.
      real(c_real) :: lEPg, VELFP(3)
! Drag force acting on each particle.
      real(c_real) :: D_FORCE(3), f_gp
! Loop bound for filter
      integer :: i,j,k
! One over cell volume
      real(c_real) :: odx, ody, odz, vol
      real(c_real) :: vel(3)
!......................................................................!

      odx = 1.0d0/dx
      ody = 1.0d0/dy
      odz = 1.0d0/dz

      vol = dx*dy*dz

! Calculate the gas phase forces acting on each particle.
      do np = 1,max_pip

         IF(.NOT.NORMAL_PARTICLE==PARTICLE_STATE(np)) CYCLE

         i = floor(des_pos_new(np,1)*odx) - 1
         j = floor(des_pos_new(np,2)*ody) - 1
         k = floor(des_pos_new(np,3)*odz) - 1

         ! Calculate the gas volume fraction, velocity, and pressure force at
         ! the particle's position.
         lEPG = EP_G(I,J,K)
         VELFP(1) = 0.5d0*(u_g(i-1,j,k) + u_g(I,J,K))
         VELFP(2) = 0.5d0*(v_g(i,j-1,k) + v_g(I,J,K))
         VELFP(3) = 0.5d0*(w_g(i,j,k-1) + w_g(I,J,K))

         ! For explicit coupling, use the drag coefficient calculated for the
         ! gas phase drag calculations.

         ! Calculate the drag coefficient.
         vel(:)=des_vel_new(np,:)
         call des_drag_gp(slo, shi, np, vel, velfp, lepg,&
                           ro_g, mu_g, f_gp, i, j, k, &
                           des_radius(np),  pvol(np), particle_phase(np))

         ! Calculate the gas-solids drag force on the particle
         D_FORCE = F_GP*(VELFP - DES_VEL_NEW(np,:))

         ! Update the contact forces (FC) on the particle to include gas
         ! pressure and gas-solids drag
         FC(np,:) = FC(np,:) + D_FORCE(:) + PVOL(np)*gradPg(I,J,K,:)

      end do

      end subroutine drag_gs_des

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: drag_gs_gas                                            !
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
      subroutine drag_gs_gas(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
         max_pip, ep_g, u_g, v_g, w_g, ro_g, mu_g, &
         f_gds, drag_bm, particle_phase, particle_state, &
         pvol, des_pos_new, des_vel_new, des_radius, dx, dy, dz)

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3)
      integer(c_int), intent(in   ) :: vlo(3),vhi(3)
      integer(c_int), intent(in   ) :: wlo(3),whi(3)
      integer(c_int), intent(in   ) :: max_pip

      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(out  ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(out  ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      real(c_real), intent(in   ) :: pvol(max_pip)
      real(c_real), intent(in   ) :: des_radius(max_pip)
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)

      integer     , intent(in   ) :: particle_state(max_pip)
      integer     , intent(in   ) :: particle_phase(max_pip)

      real(c_real), intent(in   ) :: dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! Loop counters: Particle, fluid cell, neighbor cells
      integer :: np, I, J, K
! Interpolation weight
      real(c_real) :: WEIGHT
! Interpolated gas phase quanties.
      real(c_real) :: lEPg, VELFP(3)
! Drag force (intermediate calculation)
      real(c_real) :: f_gp
! Drag sources for fluid (intermediate calculation)
      real(c_real) :: ldrag_bm(3)
! One over cell volume
      real(c_real) :: odx, ody, odz, vol
!......................................................................!
      real(c_real) :: particle_vel(3)
!......................................................................!

! Initialize fluid cell values.
      f_gds = ZERO
      drag_bm = ZERO

      odx = 1.0d0/dx
      ody = 1.0d0/dy
      odz = 1.0d0/dz

      vol = dx*dy*dz

! Calculate the gas phase forces acting on each particle.

      do np=1,max_pip

         if(nonexistent==particle_state(np)) cycle

         ! The drag force is not calculated on entering or exiting particles
         ! as their velocities are fixed and may exist in 'non fluid' cells.
         if(entering_particle==particle_state(np) .or. &
            exiting_particle==particle_state(np) .or. &
            entering_ghost==particle_state(np) .or. &
            exiting_ghost==particle_state(np)) cycle

         ! This avoids FP exceptions for some ghost particles.
         i = floor(des_pos_new(np,1)*odx) - 1
         j = floor(des_pos_new(np,2)*ody) - 1
         k = floor(des_pos_new(np,3)*odz) - 1

         ! Calculate the gas volume fraction, velocity, and at the
         ! particle's position.
         lepg = ep_g(i,j,k)
         velfp(1) = 0.5d0*(u_g(i-1,j,k) + u_g(i,j,k))
         velfp(2) = 0.5d0*(v_g(i,j-1,k) + v_g(i,j,k))
         velfp(3) = 0.5d0*(w_g(i,j,k-1) + w_g(i,j,k))

         if(lepg < epsilon(lepg)) lepg = ep_g(i,j,k)

         ! Calculate drag coefficient
         particle_vel(:) = des_vel_new(np,:)
         call des_drag_gp(slo, shi, np, particle_vel, velfp, lepg, &
                          ro_g, mu_g, f_gp, i, j, k, &
                          des_radius(np),  pvol(np), particle_phase(np))

         ldrag_bm = f_gp*des_vel_new(np,:)

         weight = one/vol

         drag_bm(i,j,k,1) = drag_bm(i,j,k,1) + ldrag_bm(1)*weight
         drag_bm(i,j,k,2) = drag_bm(i,j,k,2) + ldrag_bm(2)*weight
         drag_bm(i,j,k,3) = drag_bm(i,j,k,3) + ldrag_bm(3)*weight

         f_gds(i,j,k) = f_gds(i,j,k) + f_gp*weight

      enddo

      end subroutine drag_gs_gas

end module drag_gs_des1_module
