module calc_drag_des_module

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use discretelement, only: des_continuum_coupled
      use discretelement, only: des_explicitly_coupled
      use discretelement, only: normal_particle

      use drag_gs_des1_module, only: drag_gs_des, drag_gs_gas

  contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_DRAG_DES                                           !
!                                                                      !
!  Purpose: This subroutine is called from DES routines. It calls      !
!  functions that calcultate the drag force acting on particles. No    !
!  field variables are updated.                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     subroutine calc_drag_des(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
        max_pip, ep_g, u_g, v_g, w_g, ro_g, mu_g, gradPg, particle_state,&
        fc, drag_fc, pvol, des_pos_new, des_vel_new, des_radius,&
        particle_phase, dx, dy, dz)

        IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: wlo(3), whi(3)
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

      integer     , intent(in   ) :: particle_phase(max_pip)
      integer     , intent(in   ) :: particle_state(max_pip)

      real(c_real), intent(in   ) :: pvol(max_pip)
      real(c_real), intent(in   ) :: des_radius(max_pip)
      real(c_real), intent(in   ) :: drag_fc(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(inout) :: fc(max_pip,3)

      real(c_real), intent(in   ) :: dx, dy, dz

      integer :: II

! Apply the drag force calculated by the gas phase.
      if(des_explicitly_coupled) then

         if(des_continuum_coupled) then
            do ii = 1, max_pip
               if(normal_particle==particle_state(ii)) &
                  fc(ii,:) = fc(ii,:) + drag_fc(ii,:)
            enddo
         endif

      else

! Calculate gas-solids drag force on particle
         if(des_continuum_coupled) then
            call drag_gs_des(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
               max_pip, ep_g, u_g, v_g, w_g, ro_g, mu_g, gradpg, &
               particle_state, pvol, des_pos_new, des_vel_new, fc, des_radius, &
               particle_phase, dx, dy, dz)
         endif

      endif

      end subroutine calc_drag_des


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_DRAG_GS                                             !
!                                                                      !
!  Purpose: This subroutine is only called from the CONTINUUM side. It !
!  calls the correct routine for calculating the gas drag force.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_drag_des_2fluid(slo, shi, ulo, uhi, vlo, vhi, wlo, whi,&
         max_pip, ep_g, u_g, v_g, w_g, &
         ro_g, mu_g, f_gds, drag_bm,  particle_state, particle_phase,&
         pvol, des_pos_new, des_vel_new, des_radius, dx, dy, dz)

         IMPLICIT NONE

      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: ulo(3), uhi(3)
      integer, intent(in   ) :: vlo(3), vhi(3)
      integer, intent(in   ) :: wlo(3), whi(3)
      integer, intent(in   ) :: max_pip

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
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      integer     , intent(in   ) :: particle_state(max_pip)
      integer     , intent(in   ) :: particle_phase(max_pip)

      real(c_real), intent(in   ) :: dx, dy, dz

! Calculate gas-solids drag force.
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL DRAG_GS_GAS(slo, shi,ulo, uhi, vlo, vhi, wlo, whi, max_pip,&
            ep_g, u_g, v_g, w_g, ro_g, mu_g, f_gds, drag_bm, particle_phase,&
            particle_state, pvol, des_pos_new, des_vel_new, des_radius, &
            dx, dy, dz)
      ENDIF

      end subroutine calc_drag_des_2fluid


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_DRAG_DES_EXPLICIT                                  !
!                                                                      !
!  Purpose: This subroutine is only called from the CONTINUUM side.    !
!  Moreover, it is only called once per time step so that the total    !
!  drag force is calculated explicitly for the fluid phase and DEM     !
!  particles. This is done to ensure momentum conservation and reduce  !
!  computational overhead.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_drag_des_explicit(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, max_pip, &
         ep_g, u_g, v_g, w_g, ro_g, mu_g, f_gds, drag_bm,  &
         particle_phase, particle_state, pvol, des_pos_new, des_vel_new, &
         des_radius, dx, dy, dz)

      use comp_mean_fields_module, only: comp_mean_fields

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: wlo(3), whi(3)
      integer(c_int), intent(in   ) :: max_pip

      real(c_real), intent(inout) :: ep_g&
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

      integer, intent(in   ) :: particle_state(max_pip)
      integer, intent(in   ) :: particle_phase(max_pip)

      real(c_real), intent(in   ) :: pvol(max_pip)
      real(c_real), intent(in   ) :: des_radius(max_pip)
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)

      real(c_real), intent(in   ) :: dx, dy, dz

      ! Calculate mean fields (EPg).
      call comp_mean_fields(slo, shi, max_pip, &
         ep_g, particle_state, des_pos_new, &
         pvol, dx, dy, dz)

! Calculate gas-solids drag force on particle
      IF(DES_CONTINUUM_COUPLED) &
         call drag_gs_gas(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, max_pip, &
            ep_g, u_g, v_g, w_g, ro_g, mu_g, &
            f_gds, drag_bm, particle_phase, particle_state, pvol, &
            des_pos_new, des_vel_new, des_radius, dx, dy, dz)

      end subroutine calc_drag_des_explicit

end module calc_drag_des_module
