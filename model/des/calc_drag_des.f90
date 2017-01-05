module calc_drag_des_module

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use comp_mean_fields_module, only: comp_mean_fields
      use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
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
     SUBROUTINE CALC_DRAG_DES(ep_g, u_g, v_g, w_g, ro_g, mu_g, gradPg,&
         particle_state, fc, drag_fc, pvol, des_pos_new, &
        des_vel_new, des_radius, particle_phase, flag)

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
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)

      integer         , intent(in   ) :: particle_phase(:)
      integer         , intent(in   ) :: particle_state(:)

      real(c_real), intent(in   ) :: pvol(:)
      real(c_real), intent(in   ) :: des_radius(:)
      real(c_real), intent(in   ) :: drag_fc(:,:)
      real(c_real), intent(in   ) :: des_vel_new(:,:)
      real(c_real), intent(in   ) :: des_pos_new(:,:)
      real(c_real), intent(inout) :: fc(:,:)

      INTEGER :: II

! Apply the drag force calculated by the gas phase.
      if(des_explicitly_coupled) then

         if(des_continuum_coupled) then
            do ii = 1, size(pvol)
               if(normal_particle==particle_state(ii)) &
                  fc(ii,:) = fc(ii,:) + drag_fc(ii,:)
            enddo
         endif

      else

! Calculate gas-solids drag force on particle
         if(des_continuum_coupled) then
            call drag_gs_des(ep_g, u_g, v_g, w_g, ro_g, mu_g, &
               gradpg, flag, particle_state, pvol, des_pos_new, &
               des_vel_new, fc, des_radius, particle_phase)
         endif

      endif

      return
      END SUBROUTINE CALC_DRAG_DES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_DRAG_GS                                             !
!                                                                      !
!  Purpose: This subroutine is only called from the CONTINUUM side. It !
!  calls the correct routine for calculating the gas drag force.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DRAG_DES_2FLUID(ep_g, u_g, v_g, w_g, ro_g, mu_g, &
         f_gds, drag_bm,  particle_state, &
         particle_phase, pvol, des_pos_new, des_vel_new, des_radius)

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
      real(c_real), intent(in   ) :: des_vel_new(:,:)
      real(c_real), intent(in   ) :: des_pos_new(:,:)
      integer         , intent(in   ) :: particle_state(:)
      integer         , intent(in   ) :: particle_phase(:)

! Calculate gas-solids drag force.
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL DRAG_GS_GAS(ep_g, u_g, v_g, w_g, ro_g, mu_g, &
            f_gds, drag_bm, particle_phase, particle_state, pvol, &
            des_pos_new, des_vel_new, des_radius)
      ENDIF

      RETURN
      END SUBROUTINE CALC_DRAG_DES_2FLUID


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
      SUBROUTINE CALC_DRAG_DES_EXPLICIT(flag, ep_g, u_g, v_g, w_g, ro_g, &
         mu_g, f_gds, drag_bm,  particle_phase,  &
         particle_state, pvol, des_pos_new, des_vel_new, &
         des_radius)

      IMPLICIT NONE

      real(c_real), intent(inout) :: ep_g&
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
      integer         , intent(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)

      integer, intent(in   ) :: particle_state(:)

      integer, intent(in   ) :: particle_phase(:)

      real(c_real), intent(in   ) :: pvol(:)
      real(c_real), intent(in   ) :: des_radius(:)

      real(c_real), intent(in   ) :: des_pos_new(:,:)
      real(c_real), intent(in   ) :: des_vel_new(:,:)

! Calculate mean fields (EPg).
      CALL COMP_MEAN_FIELDS(ep_g, particle_state, des_pos_new, pvol, flag, size(des_radius))

! Calculate gas-solids drag force on particle
      IF(DES_CONTINUUM_COUPLED) &
         CALL DRAG_GS_GAS(ep_g, u_g, v_g, w_g, ro_g, mu_g, &
            f_gds, drag_bm, particle_phase, particle_state, pvol, &
            des_pos_new, des_vel_new, des_radius)

      END SUBROUTINE CALC_DRAG_DES_EXPLICIT

end module calc_drag_des_module
