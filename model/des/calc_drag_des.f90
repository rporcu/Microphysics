module calc_drag_des_module

      use comp_mean_fields_module, only: comp_mean_fields
      use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: DES_EXPLICITLY_COUPLED
      use discretelement, only: NORMAL_PARTICLE

      use drag_gs_des1_module, only: drag_gs_des1, drag_gs_gas1


      use particles_in_cell_module, only: particles_in_cell

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
     SUBROUTINE CALC_DRAG_DES(ep_g,u_g,v_g,w_g,ro_g,mu_g,gradPg,pijk,particle_state,&
        fc,drag_fc,pvol, &
        des_pos_new,des_vel_new,des_radius,particle_phase,flag,pinc)

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
      INTEGER         , INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)
      INTEGER         , INTENT(INOUT) :: pinc&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pvol,des_radius
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: drag_fc, des_vel_new, des_pos_new
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: fc
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(IN) :: particle_phase
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk

      INTEGER :: II

! Apply the drag force calculated by the gas phase.
      IF(DES_EXPLICITLY_COUPLED) THEN

         IF(DES_CONTINUUM_COUPLED) THEN
            DO II = 1, size(pvol)
               IF(NORMAL_PARTICLE==PARTICLE_STATE(II)) &
                  FC(II,:) = FC(II,:) + DRAG_FC(II,:)
            ENDDO
         ENDIF

      ELSE

! Calculate gas-solids drag force on particle
         IF(DES_CONTINUUM_COUPLED) THEN
            CALL DRAG_GS_DES1(ep_g, u_g, v_g, w_g, ro_g, mu_g, &
               gradPg, flag, pijk, particle_state, pvol, des_vel_new,&
               fc, des_radius, particle_phase)
         ENDIF

      ENDIF

      RETURN
      END SUBROUTINE CALC_DRAG_DES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_DRAG_GS                                             !
!                                                                      !
!  Purpose: This subroutine is only called from the CONTINUUM side. It !
!  calls the correct routine for calculating the gas drag force.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DRAG_DES_2FLUID(ep_g, u_g, v_g, w_g, ro_g, mu_g, flag, &
         f_gds, drag_am, drag_bm, pijk, particle_state, particle_phase, &
         pvol, des_pos_new, des_vel_new, des_radius,  pinc)

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

      INTEGER         , INTENT(INOUT) :: pinc&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)


      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pvol, des_radius
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(IN) :: particle_phase
      INTEGER, DIMENSION(:,:,:,:), INTENT(IN) :: FLAG
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk

! Calculate gas-solids drag force.
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL DRAG_GS_GAS1(ep_g, u_g, v_g, w_g, ro_g, mu_g, &
            f_gds, drag_bm, pijk, particle_phase, particle_state, pvol, &
            des_vel_new, des_radius)
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
      SUBROUTINE CALC_DRAG_DES_EXPLICIT(flag, vol_surr, ep_g, u_g, v_g, w_g, ro_g, &
         rop_g, mu_g, f_gds, drag_bm, pijk, particle_phase, iglobal_id, &
         particle_state, pmass, pvol, des_pos_new, des_vel_new, &
         des_radius, des_usr_var, pinc)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT  ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT  ) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
      INTEGER         , INTENT(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)
      DOUBLE PRECISION, INTENT(in   ) :: vol_surr&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      INTEGER,          INTENT(INOUT) :: pinc&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pvol,pmass,des_radius
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_pos_new, des_vel_new, des_usr_var
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(OUT) :: iglobal_id
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_phase
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk

! Bin particles to the fluid grid.
      CALL PARTICLES_IN_CELL(pijk, iglobal_id, particle_state, &
         des_pos_new, des_vel_new, des_radius, des_usr_var, pinc)
! Calculate mean fields (EPg).
      CALL COMP_MEAN_FIELDS(ep_g,ro_g,rop_g,pijk,particle_state,particle_phase,pmass,pvol, &
         des_pos_new,des_vel_new,des_radius,&
         des_usr_var,flag,vol_surr,iglobal_id,pinc)

! Calculate gas-solids drag force on particle
      IF(DES_CONTINUUM_COUPLED) CALL DRAG_GS_GAS1(ep_g, u_g, v_g, w_g, &
         ro_g, mu_g, f_gds, drag_bm, pijk, particle_phase, particle_state, pvol, &
         des_vel_new, des_radius)

      END SUBROUTINE CALC_DRAG_DES_EXPLICIT

end module calc_drag_des_module
