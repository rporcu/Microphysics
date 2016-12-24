MODULE MASS_INFLOW_DEM_MODULE

      use bc, only: bc_i_w, bc_j_s, bc_k_b
      use bc, only: bc_plane
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: bc_u_s, bc_v_s, bc_w_s
      use constant, only: PI
      use des_bc, only: dem_bc_poly_layout, dem_bcmi_map, dem_bcmi_ijk, dem_bcmi_ijkstart, dem_mi_time, dem_bcmi_ijkend
      use des_bc, only: dem_mi, dem_bcmi, numfrac_limit, pi_count, pi_factor
      use discretelement, only: des_explicitly_coupled, xe, yn, zt
      use discretelement, only: max_pip, dtsolid, imax_global_id
      use discretelement, only:  s_time, pip, normal_particle
      use discretelement, only: entering_particle, entering_ghost, normal_ghost, exiting_particle, exiting_ghost, nonexistent
      use param1, only: half, zero
      use constant, only: d_p0, ro_s0

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_MASS_INLET                                          !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose:  This routine fills in the necessary information for new   !
!  particles entering the system.                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MASS_INFLOW_DEM( particle_phase,  iglobal_id, particle_state, &
         des_radius, omoi, pmass, pvol, ro_sol, &
         des_vel_new, des_pos_new, ppos, omega_new, drag_fc)

      implicit none

      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: des_radius, omoi, pmass, pvol, ro_sol
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: des_vel_new, des_pos_new, ppos, omega_new, drag_fc
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(OUT) ::  iglobal_id
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_phase

      INTEGER :: IP, LS, M, NP, IJK, LC
      INTEGER :: BCV, BCV_I
      LOGICAL :: CHECK_FOR_ERRORS, OWNS

! I/J/K index of fluid cell containing the new particle.
      INTEGER :: IJKP(3)

      DOUBLE PRECISION :: DIST, POS(3)
! Random numbers -shared by all ranks-
      DOUBLE PRECISION :: RAND(3)

! HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK
      return
! HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK
      RETURN
      END SUBROUTINE MASS_INFLOW_DEM

END MODULE MASS_INFLOW_DEM_MODULE
