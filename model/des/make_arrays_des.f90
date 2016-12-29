MODULE MAKE_ARRAYS_DES_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MAKE_ARRAYS_DES(ep_g, flag, &
             particle_state, particle_phase, &
         des_radius,  ro_sol, pvol, pmass, omoi, &
         des_pos_new, des_vel_new, des_usr_var, omega_new, fc)

      USE comp_mean_fields_module, only: comp_mean_fields
      USE compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
      USE compar, only: numpes, mype
      USE constant, only: pi
      USE discretelement, only: do_nsearch, particles, max_pip, vtp_findex
      USE discretelement, only: entering_ghost, exiting_ghost, nonexistent, normal_ghost
      USE discretelement, only: print_des_data, s_time
      USE error_manager, only: err_msg, flush_err_msg, init_err_msg, finl_err_msg
      USE functions, only: ip1, jp1, kp1
      USE param1, only: zero

      USE read_par_input_module, only: read_par_input
      USE read_res0_des_module, only: read_res0_des
      USE run, only: run_type
      USE set_phase_index_module, only: set_phase_index

      USE write_des_data_module, only: write_des_data

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer, intent(inout) :: flag &
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)

      double precision, intent(  out) :: pvol(:)
      double precision, intent(  out) :: pmass(:)
      double precision, intent(  out) :: des_radius(:)
      double precision, intent(  out) :: ro_sol(:)
      double precision, intent(  out) :: omoi(:)

      double precision, intent(inout) :: fc(:,:)
      double precision, intent(  out) :: des_vel_new(:,:)
      double precision, intent(  out) :: des_pos_new(:,:)
      double precision, intent(  out) :: omega_new(:,:)
      double precision, intent(  out) :: des_usr_var(:,:)
      integer         , intent(  out) :: particle_state(:)
      integer         , intent(  out) :: particle_phase(:)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: L

      CALL INIT_ERR_MSG("MAKE_ARRAYS_DES")

! Set the initial particle data.
      IF(RUN_TYPE == 'NEW') THEN
         IF(PARTICLES /= 0) CALL READ_PAR_INPUT(particle_state, &
            des_radius, ro_sol, des_pos_new, des_vel_new)

! Initialize old values
         omega_new(:,:)   = zero

! Read the restart file.
      ELSEIF(RUN_TYPE == 'RESTART_1' .OR. RUN_TYPE == 'RESTART_2') THEN

         CALL READ_RES0_DES(   particle_state, &
            des_radius, ro_sol, des_usr_var, &
            des_pos_new, des_vel_new, omega_new)

      ELSE

         WRITE(ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1100 FORMAT('Error 1100: Unsupported RUN_TYPE for DES.')

      ENDIF

      IF(RUN_TYPE == 'RESTART_2') VTP_FINDEX=0

! setting additional particle properties now that the particles
! have been identified
      DO L = 1, MAX_PIP
! Skip 'empty' locations when populating the particle property arrays.
         if(nonexistent==particle_state(l) .or. &
            normal_ghost==particle_state(l) .or. &
            entering_ghost==particle_state(l) .or. &
            exiting_ghost==particle_state(l)) cycle
         pvol(l) = (4.0d0/3.0d0)*pi*des_radius(l)**3
         pmass(l) = pvol(l)*ro_sol(l)
         omoi(l) = 2.5d0/(pmass(l)*des_radius(l)**2) !one over moi
      ENDDO

      CALL SET_PHASE_INDEX(particle_phase,des_radius,ro_sol,particle_state)

! do_nsearch should be set before calling particle in cell
      DO_NSEARCH =.TRUE.

!      CALL NEIGHBOUR(  particle_state, des_radius, des_pos_new)

! Calculate mean fields using either interpolation or cell averaging.
      CALL COMP_MEAN_FIELDS(ep_g, particle_state, des_pos_new, pvol, flag)

      IF(RUN_TYPE /= 'RESTART_1' .AND. PRINT_DES_DATA) THEN
         S_TIME = 0.0d0 !TIME
         CALL WRITE_DES_DATA(des_radius, des_pos_new, des_vel_new, des_usr_var)
      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MAKE_ARRAYS_DES
END MODULE MAKE_ARRAYS_DES_MODULE
