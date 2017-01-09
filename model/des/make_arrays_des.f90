!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine make_arrays_des(max_pip, particle_state, particle_phase, &
      des_radius,  ro_sol, pvol, pmass, omoi, des_pos_new, des_vel_new,&
      des_usr_var, omega_new, fc, tow) &
   bind(C, name="mfix_make_arrays_des")

! Module procedures .................................................//
      use read_par_input_module, only: read_par_input
      use error_manager, only: init_err_msg, flush_err_msg

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

! Global data .......................................................//
      use discretelement, only: particles
      use discretelement, only: vtp_findex
      use error_manager, only: err_msg

      use discretelement, only: do_nsearch
      use discretelement, only: nonexistent, normal_ghost
      use discretelement, only: entering_ghost, exiting_ghost
      use constant, only: pi
      use param1, only: zero

      use run, only: run_type
      use set_phase_index_module, only: set_phase_index

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: max_pip

      real(c_real)  , intent(  out) :: pvol(max_pip)
      real(c_real)  , intent(  out) :: pmass(max_pip)
      real(c_real)  , intent(  out) :: des_radius(max_pip)
      real(c_real)  , intent(  out) :: ro_sol(max_pip)
      real(c_real)  , intent(  out) :: omoi(max_pip)

      real(c_real)  , intent(  out) :: des_vel_new(max_pip,3)
      real(c_real)  , intent(  out) :: des_pos_new(max_pip,3)
      real(c_real)  , intent(  out) :: omega_new(max_pip,3)
      real(c_real)  , intent(  out) :: des_usr_var(max_pip,1)
      integer(c_int), intent(  out) :: particle_state(max_pip)
      integer(c_int), intent(  out) :: particle_phase(max_pip)

      real(c_real), intent(  out) :: fc(max_pip,3)
      real(c_real), intent(  out) :: tow(max_pip,3)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: L


      CALL INIT_ERR_MSG("MAKE_ARRAYS_DES")

! Particle state flag
      particle_state(:) = nonexistent
      particle_phase(:) = nonexistent

! Physical properties:
      des_radius(:) = zero
      ro_sol(:) = zero
      pvol(:) = zero
      pmass(:) = zero
      omoi(:) = zero

! Particle position, velocity, etc
      des_pos_new(:,:) = zero
      des_vel_new(:,:) = zero
      omega_new(:,:) = zero

! Translation and rotational forces
      fc(:,:) = zero
      tow(:,:) = zero

! Initializing user defined array
      des_usr_var(:,:) = zero

! Set the initial particle data.
      IF(RUN_TYPE == 'NEW') THEN
         IF(PARTICLES /= 0) CALL READ_PAR_INPUT(particle_state, &
            des_radius, ro_sol, des_pos_new, des_vel_new)

! Initialize old values
         omega_new(:,:)   = zero

! Read the restart file.
      !ELSEIF(RUN_TYPE == 'RESTART_1' .OR. RUN_TYPE == 'RESTART_2') THEN

      ELSE

         WRITE(ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1100 FORMAT('Error 1100: Unsupported RUN_TYPE for DES.')

      ENDIF

      IF(RUN_TYPE == 'RESTART_2') vtp_findex = 0

      ! Set additional particle properties now that the particles
      ! have been identified
      DO L = 1, MAX_PIP

         ! Skip 'empty' locations when populating the particle property arrays.
         if(   nonexistent == particle_state(l) .or. &
              normal_ghost == particle_state(l) .or. &
            entering_ghost == particle_state(l) .or. &
             exiting_ghost == particle_state(l)) cycle

         pvol(l) = (4.0d0/3.0d0)*pi*des_radius(l)**3
         pmass(l) = pvol(l)*ro_sol(l)
         omoi(l) = 2.5d0/(pmass(l)*des_radius(l)**2) !one over moi
      ENDDO

      CALL SET_PHASE_INDEX(max_pip, particle_phase, des_radius,&
         ro_sol, particle_state)

! do_nsearch should be set before calling particle in cell
      DO_NSEARCH =.TRUE.

!      CALL NEIGHBOUR(  particle_state, des_radius, des_pos_new)

   end subroutine make_arrays_des



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: MAKE_ARRAYS_DES                                        !
!                                                                      !
!  Purpose: DES - allocating DES arrays                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine mfix_write_des_data(max_pip, des_radius, des_pos_new, &
      des_vel_new, des_usr_var) &
      bind(C, name="mfix_write_des_data")

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use discretelement, only: print_des_data, s_time
      use error_manager        , only: finl_err_msg
      use run                  , only: run_type
      use write_des_data_module, only: write_des_data

      real(c_real), intent(in   ) :: des_radius (max_pip)
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)
      real(c_real), intent(in   ) :: des_usr_var(max_pip,1)

      IF (RUN_TYPE /= 'RESTART_1' .AND. PRINT_DES_DATA) THEN
         S_TIME = 0.0d0 !TIME
         CALL WRITE_DES_DATA(max_pip, des_radius, des_pos_new,&
            des_vel_new, des_usr_var)
      ENDIF

      CALL FINL_ERR_MSG

      end subroutine mfix_write_des_data
