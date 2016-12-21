MODULE MAKE_ARRAYS_DES_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MAKE_ARRAYS_DES(ep_g,ro_g,rop_g, flag, vol_surr, &
         pijk, dg_pijk, dg_pijkprv, iglobal_id, particle_state,&
         particle_phase, neighbor_index, neighbor_index_old, &
         des_radius, ro_sol, pvol, pmass, omoi, &
         ppos, des_pos_new, des_vel_new, des_usr_var, omega_new, fc)

      USE comp_mean_fields_module, only: comp_mean_fields
      USE compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
      USE compar, only: iend1, jend1, kend1
      USE compar, only: istart2, jstart2, kstart2
      USE compar, only: numpes, mype
      USE constant, only: pi
      USE desgrid, only: desgrid_pic
      USE discretelement, only: do_nsearch, imax_global_id, pip, particles, max_pip, ighost_cnt, vtp_findex
      USE discretelement, only: entering_ghost, exiting_ghost, nonexistent, normal_ghost
      USE discretelement, only: gener_part_config, print_des_data, s_time
      USE error_manager, only: err_msg, flush_err_msg, init_err_msg, finl_err_msg
      USE functions, only: ip1, jp1, kp1
      USE generate_particles, only: GENERATE_PARTICLE_CONFIG
      USE geometry, only: vol
      USE mpi_funs_des, only: DES_PAR_EXCHANGE
      USE neighbour_module, only: neighbour
      USE param1, only: zero
      USE particles_in_cell_module, only: init_particles_in_cell, particles_in_cell, pic_search
      USE read_par_input_module, only: read_par_input
      USE read_res0_des_module, only: read_res0_des
      USE run, only: run_type, time
      USE set_filter_des_module, only: set_filter_des
      USE set_phase_index_module, only: set_phase_index
      USE stl_preproc_des, only: add_facet
      USE write_des_data_module, only: write_des_data

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer, intent(inout) :: flag &
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)
      DOUBLE PRECISION, INTENT(INOUT) :: vol_surr&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: pvol, pmass, des_radius, ro_sol, omoi
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: fc
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: des_vel_new, des_pos_new, ppos, omega_new, des_usr_var
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(OUT) :: dg_pijk, iglobal_id, dg_pijkprv
      INTEGER, DIMENSION(:), INTENT(OUT) :: neighbor_index, neighbor_index_old
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_phase
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I, J, K, L
      INTEGER :: count
      INTEGER :: I1, I2, J1, J2, K1, K2
      INTEGER :: lcurpar, lpip_all(0:numpes-1), lglobal_id

      CALL INIT_ERR_MSG("MAKE_ARRAYS_DES")

! Check interpolation input.
      CALL SET_FILTER_DES(flag)

      vol_surr(:,:,:) = ZERO

      ! Initialize vol_surr array
      DO K = KSTART2, KEND1
         DO J = JSTART2, JEND1
            DO I = ISTART2, IEND1

               I1 = I
               I2 = ip1(i)
               J1 = J
               J2 = jp1(j)
               K1 = K
               K2 = kp1(k)

               ! Looping over stencil points (node values)
               count = 0
               if(1.eq.flag(i1,j1,k1,1)) count = count + 1
               if(1.eq.flag(i2,j1,k1,1)) count = count + 1
               if(1.eq.flag(i1,j2,k1,1)) count = count + 1
               if(1.eq.flag(i2,j2,k1,1)) count = count + 1
               if(1.eq.flag(i1,j1,k2,1)) count = count + 1
               if(1.eq.flag(i2,j1,k2,1)) count = count + 1
               if(1.eq.flag(i1,j2,k2,1)) count = count + 1
               if(1.eq.flag(i2,j2,k2,1)) count = count + 1
               vol_surr(i,j,k) = dble(count) * vol

            ENDDO
         ENDDO
      ENDDO



! Set the initial particle data.
      IF(RUN_TYPE == 'NEW') THEN
         IF(PARTICLES /= 0) THEN
            IF(GENER_PART_CONFIG) THEN
               CALL GENERATE_PARTICLE_CONFIG(flag, pijk, particle_state, particle_phase, &
                  des_radius, ro_sol, &
                  des_pos_new, des_vel_new, omega_new)

            ELSE
               CALL READ_PAR_INPUT(particle_state, des_radius, ro_sol, des_pos_new, des_vel_new)
            ENDIF
         ENDIF

! Set the global ID for the particles and set the ghost cnt
         ighost_cnt = 0
         lpip_all = 0
         lpip_all(mype) = pip
         ! call global_all_sum(lpip_all)
         lglobal_id = sum(lpip_all(0:mype-1))
         imax_global_id = 0
         do lcurpar  = 1,pip
            lglobal_id = lglobal_id + 1
            iglobal_id(lcurpar) = lglobal_id
            imax_global_id = iglobal_id(pip)
         end do
         ! call global_all_max(imax_global_id)

! Initialize old values
         omega_new(:,:)   = zero

! Read the restart file.
      ELSEIF(RUN_TYPE == 'RESTART_1' .OR. RUN_TYPE == 'RESTART_2') THEN

         CALL READ_RES0_DES(dg_pijk, dg_pijkprv, iglobal_id, particle_state, &
            des_radius, ro_sol, des_usr_var, &
            des_pos_new, des_vel_new, omega_new)
         imax_global_id = maxval(iglobal_id(1:pip))
         ! call global_all_max(imax_global_id)

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
         IF(NONEXISTENT==PARTICLE_STATE(L)) CYCLE
         IF(NORMAL_GHOST==PARTICLE_STATE(L) .OR. ENTERING_GHOST==PARTICLE_STATE(L) .OR. EXITING_GHOST==PARTICLE_STATE(L)) CYCLE
         PVOL(L) = (4.0D0/3.0D0)*PI*DES_RADIUS(L)**3
         PMASS(L) = PVOL(L)*RO_SOL(L)
         OMOI(L) = 2.5D0/(PMASS(L)*DES_RADIUS(L)**2) !ONE OVER MOI
      ENDDO

      CALL SET_PHASE_INDEX(particle_phase,des_radius,ro_sol,particle_state)
      CALL INIT_PARTICLES_IN_CELL(pijk, particle_state, dg_pijk, dg_pijkprv, &
         des_usr_var, des_pos_new, des_vel_new, omega_new, fc)

! do_nsearch should be set before calling particle in cell
      DO_NSEARCH =.TRUE.
! Bin the particles to the DES grid.
      CALL DESGRID_PIC(PLOCATE=.TRUE., dg_pijkprv=dg_pijkprv, dg_pijk=dg_pijk, &
         des_pos_new=des_pos_new, particle_state=particle_state)
      CALL DES_PAR_EXCHANGE(pijk, particle_state, dg_pijk, dg_pijkprv, &
         des_usr_var, des_pos_new, des_vel_new, omega_new, fc)
      CALL PARTICLES_IN_CELL(pijk, iglobal_id, particle_state, des_pos_new, des_vel_new, des_radius, des_usr_var)

      CALL NEIGHBOUR(dg_pijk, particle_state, des_radius, des_pos_new, ppos, neighbor_index, neighbor_index_old)

! Calculate mean fields using either interpolation or cell averaging.
      CALL COMP_MEAN_FIELDS(ep_g,ro_g,rop_g,pijk,particle_state,particle_phase,pmass,pvol, &
         des_pos_new,des_vel_new,des_radius,des_usr_var,flag,vol_surr,iglobal_id)

      IF(RUN_TYPE /= 'RESTART_1' .AND. PRINT_DES_DATA) THEN
         S_TIME = TIME
         CALL WRITE_DES_DATA(des_radius, des_pos_new, des_vel_new, des_usr_var)
      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MAKE_ARRAYS_DES
END MODULE MAKE_ARRAYS_DES_MODULE
