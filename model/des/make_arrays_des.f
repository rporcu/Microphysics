!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MAKE_ARRAYS_DES(ep_g,ro_g,rop_g)

      USE comp_mean_fields_module, only: comp_mean_fields
      USE compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
      USE compar, only: iend1, jend1, kend1
      USE compar, only: istart2, jstart2, kstart2
      USE compar, only: numpes, mype
      USE constant, only: pi
      USE desgrid, only: desgrid_pic
      USE discretelement, only: entering_ghost, exiting_ghost, nonexistent, particle_state, normal_ghost, pijk
      USE discretelement, only: gener_part_config, print_des_data, s_time, iglobal_id, pvol, pmass, des_radius, ro_sol
      USE discretelement, only: omega_new, do_nsearch, imax_global_id, pip, particles, max_pip, ighost_cnt, omoi, vtp_findex
      USE discretelement, only: des_pos_new, des_vel_new
      USE error_manager, only: err_msg, flush_err_msg, init_err_msg, finl_err_msg
      USE functions, only: ip1, jp1, kp1, fluid_at
      USE generate_particles, only: GENERATE_PARTICLE_CONFIG
      USE geometry, only: vol_surr, vol
      USE mpi_funs_des, only: DES_PAR_EXCHANGE
      USE param1, only: zero
      USE particles_in_cell_module, only: init_particles_in_cell, particles_in_cell, pic_search
      USE run, only: run_type, time
      USE set_phase_index_module, only: set_phase_index
      USE stl_preproc_des, only: add_facet

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I, J, K, L
      INTEGER :: count
      INTEGER :: I1, I2, J1, J2, K1, K2
      INTEGER :: lcurpar, lpip_all(0:numpes-1), lglobal_id

      CALL INIT_ERR_MSG("MAKE_ARRAYS_DES")

! Check interpolation input.
      CALL SET_FILTER_DES

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
               if(fluid_at(i1,j1,k1)) count = count + 1
               if(fluid_at(i2,j1,k1)) count = count + 1
               if(fluid_at(i1,j2,k1)) count = count + 1
               if(fluid_at(i2,j2,k1)) count = count + 1
               if(fluid_at(i1,j1,k2)) count = count + 1
               if(fluid_at(i2,j1,k2)) count = count + 1
               if(fluid_at(i1,j2,k2)) count = count + 1
               if(fluid_at(i2,j2,k2)) count = count + 1
               vol_surr(i,j,k) = dble(count) * vol

            ENDDO
         ENDDO
      ENDDO



! Set the initial particle data.
      IF(RUN_TYPE == 'NEW') THEN
         IF(PARTICLES /= 0) THEN
            IF(GENER_PART_CONFIG) THEN
               CALL GENERATE_PARTICLE_CONFIG
            ELSE
               CALL READ_PAR_INPUT
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

         CALL READ_RES0_DES
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

      CALL SET_PHASE_INDEX(pijk)
      CALL INIT_PARTICLES_IN_CELL(pijk, particle_state, des_pos_new)

! do_nsearch should be set before calling particle in cell
      DO_NSEARCH =.TRUE.
! Bin the particles to the DES grid.
      CALL DESGRID_PIC(PLOCATE=.TRUE.)
      CALL DES_PAR_EXCHANGE
      CALL PARTICLES_IN_CELL(pijk, iglobal_id, particle_state, des_pos_new, des_vel_new)

      CALL NEIGHBOUR

! Calculate mean fields using either interpolation or cell averaging.
      CALL COMP_MEAN_FIELDS(ep_g,ro_g,rop_g,pijk,pmass,pvol,des_pos_new,des_vel_new)

      IF(RUN_TYPE /= 'RESTART_1' .AND. PRINT_DES_DATA) THEN
         S_TIME = TIME
         CALL WRITE_DES_DATA
      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MAKE_ARRAYS_DES
