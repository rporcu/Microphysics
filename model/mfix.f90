!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine MFIX(u_g, v_g, w_g, u_go, v_go, w_go, &
                p_g, p_go, pp_g, ep_g, ep_go, &
                ro_g, ro_go, rop_g, rop_go, &
                rop_ge, rop_gn, rop_gt, &
                d_e, d_n, d_t, &
                tau_u_g ,tau_v_g, tau_w_g,&
                flux_ge, flux_gn, flux_gt, &
                trD_g, lambda_g, mu_g, &
                f_gds, A_m, b_m, &
                drag_am, drag_bm, pinc, &
                flag, vol_surr, des_rop_s, &
                pijk, dg_pijk, dg_pijkprv, iglobal_id, &
                particle_state, particle_phase, des_radius, ro_sol, pvol, pmass, &
                omoi, ppos, des_pos_new, des_vel_new, des_usr_var, omega_new, des_acc_old,&
                rot_acc_old, fc, tow, wall_collision_pft)

!-----------------------------------------------
! Modules
!-----------------------------------------------

      use check_data_20_module, only: check_data_20
      use compar, only: myPE, istart3, iend3, jstart3, jend3, kstart3, kend3
      use constant, only: mmax
      use corner_module, only: get_corner_cells
      use des_allocate, only: des_allocate_arrays
      use discretelement, only: max_pip
      use error_manager, only: err_msg, finl_err_msg, flush_err_msg, init_err_msg
      use exit_mod, only: mfix_exit
      use fld_const, only: ro_g0, mu_g0
      use funits , only: dmp_log, unit_log
      use geometry, only: dx, dy, dz, ayz, axy, axz, vol, flag_mod
      use set_domain_module, only: set_domain
      use machine, only: wall_time
      use make_arrays_des_module, only: make_arrays_des
      use param1 , only: undefined, zero
      use read_res1_mod, only: read_res1
      use run, only: call_usr, run_type, dem_solids, nstep
      use run, only: dt, dt_min, dt_max, time, tstop, use_dt_prev, dt_prev
      use set_bc0_module, only: set_bc0
      use set_bc1_module, only: set_bc1
      use set_bc_dem_module, only: set_bc_dem
      use set_constprop_module, only: set_constprop
      use set_flags_module, only: set_flags1
      use set_fluidbed_p_module, only: set_fluidbed_p
      use set_ic_module, only: set_ic
      use set_ps_module, only: set_ps
      use set_ro_g_module, only: set_ro_g
      use write_out0_module, only: write_out0
      use write_out1_module, only: write_out1
      use write_out3_module, only: write_out3
      use write_res1_mod, only: write_res1
      use zero_norm_vel_module, only: zero_norm_vel
      use output_manager_module, only: init_output_vars
      use calc_coeff_module, only: calc_coeff

      use discretelement, only: neighbor_index, wall_collision_facet_id
      use discretelement, only: des_usr_var_size, neighbor_index_old
      use discretelement, only: drag_fc, nonexistent, do_old, ighost_updated

      IMPLICIT NONE

      integer         , intent(inout) :: flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)
      double precision, intent(inout) :: vol_surr&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: A_m&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,7)
      double precision, intent(inout) :: b_m&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: ep_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: ep_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: p_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: p_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: ro_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: ro_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: u_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: u_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: v_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: v_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: w_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: w_go&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: pp_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: d_e&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: d_t&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: d_n&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: mu_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: lambda_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: trD_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: tau_u_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: tau_v_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: tau_w_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: flux_gE&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: flux_gN&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: flux_gT&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: rop_gE&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_gN&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_gT&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: f_gds&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: drag_am&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: drag_bm&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,3)

      integer, intent(inout) :: pinc&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      double precision, intent(inout) :: des_rop_s&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,mmax)

      integer, intent(inout) :: pijk(max_pip,3)
      integer, intent(inout) :: dg_pijk(max_pip)
      integer, intent(inout) :: dg_pijkprv(max_pip)
      integer, intent(inout) :: iglobal_id(max_pip)
      integer, intent(inout) :: particle_state(max_pip)
      integer, intent(inout) :: particle_phase(max_pip)

      double precision, intent(inout) :: des_radius(max_pip)
      double precision, intent(inout) :: ro_sol(max_pip)
      double precision, intent(inout) :: pvol(max_pip)
      double precision, intent(inout) :: pmass(max_pip)
      double precision, intent(inout) :: omoi(max_pip)

      double precision, intent(inout) :: ppos(max_pip,3)
      double precision, intent(inout) :: des_pos_new(max_pip,3)
      double precision, intent(inout) :: des_vel_new(max_pip,3)
      double precision, intent(inout) :: des_usr_var(max_pip,1)
      double precision, intent(inout) :: omega_new(max_pip,3)

      double precision, intent(inout) :: des_acc_old(max_pip,3)
      double precision, intent(inout) :: rot_acc_old(max_pip,3)
      double precision, intent(inout) :: fc(max_pip,3)
      double precision, intent(inout) :: tow(max_pip,3)

      double precision, intent(inout) :: wall_collision_pft(3,8,max_pip)


!---------------------------------------------------------------------//
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Save TIME in input file for RESTART_2
      DOUBLE PRECISION :: TIME_SAVE
! Temporary storage for DT
      DOUBLE PRECISION :: DT_tmp

      INTEGER :: II, lb, ub

!---------------------------------------------------------------------//
      flag_mod = flag
!-----------------------------------------------

      ! This is now called from main.cpp
      ! call set_domain(flag)

      IF(DEM_SOLIDS) CALL DES_ALLOCATE_ARRAYS

      IF (DEM_SOLIDS) THEN
         DES_ROP_S(:,:,:,:) = ZERO

         lb = 1
         ub = MAX_PIP

         write(*,*)"....",ubound(iglobal_id)
         write(*,*)"....",lbound(iglobal_id)
         IGLOBAL_ID(LB:UB) = 0
         PARTICLE_STATE(LB:UB) = NONEXISTENT

! Physical properties:
         DES_RADIUS(LB:UB) = ZERO
         RO_Sol(LB:UB) = ZERO
         PVOL(LB:UB) = ZERO
         PMASS(LB:UB) = ZERO
         OMOI(LB:UB) = ZERO

! Particle position, velocity, etc
         DES_POS_NEW(LB:UB,:) = ZERO
         DES_VEL_NEW(LB:UB,:) = ZERO
         OMEGA_NEW(LB:UB,:) = ZERO

! Particle state flag
         DO II = LB, UB
            particle_state(II) = nonexistent
         ENDDO
         NEIGHBOR_INDEX(:) = 0

! DES grid bin information
         DG_PIJK(LB:UB) = -1
         DG_PIJKPRV(LB:UB) = -1
         IGHOST_UPDATED(LB:UB) = .false.

! Fluid cell bin information
         PIJK(LB:UB,:) = 0

! Translation and rotational forces
         FC(LB:UB,:) = ZERO
         TOW(LB:UB,:) = ZERO

! Collision data
         WALL_COLLISION_FACET_ID(:,LB:UB) = -1
         WALL_COLLISION_PFT(:,:,LB:UB) = ZERO

! Initializing user defined array
         IF(DES_USR_VAR_SIZE > 0) &
            DES_USR_VAR(LB:UB,:) = ZERO

! Particle center drag coefficient and explicit drag force
         DRAG_FC(LB:UB,:) = ZERO

! Higher order time integration variables.
         IF (DO_OLD) THEN
            DES_ACC_OLD(LB:UB,:) = ZERO
            ROT_ACC_OLD(LB:UB,:) = ZERO
         ENDIF
      ENDIF

! Write the initial part of the standard output file
      CALL WRITE_OUT0
!     CALL WRITE_FLAGS

! Write the initial part of the special output file(s)
      CALL WRITE_USR0

      CALL INIT_ERR_MSG('MFIX')

      DT_TMP = DT
      SELECT CASE (TRIM(RUN_TYPE))

      CASE ('NEW')
! Write the initial part of the restart files
         CALL WRITE_RES0

       CASE ('RESTART_1')
! Read the time-dependent part of the restart file
         CALL READ_RES1(ep_g,p_g,ro_g,u_g,v_g,w_g,rop_g)
         WRITE(ERR_MSG, 1010) TIME, NSTEP
         CALL FLUSH_ERR_MSG()

      CASE ('RESTART_2')
         TIME_SAVE = TIME

         CALL READ_RES1(ep_g,p_g,ro_g,u_g,v_g,w_g,rop_g)
         TIME = TIME_SAVE

         WRITE(ERR_MSG, 1010) TIME, NSTEP
         CALL FLUSH_ERR_MSG()

         CALL WRITE_RES0

! Writing the RES1 can only be done here when re-indexing is turned off
! This will be done after the cell re-indexing is done later in this file.
! This allows restarting independently of the re-indexing setting between
! the previous and current run.
         CALL WRITE_RES1(ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g)

      CASE DEFAULT
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
            ' MFIX: Do not know how to process'
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' RUN_TYPE in data file'
         call mfix_exit(myPE)

      END SELECT

      IF (DT_TMP /= UNDEFINED) THEN
         DT = MAX(DT_MIN,MIN(DT_MAX,DT))
      ELSE
         DT = DT_TMP
      ENDIF

! Set the flags for wall surfaces impermeable and identify flow
! boundaries using FLAG_E, FLAG_N, and FLAG_T
      CALL SET_FLAGS1(flag)
      flag_mod = flag

! Calculate cell volumes and face areas
      VOL = DX*DY*DZ
      AYZ = DY*DZ
      AXY = DX*DY
      AXZ = DX*DZ

! Find corner cells and set their face areas to zero
      CALL GET_CORNER_CELLS(flag)
      flag_mod = flag

! Set constant physical properties
      CALL SET_CONSTPROP(ro_g, lambda_g, mu_g, flag)

! Set initial conditions
      CALL SET_IC(ep_g, p_g, u_g, v_g, w_g, flag)

! Set point sources.
      CALL SET_PS(flag)

! Set boundary conditions
      CALL ZERO_NORM_VEL(u_g,v_g,w_g, flag)
      CALL SET_BC0(p_g,ep_g,u_g,v_g,w_g,ro_g0,flag)

! Set the pressure field for a fluidized bed
      IF (RUN_TYPE == 'NEW') CALL SET_FLUIDBED_P(p_g, ep_g, flag)

! Initialize densities.
      IF (RUN_TYPE == 'NEW') CALL SET_RO_G(ro_g,rop_g,p_g,ep_g,flag)

! Initialize time dependent boundary conditions
      CALL SET_BC1(p_g, ep_g, ro_g, rop_g, u_g, v_g, w_g, &
                   des_rop_s, flux_ge, flux_gn, flux_gt, flag)

! Check the field variable data and report errors.
      CALL CHECK_DATA_20(ep_g,p_g,ro_g,rop_g,u_g,v_g,w_g,flag)

      IF(DEM_SOLIDS) CALL MAKE_ARRAYS_DES(ep_g,ro_g,rop_g, &
         flag, vol_surr, pijk, dg_pijk, dg_pijkprv, iglobal_id, &
         particle_state, particle_phase, neighbor_index, neighbor_index_old, &
         des_radius, des_rop_s, ro_sol, pvol, pmass, omoi, &
         ppos, des_pos_new, des_vel_new, des_usr_var, omega_new, fc, pinc)

! Set the inflow/outflow BCs for DEM solids
      IF(DEM_SOLIDS) CALL SET_BC_DEM(flag)


! ######################## Moved here from time march

      CALL INIT_OUTPUT_VARS

! Parse residual strings
      CALL PARSE_RESID_STRING ()

! Call user-defined subroutine to set constants, check data, etc.
      IF (CALL_USR) CALL USR0

! Calculate all the coefficients once before entering the time loop
      CALL CALC_COEFF(flag, 2, ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, mu_g, &
         f_gds, drag_am, drag_bm, pijk, particle_phase, particle_state, &
         pvol, des_pos_new, des_vel_new, des_radius, des_rop_s, pinc)

      IF(MU_g0 == UNDEFINED) CALL CALC_MU_G(lambda_g,mu_g,mu_g0)

! Remove undefined values at wall cells for scalars
      where(rop_g == undefined) rop_g = 0.0

      CALL FINL_ERR_MSG

 1010 FORMAT('Message 1010: Read in data from .RES file for TIME = ',&
         G12.5,/'Time step number (NSTEP) =',I7)

      END subroutine MFIX
