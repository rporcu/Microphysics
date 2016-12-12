!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!
!> \mainpage Multiphase Flow with Interphase eXchanges
!!
!! MFIX is a general-purpose computer code developed at the National
!! Energy Technology Laboratory, NETL, for describing the hydrodynamics,
!! heat transfer, and chemical reactions in fluid-solid systems.
!!
!! It has been used for describing bubbling and circulating fluidized
!! beds and spouted beds. MFiX calculations give transient data on the
!! three-dimensional distribution of pressure, velocity, temperature,
!! and species mass fractions. MFiX code is based on a generally
!! accepted set of multiphase flow equations. The code is used as a
!! "test-stand" for testing and developing multiphase flow constitutive
!!  equations.
!!
!! \section Notice
!! Neither the United States Government nor any agency thereof, nor any
!! of their employees, makes any warranty, expressed or implied, or
!! assumes any legal liability or responsibility for the accuracy,
!! completeness, or usefulness of any information, apparatus, product,
!! or process disclosed or represents that its use would not infringe
!! privately owned rights.
!!
!! * MFIX is provided without any user support for applications in the
!!   user's immediate organization. It should not be redistributed in
!!   whole or in part.
!!
!! * The use of MFIX is to be acknowledged in any published paper based
!!   on computations using this software by citing the MFIX theory
!!   manual. Some of the submodels are being developed by researchers
!!   outside of NETL. The use of such submodels is to be acknowledged
!!   by citing the appropriate papers of the developers of the submodels.
!!
!! * The authors would appreciate receiving any reports of bugs or other
!!   difficulties with the software, enhancements to the software, and
!!   accounts of practical applications of this software.
!!
!! \section Disclaimer
!! This report was prepared as an account of work sponsored by an agency
!! of the United States Government. Neither the United States Government
!! nor any agency thereof, nor any of their employees, makes any
!! warranty, express or implied, or assumes any legal liability or
!! responsibility for the accuracy, completeness, or usefulness of any
!! information, apparatus, product, or process disclosed, or represents
!! that its use would not infringe privately owned rights. Reference
!! herein to any specific commercial product, process, or service by
!! trade name, trademark, manufacturer, or otherwise does not
!! necessarily constitute or imply its endorsement, recommendation, or
!! favoring by the United States Government or any agency thereof. The
!! views and opinions of authors expressed herein do not necessarily
!! state or reflect those of the United States Government or any
!! agency thereof.

      PROGRAM MFIX

!-----------------------------------------------
! Modules
!-----------------------------------------------

      use allocate_mod, only: allocate_arrays
      use check_data_20_module, only: check_data_20
      use compar, only: myPE
      use corner_module, only: get_corner_cells
      use des_allocate, only: des_allocate_arrays
      use discretelement, only: pinc, des_rop_s, max_pip
      use error_manager, only: err_msg, finl_err_msg, flush_err_msg, init_err_msg
      use exit_mod, only: mfix_exit
      use fld_const, only: ro_g0
      use funits , only: dmp_log, unit_log
      use geometry, only: flag
      use get_data_module, only: get_data
      use machine, only: wall_time
      use make_arrays_des_module, only: make_arrays_des
      use matrix, only: A_m, b_m
      use param1 , only: undefined, zero
      use read_res1_mod, only: read_res1
      use run    , only: call_usr, dt, dt_min, dt_max, time, nstep, run_type, dem_solids
      use set_bc0_module, only: set_bc0
      use set_bc1_module, only: set_bc1
      use set_bc_dem_module, only: set_bc_dem
      use set_constprop_module, only: set_constprop
      use set_flags_module, only: set_flags1
      use set_fluidbed_p_module, only: set_fluidbed_p
      use set_ic_module, only: set_ic
      use set_ps_module, only: set_ps
      use set_ro_g_module, only: set_ro_g
      use time_cpu, only: cpu_io, cpu_nlog, cpuos, time_nlog, wall0, cpu00
      use time_march_module, only: time_march
      use write_out0_module, only: write_out0
      use write_out1_module, only: write_out1
      use write_out3_module, only: write_out3
      use write_res1_mod, only: write_res1
      use zero_norm_vel_module, only: zero_norm_vel

      use discretelement, only: des_radius, ro_sol, pmass, omoi, des_pos_new, des_vel_new, omega_new, particle_state, pvol
      use discretelement, only: dg_pijk, dg_pijkprv, ighost_updated, neighbor_index, fc, tow, wall_collision_facet_id, pijk
      use discretelement, only: rot_acc_old, des_usr_var_size, particle_phase, ppos, neighbor_index_old
      use discretelement, only: wall_collision_pft, iglobal_id, drag_fc, des_acc_old, nonexistent, do_old, des_usr_var

      IMPLICIT NONE

! Fluid Variables
!---------------------------------------------------------------------//
! Void fraction
      DOUBLE PRECISION, ALLOCATABLE ::  EP_g(:,:,:), EP_go(:,:,:)
! Gas pressure
      DOUBLE PRECISION, ALLOCATABLE ::  P_g(:,:,:), P_gO(:,:,:)
! Gas density
      DOUBLE PRECISION, ALLOCATABLE ::  RO_g(:,:,:), RO_go(:,:,:)
! Macroscopic gas density
      DOUBLE PRECISION, ALLOCATABLE ::  ROP_g(:,:,:), ROP_go(:,:,:)
! x-component of gas velocity
      DOUBLE PRECISION, ALLOCATABLE ::  U_g(:,:,:), U_go(:,:,:)
! y-component of gas velocity
      DOUBLE PRECISION, ALLOCATABLE ::  V_g(:,:,:), V_go(:,:,:)
! z-component of gas velocity
      DOUBLE PRECISION, ALLOCATABLE ::  W_g(:,:,:), W_go(:,:,:)
! Gas viscosity
      DOUBLE PRECISION, ALLOCATABLE ::  MU_g(:,:,:)
! Second coefficient of viscosity
      DOUBLE PRECISION, ALLOCATABLE ::  LAMBDA_G(:,:,:)
! trace of D_g at i, j, k
      DOUBLE PRECISION, ALLOCATABLE ::  trD_g(:,:,:)
! cross terms
      DOUBLE PRECISION, ALLOCATABLE ::  TAU_U_g(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  TAU_V_g(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  TAU_W_g(:,:,:)
! Gas mass fluxes at the east, north, and top faces
      DOUBLE PRECISION, ALLOCATABLE ::  Flux_gE(:,:,:), Flux_gN(:,:,:), Flux_gT(:,:,:)
! macroscopic gas density at east, north, and top faces
      DOUBLE PRECISION, ALLOCATABLE ::  ROP_gE(:,:,:), ROP_gN(:,:,:), ROP_gT(:,:,:)
! Pressure correction equation
      DOUBLE PRECISION, ALLOCATABLE ::  Pp_g(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  d_e(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  d_n(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  d_t(:,:,:)

! drag coefficient between gas phase and discrete particles
      DOUBLE PRECISION, ALLOCATABLE :: f_gds(:,:,:)
! the coefficient add to gas momentum A matrix
      DOUBLE PRECISION, ALLOCATABLE :: drag_am(:,:,:)
! the coefficient add to gas momentum B matrix
      DOUBLE PRECISION, ALLOCATABLE :: drag_bm(:,:,:,:)

!---------------------------------------------------------------------//

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Final value of CPU time.
      DOUBLE PRECISION :: CPU1
! time used for computations.
      DOUBLE PRECISION :: CPUTIME_USED, WALLTIME_USED
! Save TIME in input file for RESTART_2
      DOUBLE PRECISION :: TIME_SAVE
! Temporary storage for DT
      DOUBLE PRECISION :: DT_tmp

      INTEGER :: II, lb, ub


! Invoke MPI initialization routines and get rank info.
      ! CALL PARALLEL_INIT

! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
! queue terminates. timestep at the beginning of execution
      CALL CPU_TIME (CPU00)
      WALL0 = WALL_TIME()

! Read input data, check data, do computations for IC and BC locations
! and flows, and set geometry parameters such as X, X_E, DToDX, etc.
      call get_data(ro_sol,flag)

! Allocate array storage.
      CALL ALLOCATE_ARRAYS(A_m, B_m,ep_g,p_g,ro_g,rop_g,u_g,v_g,w_g,&
         ep_go,p_go,ro_go,rop_go,u_go,v_go,w_go,d_e,d_n,d_t,pp_g,&
         mu_g,lambda_g,trD_g,tau_u_g,tau_v_g,tau_w_g,flux_ge,&
         flux_gn,flux_gt,rop_ge,rop_gn,rop_gt, f_gds, drag_am, drag_bm)

      IF(DEM_SOLIDS) CALL DES_ALLOCATE_ARRAYS
      IF (DEM_SOLIDS) THEN
         PINC(:,:,:) = 0
         DES_ROP_S(:,:,:,:) = ZERO

         lb = 1
         ub = MAX_PIP

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

#ifdef MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

      IF (DT_TMP /= UNDEFINED) THEN
         DT = MAX(DT_MIN,MIN(DT_MAX,DT))
      ELSE
         DT = DT_TMP
      ENDIF

! Set the flags for wall surfaces impermeable and identify flow
! boundaries using FLAG_E, FLAG_N, and FLAG_T
      CALL SET_FLAGS1(flag)

! Calculate cell volumes and face areas
      CALL SET_GEOMETRY1

! Find corner cells and set their face areas to zero
      CALL GET_CORNER_CELLS(flag)

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
         flux_ge, flux_gn, flux_gt, flag)

! Check the field variable data and report errors.
      CALL CHECK_DATA_20(ep_g,p_g,ro_g,rop_g,u_g,v_g,w_g,flag)

      IF(DEM_SOLIDS) CALL MAKE_ARRAYS_DES(ep_g,ro_g,rop_g, &
         flag, pijk, dg_pijk, dg_pijkprv, iglobal_id, particle_state, particle_phase, neighbor_index, neighbor_index_old, &
         des_radius, ro_sol, pvol, pmass, omoi, &
         ppos, des_pos_new, des_vel_new, des_usr_var, omega_new, fc)

! Set the inflow/outflow BCs for DEM solids
      IF(DEM_SOLIDS) CALL SET_BC_DEM(flag)

! Initializations for CPU time calculations in iterate
      CPUOS = 0.
      CALL CPU_TIME (CPU1)
      CPU_NLOG = CPU1
      TIME_NLOG = TIME - DT

! Find the solution of the equations from TIME to TSTOP at
! intervals of DT
      call time_march(u_g, v_g, w_g, u_go, v_go, w_go, &
         p_g, p_go, pp_g, ep_g, ep_go, &
         ro_g, ro_go, rop_g, rop_go, &
         rop_ge, rop_gn, rop_gt, d_e, d_n, d_t, &
         tau_u_g, tau_v_g, tau_w_g,&
         flux_ge, flux_gn, flux_gt, trd_g, lambda_g, mu_g, &
         f_gds, drag_am, drag_bm, flag, &
         pijk, dg_pijk, dg_pijkprv, iglobal_id, particle_state, particle_phase, &
         des_radius, ro_sol, pvol, pmass, omoi, neighbor_index, neighbor_index_old, &
         ppos, des_pos_new, des_vel_new, des_usr_var, omega_new, des_acc_old, rot_acc_old, fc, tow, wall_collision_pft)

! Call user-defined subroutine after time-loop.
      IF (CALL_USR) CALL USR3(u_g, v_g, w_g, p_g)

! Compute the CPU time and write it out in the .OUT file.
      CPUTIME_USED = 0.0d0
      WALLTIME_USED = WALL_TIME() - WALL0
      CALL WRITE_OUT3 (CPUTIME_USED, WALLTIME_USED, CPU_IO)

! Finalize and terminate MPI
      ! call parallel_fin

      CALL FINL_ERR_MSG

      STOP

 1010 FORMAT('Message 1010: Read in data from .RES file for TIME = ',&
         G12.5,/'Time step number (NSTEP) =',I7)

      END PROGRAM MFIX
