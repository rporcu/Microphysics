module des_time_march_module

   use particles_in_cell_module, only: particles_in_cell
   use mass_outflow_dem_module, only: mass_outflow_dem

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_MARCH                                       !
!     Author: Jay Boyalakuntla                        Date: 21-Jun-04  !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_TIME_MARCH(ep_g, p_g, u_g, v_g, w_g, ro_g, rop_g, mu_g, &
         pijk, dg_pijk, dg_pijkprv, iglobal_id, particle_state, particle_phase, neighbor_index, neighbor_index_old, &
         des_radius, ro_sol, pvol, pmass, omoi, des_usr_var, &
         ppos, des_pos_new, des_vel_new, omega_new, des_acc_old, rot_acc_old, fc, tow, wall_collision_pft)

      USE neighbour_module, only: neighbour
      use calc_drag_des_module, only: calc_drag_des
      use calc_epg_des_module, only: calc_epg_des
      use calc_force_dem_module, only: calc_force_dem
      use calc_pg_grad_module, only: calc_pg_grad
      use cfnewvalues_module, only: cfnewvalues
      use comp_mean_fields_module, only: comp_mean_fields
      use compar, only: iend3, jend3, kend3
      use compar, only: istart3, jstart3, kstart3
      use compar, only: numpes
      use des_allocate, only: particle_grow
      use des_bc, only: DEM_BCMI, DEM_BCMO
      use desgrid, only: desgrid_pic
      use discretelement, only: des_continuum_coupled, des_explicitly_coupled, des_periodic_walls, dtsolid, ighost_cnt
      use discretelement, only: drag_fc
      use discretelement, only: pip, s_time, do_nsearch, neighbor_search_n
      use drag_gs_des1_module, only: drag_gs_des1
      use error_manager, only: err_msg, init_err_msg, finl_err_msg, ival, flush_err_msg
      use machine, only:  wall_time
      use mass_inflow_dem_module, only: mass_inflow_dem
      use mpi_funs_des, only: DES_PAR_EXCHANGE
      use output_manager_module, only: output_manager
      use param1, only: zero
      use run, only: CALL_USR
      use run, only: NSTEP
      use run, only: TIME, TSTOP, DT

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:), INTENT(INOUT) :: wall_collision_pft
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: pvol, pmass, des_radius, ro_sol, omoi
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: des_acc_old, rot_acc_old, fc, tow
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: des_vel_new, des_pos_new, ppos, omega_new, des_usr_var
      INTEGER, DIMENSION(:), INTENT(INOUT) :: NEIGHBOR_INDEX, NEIGHBOR_INDEX_OLD
      INTEGER(KIND=1), DIMENSION(:), INTENT(INOUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(OUT) :: dg_pijk, iglobal_id
      INTEGER, DIMENSION(:), INTENT(OUT) :: dg_pijkprv
      INTEGER, DIMENSION(:), INTENT(OUT) :: particle_phase
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: pijk

!------------------------------------------------
! Local variables
!------------------------------------------------
! Total number of particles
      INTEGER, SAVE :: NP=0

! time step loop counter index
      INTEGER :: NN

! loop counter index for any initial particle settling incoupled cases
      INTEGER :: FACTOR

! Temporary variables when des_continuum_coupled is T to track
! changes in solid time step
      DOUBLE PRECISION :: TMP_DTS, DTSOLID_TMP

! Numbers to calculate wall time spent in DEM calculations.
      DOUBLE PRECISION :: TMP_WALL
! Pressure gradient
      DOUBLE PRECISION, allocatable :: gradPg(:,:,:,:)
!.......................................................................!

      allocate (gradPg(istart3:iend3, jstart3:jend3, kstart3:kend3,3) )

! In case of restarts assign S_TIME from MFIX TIME
      S_TIME = TIME
      TMP_DTS = ZERO
      DTSOLID_TMP = ZERO
      TMP_WALL = WALL_TIME()

! Initialize time stepping variables for coupled gas/solids simulations.
      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DT.GE.DTSOLID) THEN
            FACTOR = CEILING(real(DT/DTSOLID))
         ELSE
            FACTOR = 1
            DTSOLID_TMP = DTSOLID
            DTSOLID = DT
         ENDIF

! Initialize time stepping variable for pure granular simulations.
      ELSE
         FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID))
         DT = DTSOLID
         CALL OUTPUT_MANAGER(ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g, &
         iglobal_id, particle_state, des_radius, ro_sol, des_pos_new, des_vel_new, des_usr_var, omega_new, &
         .FALSE., .FALSE.)
      ENDIF   ! end if/else (des_continuum_coupled)

      NP = PIP - IGHOST_CNT
      ! CALL GLOBAL_ALL_SUM(NP)

      IF(DES_CONTINUUM_COUPLED) THEN
         WRITE(ERR_MSG, 1000) trim(iVal(factor)), trim(iVAL(NP))
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)
      ELSE
         WRITE(ERR_MSG, 1100) TIME, DTSOLID, trim(iVal(factor))
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)
      ENDIF
 1000 FORMAT(/'DEM NITs: ',A,3x,'Total PIP: ', A)
 1100 FORMAT(/'Time: ',g12.5,3x,'DT: ',g12.5,3x,'DEM NITs: ',A)

      IF(CALL_USR) CALL USR0_DES

      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DES_EXPLICITLY_COUPLED) THEN
            CALL DRAG_GS_DES1(ep_g, u_g, v_g, w_g, ro_g, mu_g, gradPg, pijk, particle_state, pvol, des_vel_new, fc, des_radius)
         ENDIF
         CALL CALC_PG_GRAD(p_g, gradPg, pijk)
      ENDIF


! Main DEM time loop
!----------------------------------------------------------------->>>
      DO NN = 1, FACTOR

         IF(DES_CONTINUUM_COUPLED) THEN
! If the current time in the discrete loop exceeds the current time in
! the continuum simulation, exit the discrete loop
            IF(S_TIME.GE.(TIME+DT)) EXIT
! If next time step in the discrete loop will exceed the current time
! in the continuum simulation, modify the discrete time step so final
! time will match
            IF((S_TIME+DTSOLID).GT.(TIME+DT)) THEN
               TMP_DTS = DTSOLID
               DTSOLID = TIME + DT - S_TIME
            ENDIF
         ENDIF

! Calculate forces acting on particles (collisions, drag, etc).
         CALL CALC_FORCE_DEM(particle_phase, particle_state, &
            des_radius, des_pos_new, des_vel_new, omega_new, fc, tow, wall_collision_pft)
! Calculate or distribute fluid-particle drag force.
         CALL CALC_DRAG_DES(ep_g,u_g,v_g,w_g,ro_g,mu_g,gradPg,pijk,particle_state,fc,drag_fc,pvol, &
            des_pos_new,des_vel_new,des_radius)

! Call user functions.
         IF(CALL_USR) CALL USR1_DES
! Update position and velocities
         CALL CFNEWVALUES(particle_state, des_radius, omoi, ppos, des_pos_new, des_vel_new, omega_new, fc, tow, &
            des_acc_old, rot_acc_old)

! Set DO_NSEARCH before calling DES_PAR_EXCHANGE.
         DO_NSEARCH = (NN == 1 .OR. MOD(NN,NEIGHBOR_SEARCH_N) == 0)

! Add/Remove particles to the system via flow BCs.
         IF(DEM_BCMI > 0) CALL MASS_INFLOW_DEM(PIJK, particle_phase, DG_PIJK, iglobal_id, PARTICLE_STATE, &
            DES_RADIUS, OMOI, PMASS, PVOL, RO_SOL, &
            DES_VEL_NEW, DES_POS_NEW, PPOS, OMEGA_NEW, DRAG_FC)
         IF(DEM_BCMO > 0) CALL MASS_OUTFLOW_DEM(DO_NSEARCH, particle_phase, iglobal_id, particle_state, &
            des_radius, omoi, pmass, pvol, ro_sol, &
            des_vel_new, des_pos_new, ppos, omega_new, fc, tow)

! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
         IF (DO_NSEARCH .OR. (numPEs>1) .OR. DES_PERIODIC_WALLS) THEN
            CALL DESGRID_PIC(.TRUE., dg_pijkprv=dg_pijkprv, dg_pijk=dg_pijk, &
               des_pos_new=des_pos_new, particle_state=particle_state)
            CALL DES_PAR_EXCHANGE(des_pos_new, dg_pijk, dg_pijkprv, particle_state)
         ENDIF

         IF(DO_NSEARCH) CALL NEIGHBOUR(dg_pijk, particle_state, des_radius, des_pos_new, ppos, neighbor_index, neighbor_index_old)

! Explicitly coupled simulations do not need to rebin particles to
! the fluid grid every time step. However, this implies that the
! fluid cell information and interpolation weights become stale.
         IF(DES_CONTINUUM_COUPLED .AND. &
            .NOT.DES_EXPLICITLY_COUPLED) THEN
! Bin particles to fluid grid.
            CALL PARTICLES_IN_CELL(pijk, iglobal_id, particle_state, des_pos_new, des_vel_new, des_radius, des_usr_var)
! Calculate mean fields (EPg).
            CALL COMP_MEAN_FIELDS(ep_g,ro_g,rop_g,pijk,particle_phase,pmass,pvol, &
               des_pos_new,des_vel_new,des_radius,des_usr_var,iglobal_id)
         ENDIF

! Update time to reflect changes
         S_TIME = S_TIME + DTSOLID

! The following section targets data writes for DEM only cases:
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
! Keep track of TIME and number of steps for DEM simulations
            TIME = S_TIME
            NSTEP = NSTEP + 1
! Call the output manager to write RES data.
            CALL OUTPUT_MANAGER(ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g, &
               iglobal_id, particle_state, des_radius, ro_sol, des_pos_new, des_vel_new, des_usr_var, omega_new, &
               .FALSE., .FALSE.)
         ENDIF  ! end if (.not.des_continuum_coupled)

         IF(CALL_USR) CALL USR2_DES

      ENDDO ! end do NN = 1, FACTOR

! END DEM time loop
!-----------------------------------------------------------------<<<

      IF(CALL_USR) CALL USR3_DES

       CALL CALC_EPG_DES(ep_g,ro_g,rop_g,des_pos_new,des_vel_new,des_radius,des_usr_var,iglobal_id)

! When coupled, and if needed, reset the discrete time step accordingly
      IF(DT.LT.DTSOLID_TMP) THEN
         DTSOLID = DTSOLID_TMP
      ENDIF

      IF(TMP_DTS.NE.ZERO) THEN
         DTSOLID = TMP_DTS
         TMP_DTS = ZERO
      ENDIF

      deallocate(gradPg)

      IF(.NOT.DES_CONTINUUM_COUPLED)THEN
         WRITE(ERR_MSG,"('<---------- END DES_TIME_MARCH ----------')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ELSE
         ! call send_recv(ep_g,2)
         ! call send_recv(rop_g,2)

         TMP_WALL = WALL_TIME() - TMP_WALL
         IF(TMP_WALL > 1.0d-10) THEN
            WRITE(ERR_MSG, 9000) trim(iVal(dble(FACTOR)/TMP_WALL))
         ELSE
            WRITE(ERR_MSG, 9000) '+Inf'
         ENDIF
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)

 9000 FORMAT('    NITs/SEC = ',A)

      ENDIF

      END SUBROUTINE DES_TIME_MARCH

end module des_time_march_module
