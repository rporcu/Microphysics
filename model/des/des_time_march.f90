module des_time_march_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

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
             particle_state, particle_phase, &
             des_radius,  ro_sol, pvol, pmass, omoi, des_usr_var, &
             des_pos_new, des_vel_new, omega_new, des_acc_old, rot_acc_old, &
             drag_fc, fc, tow, pairs, pair_count, flag, &
             time, dt, dx, dy, dz, nstep) &
             bind(C, name="mfix_des_time_march")

      USE calc_collision_wall, only: calc_dem_force_with_wall_stl
      use calc_drag_des_module, only: calc_drag_des
      use calc_force_dem_module, only: calc_force_dem
      use calc_pg_grad_module, only: calc_pg_grad
      use cfnewvalues_module, only: cfnewvalues
      use comp_mean_fields_module, only: comp_mean_fields
      use compar, only: iend3, jend3, kend3
      use compar, only: istart3, jstart3, kstart3
      use discretelement, only: des_continuum_coupled, des_explicitly_coupled
      use discretelement, only: dtsolid
      use discretelement, only: pip, s_time, do_nsearch
      use drag_gs_des1_module, only: drag_gs_des
      use error_manager, only: err_msg, init_err_msg, finl_err_msg, ival, flush_err_msg
      use machine, only:  wall_time
      use output_manager_module, only: output_manager
      use param1, only: zero
      use run, only: CALL_USR
      use run, only: TSTOP
      use discretelement, only: max_pip

      IMPLICIT NONE

      real(c_real), intent(inout) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(inout) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), intent(in   ) :: mu_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      integer(c_int), intent(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)

      real(c_real)  , intent(inout) :: time, dt
      real(c_real)  , intent(in   ) :: dx, dy, dz
      integer(c_int), intent(inout) :: nstep

      real(c_real), intent(inout) :: pvol(max_pip)
      real(c_real), intent(inout) :: pmass(max_pip)
      real(c_real), intent(inout) :: des_radius(max_pip)
      real(c_real), intent(inout) :: ro_sol(max_pip)
      real(c_real), intent(inout) :: omoi(max_pip)

      real(c_real), intent(inout) :: des_pos_new(max_pip,3)
      real(c_real), intent(inout) :: des_vel_new(max_pip,3)

      real(c_real), intent(inout) :: des_acc_old(max_pip,3)
      real(c_real), intent(inout) :: rot_acc_old(max_pip,3)
      real(c_real), intent(inout) :: drag_fc(max_pip,3)
      real(c_real), intent(inout) :: fc(max_pip,3)
      real(c_real), intent(inout) :: tow(max_pip,3)

      real(c_real), intent(inout) :: omega_new(max_pip,3)
      real(c_real), intent(inout) :: des_usr_var(max_pip,1)

      integer(c_int), intent(inout) :: particle_state(max_pip)
      integer(c_int), intent(inout) :: particle_phase(max_pip)

      integer(c_int), intent(inout) :: pairs(6*max_pip,2)
      integer(c_int), intent(inout) :: pair_count

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
      real(c_real) :: TMP_DTS, DTSOLID_TMP

! Numbers to calculate wall time spent in DEM calculations.
      real(c_real) :: TMP_WALL
! Pressure gradient
      real(c_real), allocatable :: gradPg(:,:,:,:)
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
         CALL OUTPUT_MANAGER(time, dt, nstep,ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g, &
             particle_state, des_radius, ro_sol, des_pos_new,&
            des_vel_new, des_usr_var, omega_new, 0, 0)
      ENDIF   ! end if/else (des_continuum_coupled)

      NP = PIP
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
            CALL DRAG_GS_DES(ep_g, u_g, v_g, w_g, ro_g, mu_g, gradPg, &
               flag, particle_state, pvol, des_pos_new, des_vel_new, fc, &
               des_radius,  particle_phase)
         ENDIF
         CALL CALC_PG_GRAD(p_g, gradPg,  particle_state, des_pos_new, &
                           pvol, drag_fc, flag, dx, dy, dz)
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

! Calculate forces from particle-wall collisions
         CALL CALC_DEM_FORCE_WITH_WALL_STL(particle_phase, particle_state,  &
            des_radius, des_pos_new, des_vel_new, omega_new, fc, tow)

! Calculate pairs of colliding particles
         CALL CALC_COLLISIONS(pairs, pair_count, particle_state, des_radius, des_pos_new)

! Calculate forces from particle-particle collisions
         CALL CALC_FORCE_DEM(particle_phase, des_radius, des_pos_new, des_vel_new, omega_new, pairs, pair_count, fc, tow)

! Calculate or distribute fluid-particle drag force.
         CALL CALC_DRAG_DES(ep_g,u_g,v_g,w_g,ro_g,mu_g,gradPg,particle_state,&
            fc,drag_fc,pvol, &
            des_pos_new,des_vel_new,des_radius,particle_phase,flag)

! Call user functions.
         IF(CALL_USR) CALL USR1_DES
! Update position and velocities
         CALL CFNEWVALUES(particle_state, des_radius, pmass, omoi, &
            des_pos_new, des_vel_new, omega_new, fc, tow, &
            des_acc_old, rot_acc_old)

! Set DO_NSEARCH before calling DES_PAR_EXCHANGE.
         DO_NSEARCH = .TRUE.

! Add/Remove particles to the system via flow BCs.
!         IF(DEM_BCMI > 0) CALL MASS_INFLOW_DEM
!         IF(DEM_BCMO > 0) CALL MASS_OUTFLOW_DEM

! Calculate mean fields (EPg).
         CALL COMP_MEAN_FIELDS(ep_g, particle_state, des_pos_new, pvol, flag, size(pvol))

         ! IF(DO_NSEARCH) CALL NEIGHBOUR(  particle_state, des_radius,&
         !    des_pos_new, neighbor_index, neighbor_index_old)


! Update time to reflect changes
         S_TIME = S_TIME + DTSOLID

! The following section targets data writes for DEM only cases:
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
! Keep track of TIME and number of steps for DEM simulations
            TIME = S_TIME
            NSTEP = NSTEP + 1
! Call the output manager to write RES data.
            CALL OUTPUT_MANAGER(time, dt, nstep,ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g, &
                particle_state, des_radius, ro_sol, &
               des_pos_new, des_vel_new, des_usr_var, omega_new, 0, 0)
         ENDIF  ! end if (.not.des_continuum_coupled)

         IF(CALL_USR) CALL USR2_DES(des_pos_new, des_vel_new, omega_new)

      ENDDO ! end do NN = 1, FACTOR

! END DEM time loop
!-----------------------------------------------------------------<<<

      IF(CALL_USR) CALL USR3_DES(des_pos_new, des_vel_new, omega_new)

! When coupled, and if needed, reset the discrete time step accordingly
      IF(DT.LT.DTSOLID_TMP) THEN
         DTSOLID = DTSOLID_TMP
      ENDIF

      IF(abs(TMP_DTS) > ZERO) THEN
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


   CONTAINS

      SUBROUTINE CALC_COLLISIONS(pairs, pair_count, particle_state, des_radius, des_pos_new)

         USE discretelement, only: nonexistent
         USE param1, only: small_number

         IMPLICIT NONE

      integer, intent(out), dimension(:,:) :: pairs
      integer, intent(out) :: pair_count

      real(c_real), intent(in) :: des_pos_new(:,:)
      real(c_real), intent(in) :: des_radius(:)
      integer, intent(in) :: particle_state(:)

      INTEGER :: i, ll
      real(c_real) :: rad
      real(c_real) :: DIST(3), DIST_MAG, POS(3)

      pair_count = 0

      DO LL = 1, PIP-1

         IF(NONEXISTENT==PARTICLE_STATE(LL)) CYCLE
         pos = DES_POS_NEW(LL,:)
         rad = DES_RADIUS(LL)

         DO I = LL+1, PIP
            IF(NONEXISTENT==PARTICLE_STATE(I)) CYCLE

            DIST(:) = DES_POS_NEW(I,:) - POS(:)
            DIST_MAG = dot_product(DIST,DIST)

            IF(DIST_MAG < (rad + DES_RADIUS(I) - SMALL_NUMBER)**2) THEN
               pair_count = pair_count + 1
               pairs(pair_count, 1) = ll
               pairs(pair_count, 2) = i
            ENDIF
         ENDDO
      ENDDO
      END SUBROUTINE CALC_COLLISIONS

      END SUBROUTINE DES_TIME_MARCH

end module des_time_march_module
