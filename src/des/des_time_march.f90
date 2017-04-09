module des_time_march_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: des_time_march                                       !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine des_time_march(max_pip, slo, shi, ulo, uhi, vlo, vhi,&
         wlo, whi, lo, hi, domlo, domhi, ep_g, p_g, u_g, v_g, w_g, ro_g, mu_g, &
         particle_state, particle_phase, &
         des_radius,  pvol, pmass, omoi, des_usr_var, &
         des_pos_new, des_vel_new, omega_new, des_acc_old, rot_acc_old, &
         drag_fc, fc, tow, pairs, pair_count, &
         time, dt, dx, dy, dz, xlength, ylength, zlength, nstep) &
         bind(C, name="mfix_des_time_march")

      use calc_collision_wall, only: calc_dem_force_with_wall_stl
      use calc_drag_des_module, only: calc_drag_des
      use calc_force_dem_module, only: calc_force_dem
      use calc_pg_grad_module, only: calc_pg_grad
      use comp_mean_fields_module, only: comp_mean_fields
      use cfnewvalues_module, only: cfnewvalues
      use discretelement, only: des_continuum_coupled, des_explicitly_coupled
      use discretelement, only: dtsolid
      use discretelement, only: s_time, do_nsearch
      use drag_gs_des1_module, only: drag_gs_des
      use error_manager, only: err_msg, init_err_msg, finl_err_msg, ival, flush_err_msg
      use machine, only:  wall_time
      use output_manager_module, only: output_manager
      use param1, only: zero
      use run, only: call_usr
      use run, only: TSTOP

      implicit none

      integer(c_int), intent(in   ) :: max_pip
      integer(c_int), intent(in   ) :: slo(3), shi(3)
      integer(c_int), intent(in   ) :: ulo(3), uhi(3)
      integer(c_int), intent(in   ) :: vlo(3), vhi(3)
      integer(c_int), intent(in   ) :: wlo(3), whi(3)
      integer(c_int), intent(in   ) ::  lo(3),  hi(3)
      integer(c_int), intent(in   ) :: domlo(3), domhi(3)

      real(c_real), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real)  , intent(inout) :: time, dt
      real(c_real)  , intent(in   ) :: dx, dy, dz
      real(c_real)  , intent(in   ) :: xlength, ylength, zlength
      integer(c_int), intent(inout) :: nstep

      real(c_real), intent(inout) :: pvol(max_pip)
      real(c_real), intent(inout) :: pmass(max_pip)
      real(c_real), intent(inout) :: des_radius(max_pip)
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
! time step loop counter index
      integer :: NN

! loop counter index for any initial particle settling incoupled cases
      integer :: FACTOR

! Temporary variables when des_continuum_coupled is T to track
! changes in solid time step
      real(c_real) :: TMP_DTS, DTSOLID_TMP

! Numbers to calculate wall time spent in DEM calculations.
      real(c_real) :: TMP_WALL
! Pressure gradient
      real(c_real), allocatable :: gradPg(:,:,:,:)
!.......................................................................!

      allocate (gradPg (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3) )

! In case of restarts assign S_TIME from MFIX TIME
      S_TIME = TIME
      TMP_DTS = ZERO
      DTSOLID_TMP = ZERO
      TMP_WALL = WALL_TIME()

! Initialize time stepping variables for coupled gas/solids simulations.
      IF(des_continuum_coupled) THEN
         IF(DT.GE.DTSOLID) THEN
            FACTOR = CEILING(real(DT/DTSOLID))
         ELSE
            FACTOR = 1
            DTSOLID_TMP = DTSOLID
            DT = DTSOLID
         ENDIF

! Initialize time stepping variable for pure granular simulations.
      ELSE
         FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID))
         DT = DTSOLID
         call output_manager(max_pip, time, dt, &
            xlength, ylength, zlength, nstep, &
            particle_state, des_radius, &
            des_pos_new, des_vel_new, des_usr_var, omega_new, 0)
      ENDIF   ! end if/else (des_continuum_coupled)


      IF(des_continuum_coupled) THEN
         WRITE(ERR_MSG, 1000) trim(iVal(factor))
         call FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)
      ELSE
         WRITE(ERR_MSG, 1100) TIME, DTSOLID, trim(iVal(factor))
         call FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)
      ENDIF

 1000 FORMAT(/'DEM NITs: ',A)
 1100 FORMAT(/'Time: ',g12.5,3x,'DT: ',g12.5,3x,'DEM NITs: ',A)

      IF(call_usr) call USR0_DES

      IF(des_continuum_coupled) THEN
         IF(DES_EXPLICITLY_COUPLED) THEN
            call drag_gs_des(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, max_pip, &
               ep_g, u_g, v_g, w_g, ro_g, mu_g, &
               gradPg, particle_state, pvol, des_pos_new, &
               des_vel_new, fc, des_radius,  particle_phase, dx, dy, dz)
         ENDIF
         call calc_pg_grad(slo, shi, lo, hi, max_pip, &
                           p_g, gradPg,  particle_state, des_pos_new, &
                           pvol, drag_fc, dx, dy, dz, &
                           xlength, ylength, zlength, domlo, domhi)
      ENDIF


! Main DEM time loop
!----------------------------------------------------------------->>>
      DO NN = 1, FACTOR

         IF(des_continuum_coupled) THEN
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
         call CALC_DEM_FORCE_WITH_WALL_STL(particle_phase, particle_state,  &
            des_radius, des_pos_new, des_vel_new, omega_new, fc, tow,&
            xlength, ylength, zlength)

! Calculate pairs of colliding particles
         call CALC_COLLISIONS(max_pip, pairs, pair_count, particle_state, &
            des_radius, des_pos_new)

! Calculate forces from particle-particle collisions
         call CALC_FORCE_DEM(particle_phase, des_radius, des_pos_new, &
            des_vel_new, omega_new, pairs, pair_count, fc, tow)

! Calculate or distribute fluid-particle drag force.
         call calc_drag_des(slo,shi, ulo, uhi, vlo, vhi, wlo, whi,max_pip,&
            ep_g,u_g,v_g,w_g,ro_g,mu_g, gradPg, particle_state, fc,&
            drag_fc,pvol, des_pos_new,des_vel_new,&
            des_radius,particle_phase,dx, dy, dz)

! Call user functions.
         IF(call_usr) call USR1_DES
! Update position and velocities
         call CFNEWVALUES(max_pip, particle_state, pmass, omoi, &
            des_pos_new, des_vel_new, omega_new, fc, tow, &
            des_acc_old, rot_acc_old)

! Set DO_NSEARCH before calling DES_PAR_EXCHANGE.
         DO_NSEARCH = .TRUE.

! Add/Remove particles to the system via flow BCs.
!         IF(DEM_BCMI > 0) call MASS_INFLOW_DEM
!         IF(DEM_BCMO > 0) call MASS_OUTFLOW_DEM

! Calculate mean fields (EPg).
         call comp_mean_fields(slo, shi, max_pip, &
            ep_g, particle_state, des_pos_new, pvol, &
            dx, dy, dz)

         ! IF(DO_NSEARCH) call NEIGHBOUR(  particle_state, des_radius,&
         !    des_pos_new, neighbor_index, neighbor_index_old)


! Update time to reflect changes
         S_TIME = S_TIME + DTSOLID

! The following section targets data writes for DEM only cases:
         IF(.NOT.des_continuum_coupled) THEN
! Keep track of TIME and number of steps for DEM simulations
            TIME = S_TIME
            NSTEP = NSTEP + 1

            ! Call the output manager to write RES data.
            call output_manager(max_pip, time, dt, &
               xlength, ylength, zlength, nstep, &
               particle_state, des_radius, &
               des_pos_new, des_vel_new, des_usr_var, omega_new, 0)
         ENDIF  ! end if (.not.des_continuum_coupled)

         IF(call_usr) call USR2_DES(max_pip, des_pos_new, des_vel_new, omega_new)

      ENDDO ! end do NN = 1, FACTOR

! END DEM time loop
!-----------------------------------------------------------------<<<

      IF(call_usr) call USR3_DES(max_pip, des_pos_new, des_vel_new, omega_new)

! When coupled, and if needed, reset the discrete time step accordingly
      IF(DT.LT.DTSOLID_TMP) THEN
         DTSOLID = DTSOLID_TMP
      ENDIF

      IF(abs(TMP_DTS) > ZERO) THEN
         DTSOLID = TMP_DTS
         TMP_DTS = ZERO
      ENDIF

      deallocate(gradPg)

      IF(.NOT.des_continuum_coupled)THEN
         WRITE(ERR_MSG,"('<---------- END DES_TIME_MARCH ----------')")
         call FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ELSE
         ! call send_recv(ep_g,2)
         ! call send_recv(rop_g,2)

         TMP_WALL = WALL_TIME() - TMP_WALL
         IF(TMP_WALL > 1.0d-10) THEN
            WRITE(ERR_MSG, 9000) trim(iVal(dble(FACTOR)/TMP_WALL))
         ELSE
            WRITE(ERR_MSG, 9000) '+Inf'
         ENDIF
         call FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)

 9000 FORMAT('    NITs/SEC = ',A)

      ENDIF


   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: calc_collisions                                          !
! Purpose: Build a list of collision pairs.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_COLLISIONS(max_pip, pairs, pair_count, &
      particle_state, des_radius, des_pos_new)

      use discretelement, only: nonexistent
      use param1, only: small_number

      implicit none

! Dummy arguments ......................................................
      integer, intent(in) :: max_pip

      integer     , intent(  out) :: pair_count
      integer     , intent(  out) :: pairs(6*max_pip,2)

      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_radius(max_pip)
      integer     , intent(in   ) :: particle_state(max_pip)

! Local variables ......................................................
      integer :: i, ll
      real(c_real) :: rad
      real(c_real) :: DIST(3), DIST_MAG, POS(3)

      pair_count = 0

      DO LL = 1, max_PIP-1

         IF(NONEXISTENT==PARTICLE_STATE(LL)) CYCLE
         pos = DES_POS_NEW(LL,:)
         rad = DES_RADIUS(LL)

         DO I = LL+1, max_PIP
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
      end subroutine calc_collisions

      end subroutine des_time_march

end module des_time_march_module
