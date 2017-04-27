module des_time_march_module

  use amrex_fort_module, only: c_real => amrex_real
  use iso_c_binding ,    only: c_int

  implicit none

  private
  public des_time_march

  ! Let's make some variables module variables to make the transition
  ! to C++ easier ( even though it is really a nasty thing to do ).
  integer(c_int), save  :: factor      
  real(c_real),   save  :: dtsolid_tmp 

  ! Temporary variables when des_continuum_coupled is T to track
  ! changes in solid time step
  real(c_real),   save  :: TMP_DTS

  ! Pressure gradient
  real(c_real),   allocatable, save :: gradPg(:,:,:,:), fc(:,:), tow(:,:)

  ! Define interface for external functions so there will be no 
  ! warning at compile time
  interface

     subroutine usr0_des()
     end subroutine usr0_des

     subroutine usr1_des()
     end subroutine usr1_des

     subroutine usr2_des(np, pos, vel, omega)
       import  c_real
       integer     , intent(in)    :: np
       real(c_real), intent(in)    :: pos(np,3), vel(np,3)
       real(c_real), intent(inout) :: omega(np,3)
     end subroutine usr2_des

     subroutine usr3_des(np, pos, vel, omega)
       import  c_real
       integer     , intent(in)    :: np
       real(c_real), intent(in)    :: pos(np,3), vel(np,3), omega(np,3)
     end subroutine usr3_des

  end interface


contains



  subroutine des_init_time_loop( np, slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
       & lo, hi, domlo, domhi, ep_g, p_g, u_g, v_g, w_g, ro_g, mu_g,     &
       & state, phase, radius, vol, pos, vel, omega, drag,&
       & time, dt, dx, dy, dz, xlength, ylength, zlength, nstep) &
       bind(C, name="mfix_des_init_time_loop")

    use discretelement,          only: dtsolid, s_time
    use calc_drag_des_module,    only: calc_drag_des
    use calc_pg_grad_module,     only: calc_pg_grad
    use drag_gs_des1_module,     only: drag_gs_des
    use error_manager,           only: err_msg, init_err_msg, finl_err_msg, ival, flush_err_msg
    use output_manager_module,   only: output_manager
    use discretelement,          only: des_continuum_coupled, des_explicitly_coupled
    use run,                     only: call_usr, tstop
    use param,                   only: zero

    integer(c_int), intent(in)    :: np, slo(3), shi(3), ulo(3), uhi(3)
    integer(c_int), intent(in)    :: vlo(3), vhi(3), wlo(3), whi(3)
    integer(c_int), intent(in)    :: lo(3), hi(3), domlo(3), domhi(3)
    real(c_real),   intent(in)    :: xlength, ylength, zlength, dx, dy, dz
    real(c_real),   intent(inout) :: time, dt
    real(c_real),   intent(in)    :: &
         p_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
         u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
         v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
         w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
         ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
         mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(c_real), intent(inout) :: vol(np), radius(np)
    real(c_real), intent(inout) :: drag(np,3), omega(np,3)
    real(c_real), intent(inout) :: pos(np,3), vel(np,3)

    real(c_real), intent(inout) :: &
         ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    integer(c_int), intent(inout) :: nstep, state(np), phase(np)


    allocate (gradPg (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3) )
    allocate ( tow(np,3), fc(np,3) )

    tow        = zero
    fc         = zero
    
    
    ! In case of restarts assign S_TIME from MFIX TIME
    S_TIME = TIME
    TMP_DTS = ZERO    

    ! Initialize time stepping variables for coupled gas/solids simulations.
    if ( des_continuum_coupled ) then

       if ( DT >= DTSOLID ) then
          FACTOR = ceiling(real(DT/DTSOLID))
       else
          FACTOR = 1
          DTSOLID_TMP = DTSOLID
          DT = DTSOLID
       end if
       
       ! Initialize time stepping variable for pure granular simulations.
    else

       FACTOR = ceiling(real((TSTOP-TIME)/DTSOLID))
       DT = DTSOLID
       call output_manager(np, time, dt, xlength, ylength, zlength, nstep, &
            state, radius, pos, vel, omega, 0)

    end if   

    if (des_continuum_coupled) then
       write(ERR_MSG, 1000) trim(iVal(factor))
       call FLUSH_ERR_MSG(HEADER=.false., FOOTER=.false., LOG=.false.)
    else
       write(ERR_MSG, 1100) TIME, DTSOLID, trim(iVal(factor))
       call FLUSH_ERR_MSG(HEADER=.false., FOOTER=.false., LOG=.false.)
    endif

1000 FORMAT(/'DEM NITs: ',A)
1100 FORMAT(/'Time: ',g12.5,3x,'DT: ',g12.5,3x,'DEM NITs: ',A)

    IF(call_usr) call USR0_DES

    IF(des_continuum_coupled) THEN
       IF(DES_EXPLICITLY_COUPLED) THEN
          call drag_gs_des(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, np, &
               ep_g, u_g, v_g, w_g, ro_g, mu_g, &
               gradPg, state, vol, pos, vel, fc, radius, phase,&
               dx, dy, dz)
       ENDIF
       call calc_pg_grad(slo, shi, lo, hi, np, &
            p_g, gradPg,  state, pos, vol, drag, dx, dy, dz, &
            xlength, ylength, zlength, domlo, domhi)
    ENDIF

    

  end subroutine des_init_time_loop


  subroutine des_finalize_time_loop( np, dt, pos, vel, omega ) &
       bind(C, name="mfix_des_finalize_time_loop")

    use discretelement, only: dtsolid
    use param,          only: zero
    use run,            only: call_usr

    integer(c_int), intent(in)    :: np
    real(c_real),   intent(inout) :: dt, pos(np,3), vel(np,3), omega(np,3)    

    if ( call_usr ) call usr3_des(np, pos, vel, omega)

    ! When coupled, and if needed, reset the discrete time step accordingly
    if ( DT < DTSOLID_TMP ) then
       DTSOLID = DTSOLID_TMP
    endif

    if ( abs(TMP_DTS) > ZERO ) then
       DTSOLID = TMP_DTS
       TMP_DTS = ZERO
    endif

    deallocate( gradPg, tow, fc )

  end subroutine des_finalize_time_loop

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !     Subroutine: des_time_march                                       !
  !                                                                      !
  !     Purpose: Main DEM driver routine                                 !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine des_time_march( np, slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
       & lo, hi, domlo, domhi, ep_g, p_g, u_g, v_g, w_g, ro_g, mu_g,     &
       & state, phase, radius, vol, mass, omoi, pos, vel, omega, acc, alpha, drag,&
       & time, dt, dx, dy, dz, xlength, ylength, zlength, nstep) &
       bind(C, name="mfix_des_time_march")

    use calc_collision_wall,     only: calc_dem_force_with_wall_stl
    use calc_drag_des_module,    only: calc_drag_des
    use calc_force_dem_module,   only: calc_force_dem
    use calc_pg_grad_module,     only: calc_pg_grad
    use comp_mean_fields_module, only: comp_mean_fields
    use cfnewvalues_module,      only: cfnewvalues
    use discretelement,          only: dtsolid, s_time, des_continuum_coupled
    use drag_gs_des1_module,     only: drag_gs_des
    use error_manager,           only: init_err_msg, finl_err_msg, ival, flush_err_msg
    use output_manager_module,   only: output_manager

    integer(c_int), intent(in)    :: np, slo(3), shi(3), ulo(3), uhi(3)
    integer(c_int), intent(in)    :: vlo(3), vhi(3), wlo(3), whi(3)
    integer(c_int), intent(in)    :: lo(3), hi(3), domlo(3), domhi(3)
    real(c_real),   intent(in)    :: xlength, ylength, zlength, dx, dy, dz
    real(c_real),   intent(inout) :: time, dt
    real(c_real),   intent(in)    :: &
         p_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
         u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
         v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
         w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
         ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
         mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(c_real), intent(inout) :: vol(np), mass(np), radius(np), omoi(np)
    real(c_real), intent(inout) :: alpha(np,3), drag(np,3), omega(np,3)
    real(c_real), intent(inout) :: pos(np,3), vel(np,3), acc(np,3)

    real(c_real), intent(inout) :: &
         ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    integer(c_int), intent(inout) :: nstep, state(np), phase(np)

    !------------------------------------------------
    ! Local variables
    !------------------------------------------------
    ! time step loop counter index
    integer :: NN

    call des_init_time_loop( np, slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
       & lo, hi, domlo, domhi, ep_g, p_g, u_g, v_g, w_g, ro_g, mu_g,     &
       & state, phase, radius, vol, pos, vel, omega, drag,&
       & time, dt, dx, dy, dz, xlength, ylength, zlength, nstep) 

    ! Main DEM time loop
    !----------------------------------------------------------------->>>
    do NN = 1, FACTOR

       if ( des_continuum_coupled ) then
          ! If the current time in the discrete loop exceeds the current time in
          ! the continuum simulation, exit the discrete loop
          if ( S_TIME > (TIME+DT) ) exit
          ! If next time step in the discrete loop will exceed the current time
          ! in the continuum simulation, modify the discrete time step so final
          ! time will match
          if ( (S_TIME+DTSOLID) > (TIME+DT) ) then
             TMP_DTS = DTSOLID
             DTSOLID = TIME + DT - S_TIME
          end if
       end if

       call des_time_loop_ops( np, slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
            & ep_g, u_g, v_g, w_g, ro_g, mu_g,     &
            & state, phase, radius, vol, mass, omoi, pos, vel, omega, acc, alpha, drag,&
            & time, dt, dx, dy, dz, xlength, ylength, zlength, nstep) 

    ENDDO 

    call des_finalize_time_loop( np, dt, pos, vel, omega )

  end subroutine des_time_march

  !
  ! This subroutine performs a single time step
  ! 
  subroutine des_time_loop_ops( np, slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
       & ep_g, u_g, v_g, w_g, ro_g, mu_g,     &
       & state, phase, radius, vol, mass, omoi, pos, vel, omega, acc, alpha, drag,&
       & time, dt, dx, dy, dz, xlength, ylength, zlength, nstep)  &
       bind(C, name="mfix_des_time_loop_ops")

    use calc_collision_wall,     only: calc_dem_force_with_wall_stl
    use calc_drag_des_module,    only: calc_drag_des
    use calc_force_dem_module,   only: calc_force_dem
    use calc_pg_grad_module,     only: calc_pg_grad
    use comp_mean_fields_module, only: comp_mean_fields
    use cfnewvalues_module,      only: cfnewvalues
    use discretelement,          only: dtsolid, s_time, des_continuum_coupled
    use drag_gs_des1_module,     only: drag_gs_des
    use error_manager,           only: init_err_msg, finl_err_msg, ival, flush_err_msg
    use output_manager_module,   only: output_manager
    use run,                     only: call_usr

    integer(c_int), intent(in)    :: np, slo(3), shi(3), ulo(3), uhi(3)
    integer(c_int), intent(in)    :: vlo(3), vhi(3), wlo(3), whi(3)
    real(c_real),   intent(in)    :: xlength, ylength, zlength, dx, dy, dz
    real(c_real),   intent(inout) :: time, dt
    real(c_real),   intent(in)    :: &
         u_g(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3)), &
         v_g(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3)), &
         w_g(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3)), &
         ro_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)), &
         mu_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(c_real), intent(inout) :: vol(np), mass(np), radius(np), omoi(np)
    real(c_real), intent(inout) :: alpha(np,3), drag(np,3), omega(np,3)
    real(c_real), intent(inout) :: pos(np,3), vel(np,3), acc(np,3)

    real(c_real), intent(inout) :: &
         ep_g(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    integer(c_int), intent(inout) :: nstep, state(np), phase(np)


    ! calculate forces from particle-wall collisions
    call calc_dem_force_with_wall_stl( phase, state, radius, &
         pos, vel, omega, fc, tow, xlength, ylength, zlength )

    ! calculate forces from particle-particle collisions
    call calc_force_dem(phase, radius, pos, vel, omega, state, & 
         fc, tow)

    ! calculate or distribute fluid-particle drag force.
    call calc_drag_des(slo,shi, ulo, uhi, vlo, vhi, wlo, whi,np,&
         ep_g,u_g,v_g,w_g,ro_g,mu_g, gradpg, state, fc,&
         drag, vol, pos, vel,&
         radius, phase,dx, dy, dz)

    ! call user functions.
    if ( call_usr ) call usr1_des

    ! update position and velocities
    call cfnewvalues(np, state, mass, omoi, &
         pos, vel, omega, fc, tow, &
         acc, alpha)


    ! calculate mean fields (epg).
    call comp_mean_fields(slo, shi, np, &
         ep_g, state, pos, vol, &
         dx, dy, dz)


    ! update time to reflect changes
    s_time = s_time + dtsolid

    ! the following section targets data writes for dem only cases:
    if(.not.des_continuum_coupled) then
       ! keep track of time and number of steps for dem simulations
       time = s_time
       nstep = nstep + 1

       ! call the output manager to write res data.
       call output_manager(np, time, dt, &
            xlength, ylength, zlength, nstep, &
            state, radius, &
            pos, vel, omega, 0)
    endif  ! end if (.not.des_continuum_coupled)

    if(call_usr) call usr2_des(np, pos, vel, omega)

  end subroutine des_time_loop_ops



end module des_time_march_module
