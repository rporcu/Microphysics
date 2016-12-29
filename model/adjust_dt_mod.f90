module adjust_dt

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ADJUST_DT(IER, NIT)                                    !
!  Author: M. Syamlal                                 Date: FEB-10-97  !
!                                                                      !
!  Purpose: Automatically adjust time step.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   integer function adjustdt (ier, nit, dt) &
      bind(C, name="mfix_adjustdt")

      use iso_c_binding, only: c_double, c_int

! Global Variables:
!---------------------------------------------------------------------//
! User defined aximum number of iterations
      use leqsol, only: MAX_NIT
! User defined: min, max DT and adjustment factor
      use run, only: DT_MIN, DT_MAX, DT_FAC
! Current DT (1/DT) and direction of last change (+/-)
      use run, only: DT_DIR

      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO, ONE, UNDEFINED

! Module proceedures:
!---------------------------------------------------------------------//
! Routine to break successive time step reductions.
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg
      use error_manager, only: init_err_msg, ival

      use calc_coeff_module, only: calc_coeff_all

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! integer flag: 0=good, 100=initialize, otherwise bad.
      integer(c_int), intent(inout) :: ier
! number of iterations for current time step
      integer(c_int), intent(in) :: nit
! fluid solver time-step size
      real(c_double), intent(inout) :: dt

! Local Variables:
!---------------------------------------------------------------------//
! Number of steps in between DT adjustments.
      INTEGER, PARAMETER :: STEPS_MIN = 5
! Number of time steps since last DT adjustment
      INTEGER, SAVE :: STEPS_TOT=0
! number of iterations since last DT adjustment
      INTEGER, SAVE :: NIT_TOT=0
! Iterations per second for last dt
      DOUBLE PRECISION, SAVE :: NIToS=0.0
! Current number of iterations per second
      DOUBLE PRECISION :: NITOS_NEW
!......................................................................!

! Initialize the function result.
      adjustdt = 0

! Steady-state simulation.
      if (dt==undefined .or. dt<=zero) return

! Iterate successfully converged.
!---------------------------------------------------------------------//
      if(ier == 0) then

! Calculate a new DT every STEPS_MIN time steps.
         if(steps_tot >= steps_min) then
            nitos_new = dble(nit_tot)/(steps_tot*dt)
            if (nitos_new > nitos) dt_dir = dt_dir*(-1)
            steps_tot = 0
            nitos = nitos_new
            nit_tot = 0
            if (dt_dir > 0) then
               if(nit < max_nit) dt = min(dt_max,dt/dt_fac)
            else
               dt = dt*dt_fac
            endif

! Write the convergence stats to the screen/log file.
            WRITE(ERR_MSG,"('DT=',g11.4,3x,'NIT/s=',A)")  &
               DT, trim(iVal(nint(NITOS)))
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., &
               FOOTER=.FALSE., LOG=.FALSE.)

         else
            steps_tot = steps_tot + 1
            nit_tot = nit_tot + nit
         endif
! No need to iterate again
         adjustdt = 0

! Iterate failed to converge.
!---------------------------------------------------------------------//
      else

! Clear the error flag.
         IER = 0

! Reset counters.
         steps_tot = 0
         nitos = 0.
         nit_tot = 0

! Reduce the step size.
         dt = dt*dt_fac

! The step size has decreased to the minimum.
         if (dt_fac >= one) then

            write(err_msg,"(3x,a)") &
               'DT_FAC >= 1. Recovery not possible!'
            call flush_err_msg(abort=.true., &
               header=.false., footer=.false.)

         elseif (dt > dt_min) then

            WRITE(ERR_MSG,"(3X,'Recovered: Dt=',G12.5,' :-)')") DT
            call flush_err_msg(header=.false., footer=.false.)

!             ep_g = ep_go
!             p_g =  p_go
!             ro_g = ro_go
!             rop_g = rop_go
!             U_g =  U_go
!             V_g =  V_go
!             W_g =  W_go

! ! Recalculate all coefficients
!             CALL CALC_COEFF_ALL(ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, &
!                mu_g, f_gds, drag_bm,  particle_phase,  &
!                particle_state, pvol, des_pos_new, des_vel_new, des_radius,  &
!                flag)
! Iterate again with new dt
            adjustdt = 1

! Set the return flag stop iterating.
         ELSE

! Prevent DT from dropping below DT_MIN.
            adjustdt = 0
         ENDIF

      ENDIF

      RETURN
      END FUNCTION ADJUSTDT

      END MODULE ADJUST_DT
