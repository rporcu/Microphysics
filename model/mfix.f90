!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine mfix1(time, dt, u_g, v_g, w_g, &
   p_g, ep_g, ro_g, rop_g, &
   d_e, d_n, d_t, &
   flux_ge, flux_gn, flux_gt, &
   trD_g, lambda_g, mu_g, flag) &
   bind(C, name="mfix_main1")

!-----------------------------------------------
! Modules
!-----------------------------------------------

      use calc_coeff_module, only: calc_coeff
      use compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use corner_module, only: get_corner_cells
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg
      use exit_mod, only: mfix_exit
      use geometry, only: dx, dy, dz, ayz, axy, axz, vol, flag_mod
      use iso_c_binding, only: c_double, c_int
      use machine, only: wall_time
      use output_manager_module, only: init_output_vars
      use param1 , only: is_defined, is_undefined, undefined
      use parse_resid_string_module, only: parse_resid_string
      use run, only: run_type
      use set_bc0_module, only: set_bc0
      use set_constprop_module, only: set_constprop
      use set_flags_module, only: set_flags1
      use set_ps_module, only: set_ps
      use write_out0_module, only: write_out0
      use write_out1_module, only: write_out1
      use write_out3_module, only: write_out3
      use zero_norm_vel_module, only: zero_norm_vel
      use init_fluid_module, only: init_fluid

      use fld_const, only: ro_g0

      IMPLICIT NONE

      real(c_double), intent(inout) :: time, dt

      integer(c_int), intent(inout) :: flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

      real(c_double), intent(inout) :: ep_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: p_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: ro_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: rop_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      real(c_double), intent(inout) :: u_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: v_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: w_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      real(c_double), intent(inout) :: d_e&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: d_t&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: d_n&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      real(c_double), intent(inout) :: mu_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: lambda_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: trD_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      real(c_double), intent(inout) :: flux_gE&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: flux_gN&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: flux_gT&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

!---------------------------------------------------------------------//
      flag_mod = flag

      CALL INIT_OUTPUT_VARS(time, dt)

! Parse residual strings
      CALL PARSE_RESID_STRING ()

! Write the initial part of the standard output file
      CALL WRITE_OUT0(time, dt)

! Write the initial part of the special output file(s)
      CALL WRITE_USR0

      CALL INIT_ERR_MSG('MFIX')

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

! Set point sources.
      CALL SET_PS(flag)


! Initialize fluid variables .........................................

! Set boundary conditions
      CALL ZERO_NORM_VEL(u_g,v_g,w_g, flag)
      CALL SET_BC0(p_g,ep_g,u_g,v_g,w_g,ro_g0,flag)

      call init_fluid(ep_g, ro_g, rop_g, p_g, u_g, v_g, w_g, &
         mu_g, lambda_g, flag)

!.....................................................................


      end subroutine mfix1
