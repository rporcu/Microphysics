!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine mfix1(slo, shi, lo, hi, time, dt, u_g, v_g, w_g, &
   p_g, ep_g, ro_g, rop_g, &
   d_e, d_n, d_t, &
   flux_ge, flux_gn, flux_gt, &
   trD_g, lambda_g, mu_g, flag, &
   dx, dy, dz) &
   bind(C, name="mfix_main1")

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use calc_coeff_module, only: calc_coeff
      use corner_module, only: get_corner_cells
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg
      use exit_mod, only: mfix_exit
      use geometry, only: flag_mod
      use machine, only: wall_time
      use output_manager_module, only: init_output_vars
      use param1 , only: is_defined, is_undefined
      use parse_resid_string_module, only: parse_resid_string
      use set_bc0_module, only: set_bc0
      use set_flags_module, only: set_flags1
      use set_ps_module, only: set_ps
      use write_out0_module, only: write_out0
      use write_out1_module, only: write_out1
      use write_out3_module, only: write_out3
      use zero_norm_vel_module, only: zero_norm_vel

      use fld_const, only: ro_g0

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      real(c_real), intent(inout) :: time, dt
      real(c_real), intent(in   ) :: dx, dy, dz

      integer(c_int), intent(inout) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(inout) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(inout) :: d_e&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: d_t&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: d_n&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(inout) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: lambda_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: trD_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(inout) :: flux_gE&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_gN&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: flux_gT&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

!---------------------------------------------------------------------//
      flag_mod = flag

      call init_output_vars(time, dt)

! Parse residual strings
      call parse_resid_string ()

! Write the initial part of the standard output file
      call write_out0(time, dt, dx, dy, dz)

! Write the initial part of the special output file(s)
      call write_usr0

      call init_err_msg('MFIX')

! Set the flags for wall surfaces impermeable and identify flow
! boundaries using FLAG_E, FLAG_N, and FLAG_T
      call set_flags1(slo,shi,flag)
      flag_mod = flag

      ! Find corner cells and set their face areas to zero
      call get_corner_cells(slo,shi,lo,hi,flag)
      flag_mod = flag

! Set point sources.
      call set_ps(slo,shi,lo,hi,flag,dx,dy,dz)

      ! Set normal velocities to zero as appropriate
      call zero_norm_vel(slo,shi,lo,hi,u_g,v_g,w_g,flag)

      ! Set boundary conditions
      call set_bc0(slo,shi,lo,hi,p_g,ep_g,u_g,v_g,w_g,ro_g0,flag)

      end subroutine mfix1
