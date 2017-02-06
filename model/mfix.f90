!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine mfix1(slo, shi, lo, hi, time, dt, u_g, v_g, w_g, &
                 p_g, ep_g, &
                 bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                 bc_klo_type, bc_khi_type, flag, dx, dy, dz) &
   bind(C, name="mfix_main1")

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use calc_coeff_module, only: calc_coeff
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg
      use exit_mod, only: mfix_exit
      use geometry, only: domlo, domhi
      use machine, only: wall_time
      use output_manager_module, only: init_output_vars
      use param1 , only: is_defined, is_undefined
      use parse_resid_string_module, only: parse_resid_string
      use set_bc0_module, only: set_bc0
      use set_ps_module, only: set_ps
      use write_out0_module, only: write_out0
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

      real(c_real), intent(inout) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (slo(1):shi(1),slo(2):shi(2),2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (slo(1):shi(1),slo(2):shi(2),2)
!---------------------------------------------------------------------//
      call init_output_vars(time, dt)

      ! Parse residual strings
      call parse_resid_string ()

      ! Write the initial part of the standard output file
      call write_out0(time, dt, dx, dy, dz)

      ! Write the initial part of the special output file(s)
      call write_usr0

      call init_err_msg('MFIX')

      ! Set point sources.
      call set_ps(slo,shi,flag,dx,dy,dz)

      ! Set normal velocities to zero as appropriate
      call zero_norm_vel(slo,shi,u_g,v_g,w_g,flag)

      ! Set boundary conditions
      call set_bc0(slo,shi,p_g,ep_g,u_g,v_g,w_g,ro_g0, &
                   bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                   bc_klo_type, bc_khi_type, flag)


      end subroutine mfix1
