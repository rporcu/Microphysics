!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine mfix1(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
                 time, dt, u_g, v_g, w_g, &
                 p_g, ep_g, &
                 bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                 bc_klo_type, bc_khi_type, dx, dy, dz) &
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

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

      real(c_real), intent(inout) :: time, dt
      real(c_real), intent(in   ) :: dx, dy, dz

      real(c_real), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(inout) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (domlo(2)-2:domhi(2)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (domlo(1)-2:domhi(1)+2,domlo(3)-2:domhi(3)+2,2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (domlo(1)-2:domhi(1)+2,domlo(2)-2:domhi(2)+2,2)
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
      call set_ps(dx,dy,dz)

      ! Set normal velocities to zero as appropriate
      call zero_norm_vel(slo,shi,ulo,uhi,vlo,vhi,wlo,whi,u_g,v_g,w_g,&
                         bc_ilo_type, bc_ihi_type, &
                         bc_jlo_type, bc_jhi_type, &
                         bc_klo_type, bc_khi_type)

      ! Set boundary conditions
      call set_bc0(slo,shi,ulo,uhi,vlo,vhi,wlo,whi,p_g,ep_g,u_g,v_g,w_g,&
                   bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                   bc_klo_type, bc_khi_type)


      end subroutine mfix1
