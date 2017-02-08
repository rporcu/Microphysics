module calc_coeff_module

      use calc_drag_des_module, only: calc_drag_des_explicit, calc_drag_des_2fluid
      use calc_tau_u_g_module, only: calc_tau_u_g
      use calc_tau_v_g_module, only: calc_tau_v_g
      use calc_tau_w_g_module, only: calc_tau_w_g
      use calc_trd_g_module, only: calc_trd_g
      use physical_prop_module, only: physical_prop

      use bl_fort_module, only: c_real
      use iso_c_binding , only: c_int

  contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      !
!  Subroutine: CALC_COEFF_ALL                                          !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, exchange rates, and reaction rates.  !
!                                                                      !
!  Author: M. Syamlal                                 Date: 25-AUG-05  !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Variables referenced:                                               !
!  Variables modified:                                                 !
!  Local variables:                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_coeff_all(slo, shi, &
        ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
        max_pip, &
        ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, &
        mu_g, f_gds, drag_bm,  particle_phase,  &
        particle_state, pvol, des_pos_new, des_vel_new, des_radius,  &
        flag, dx, dy, dz) bind(C, name="calc_coeff_all")

! Global variables:
!-----------------------------------------------------------------------

      ! Flag for explcit coupling between the fluid and particles.
      use discretelement, only: DES_EXPLICITLY_COUPLED

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: max_pip

      real(c_real), intent(in   ) ::  p_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_g&
            (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
            (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: mu_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(inout) :: ro_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ep_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: f_gds&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: drag_bm&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      real(c_real), intent(in   ) :: pvol(max_pip)
      real(c_real), intent(in   ) :: des_radius(max_pip)

      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)

      real(c_real), intent(in   ) :: dx, dy, dz

      integer(c_int), intent(in   ) :: particle_state(max_pip)
      integer(c_int), intent(in   ) :: particle_phase(max_pip)

!-----------------------------------------------------------------------

      ! Calculate all physical properties, transport properties,
      ! and exchange rates.
      call calc_coeff(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
         max_pip, 2, ro_g, p_g, ep_g, &
         rop_g, u_g, v_g, w_g, mu_g, f_gds, drag_bm,  particle_phase, &
         particle_state, pvol, des_pos_new, des_vel_new, des_radius, &
         dx, dy, dz)

      if (des_explicitly_coupled) &
         call calc_drag_des_explicit(&
            slo, shi, max_pip, flag, ep_g, u_g, v_g, w_g, ro_g, mu_g, f_gds, &
            drag_bm, particle_phase,  particle_state, pvol, &
            des_pos_new, des_vel_new, des_radius, dx, dy, dz)

      end subroutine calc_coeff_all

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_COEFF                                              !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, and exchange rates.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_coeff(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, lo, hi, &
         max_pip, plevel, &
         ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, mu_g, f_gds, drag_bm,&
         particle_phase, particle_state, pvol, des_pos_new, des_vel_new,&
         des_radius, dx, dy, dz)&
        bind(C, name="calc_coeff")

      use discretelement, only: des_explicitly_coupled
      use discretelement, only: des_continuum_coupled

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
! Level to calculate physical properties.
! 0) Only density
! 1) Everything but density
! 2) All physical properties
      integer(c_int), intent(in   ) :: slo(3), shi(3), lo(3), hi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: max_pip

      integer(c_int), intent(in   ) :: plevel

      real(c_real), intent(in   ) ::  p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_g&
            (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
            (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(inout) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(  out) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      real(c_real), intent(in   ) :: pvol(max_pip)
      real(c_real), intent(in   ) :: des_radius(max_pip)
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)

      integer(c_int), intent(in   ) :: particle_state(max_pip)
      integer(c_int), intent(in   ) :: particle_phase(max_pip)

      real(c_real)  , intent(in   ) :: dx, dy, dz

! Calculate physical properties: (density, specific heat, diameter)
      call physical_prop(slo, shi, lo, hi, plevel, ro_g, p_g, ep_g, rop_g)

! Calculate interphase coeffs: (momentum and energy)
      if (des_continuum_coupled .and. .not.des_explicitly_coupled)   &
         call calc_drag_des_2fluid(slo, shi, max_pip, ep_g,  &
         u_g, v_g, w_g, ro_g, mu_g, f_gds, drag_bm, particle_state,  &
         particle_phase, pvol, des_pos_new, des_vel_new, des_radius, &
         dx, dy, dz)

      end subroutine calc_coeff

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_TRD_AND_TAU                                        !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, and exchange rates.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_trd_and_tau(slo,shi,ulo,uhi,vlo,vhi,wlo,whi,lo,hi,&
         tau_u_g,tau_v_g,tau_w_g,trd_g,&
         u_g,v_g,w_g,lambda_g,mu_g,dx,dy,dz) &
        bind(C, name="calc_trd_and_tau")

      implicit none

      integer(c_int), intent(in ) :: slo(3),shi(3)
      integer(c_int), intent(in ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in ) ::  lo(3), hi(3)

      ! Stress tensor cross terms.
      real(c_real), intent(inout) :: tau_u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: tau_v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: tau_w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: trd_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: u_g&
            (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: v_g&
            (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(in   ) :: lambda_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), intent(in   ) :: dx,dy,dz

      ! Calculate the trace of the stress tensor (gas phase; m=0)
      call calc_trd_g(slo,shi,lo,hi,trd_g,u_g,v_g,w_g,dx,dy,dz)

      ! Calculate the cross terms of the stress tensor (gas phase; m=0)
      call calc_tau_u_g (slo,shi,lo,hi,tau_u_g,trd_g,u_g,v_g,w_g,lambda_g,mu_g,dx,dy,dz)
      call calc_tau_v_g (slo,shi,lo,hi,tau_v_g,trd_g,u_g,v_g,w_g,lambda_g,mu_g,dx,dy,dz)
      call calc_tau_w_g (slo,shi,lo,hi,tau_w_g,trd_g,u_g,v_g,w_g,lambda_g,mu_g,dx,dy,dz)

      end subroutine calc_trd_and_tau

end module calc_coeff_module
