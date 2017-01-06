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
      subroutine calc_coeff_all(slo, shi, lo, hi, &
        ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, &
        mu_g, f_gds, drag_bm,  particle_phase,  &
        particle_state, pvol, des_pos_new, des_vel_new, des_radius,  &
        flag, dx, dy, dz) bind(C, name="calc_coeff_all")

! Global variables:
!-----------------------------------------------------------------------
      use discretelement, only: max_pip

      ! Flag for explcit coupling between the fluid and particles.
      use discretelement, only: DES_EXPLICITLY_COUPLED

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      real(c_real), intent(inout) :: ro_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) ::  p_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ep_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: mu_g&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(OUT  ) :: f_gds&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(OUT  ) :: drag_bm&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      integer(c_int), intent(in   ) :: flag&
            (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: pvol(max_pip)
      real(c_real), intent(in   ) :: des_radius(max_pip)

      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)

      real(c_real), intent(in   ) :: dx, dy, dz

      integer(c_int), intent(in   ) :: particle_state(max_pip)
      integer(c_int), intent(in   ) :: particle_phase(max_pip)

!-----------------------------------------------------------------------
      write(6,*)'slo',slo; flush(6)
      write(6,*)'shi',shi; flush(6)

      ! Calculate all physical properties, transport properties,
      ! and exchange rates.
      CALL CALC_COEFF(slo, shi, lo, hi, flag, 2, ro_g, p_g, ep_g, rop_g, &
         u_g, v_g, w_g, mu_g, f_gds, drag_bm,  particle_phase, &
         particle_state, pvol, des_pos_new, des_vel_new, des_radius, &
         dx, dy, dz)

      if (des_explicitly_coupled) call calc_drag_des_explicit(&
         slo, shi, flag, ep_g, u_g, v_g, w_g, ro_g, mu_g, f_gds, &
         drag_bm, particle_phase,  particle_state, pvol, &
         des_pos_new, des_vel_new, des_radius, dx, dy, dz)

      end subroutine calc_coeff_all

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_COEFF                                              !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, and exchange rates.                  !
!                                                                      !
!  Author: M. Syamlal                                 Date: 25-AUG-05  !
!  Reviewer:                                          Date:            !
!                                                                      !
!                                                                      !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Variables referenced:                                               !
!  Variables modified:                                                 !
!  Local variables:                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_coeff(slo, shi, lo, hi, flag, plevel, ro_g, p_g, ep_g, &
         rop_g, u_g, v_g, w_g, mu_g, f_gds, drag_bm,  particle_phase, &
         particle_state, pvol, des_pos_new, des_vel_new, des_radius, &
         dx, dy, dz)&
        bind(C, name="calc_coeff")

      use discretelement, only: max_pip
      use discretelement, only: DES_EXPLICITLY_COUPLED
      use discretelement, only: DES_CONTinUUM_COUPLED

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
! Level to calculate physical properties.
! 0) Only density
! 1) Everything but density
! 2) All physical properties
      integer(c_int), intent(in   ) :: slo(3), shi(3), lo(3), hi(3)
      integer(c_int), intent(in   ) :: plevel

      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(inout) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) ::  p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: rop_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(out  ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(out  ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)

      real(c_real), intent(in   ) :: pvol(max_pip)
      real(c_real), intent(in   ) :: des_radius(max_pip)
      real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_real), intent(in   ) :: des_vel_new(max_pip,3)

      integer(c_int), intent(in   ) :: particle_state(max_pip)
      integer(c_int), intent(in   ) :: particle_phase(max_pip)

      real(c_real)  , intent(in   ) :: dx, dy, dz

! Calculate physical properties: (density, specific heat, diameter)
      call physical_prop(slo, shi, plevel, ro_g, p_g, ep_g, rop_g, flag)

! Calculate interphase coeffs: (momentum and energy)
      if (des_continuum_coupled .and. .not.des_explicitly_coupled)   &
         call calc_drag_des_2fluid(slo, shi, ep_g, u_g, v_g, w_g,    &
         ro_g, mu_g, f_gds, drag_bm, particle_state, particle_phase, &
         pvol, des_pos_new, des_vel_new, des_radius, dx, dy, dz)

      end subroutine calc_coeff

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_TRD_AND_TAU                                        !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, and exchange rates.                  !
!                                                                      !
!  Author: M. Syamlal                                 Date: 25-AUG-05  !
!  Reviewer:                                          Date:            !
!                                                                      !
!                                                                      !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Variables referenced:                                               !
!  Variables modified:                                                 !
!  Local variables:                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_trd_and_tau(tau_u_g,tau_v_g,tau_w_g,trd_g,&
         ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag,dx,dy,dz) &
        bind(C, name="calc_trd_and_tau")

      use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3
      use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3

      implicit none

      ! Stress tensor cross terms.
      real(c_real), intent(inout) :: tau_u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), intent(inout) :: tau_v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), intent(inout) :: tau_w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), intent(inout) :: trd_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)

      real(c_real), intent(in   ) :: ep_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), intent(in   ) :: u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), intent(in   ) :: v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), intent(in   ) :: w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), intent(in   ) :: lambda_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_real), intent(in   ) :: mu_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      integer(c_int), intent(in   ) :: flag&
            (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

      real(c_real), intent(in   ) :: dx,dy,dz

      ! Calculate the trace of the stress tensor (gas phase; m=0)
      call calc_trd_g(trd_g,u_g,v_g,w_g,flag,dx,dy,dz)

      ! Calculate the cross terms of the stress tensor (gas phase; m=0)
      call calc_tau_u_g (tau_u_g,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag,dx,dy,dz)
      call calc_tau_v_g (tau_v_g,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag,dx,dy,dz)
      call calc_tau_w_g (tau_w_g,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag,dx,dy,dz)

      end subroutine calc_trd_and_tau

end module calc_coeff_module
