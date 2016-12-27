module calc_coeff_module

      use calc_drag_des_module, only: calc_drag_des_explicit, calc_drag_des_2fluid
      use calc_tau_u_g_module, only: calc_tau_u_g
      use calc_tau_v_g_module, only: calc_tau_v_g
      use calc_tau_w_g_module, only: calc_tau_w_g
      use calc_trd_g_module, only: calc_trd_g
      use physical_prop_module, only: physical_prop
      use iso_c_binding, only: c_double, c_int

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
     SUBROUTINE CALC_COEFF_ALL(ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, &
        mu_g, f_gds, drag_bm,  particle_phase,  &
        particle_state, pvol, des_pos_new, des_vel_new, des_radius,  &
        flag) bind(C, name="mfix_calc_coeff_all")

! Global variables:
!-----------------------------------------------------------------------
      use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3
      use discretelement, only: max_pip

      ! Flag for explcit coupling between the fluid and particles.
      use discretelement, only: DES_EXPLICITLY_COUPLED

      implicit none

      real(c_double), INTENT(INOUT) :: ro_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(IN   ) ::  p_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(INOUT) :: ep_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(INOUT) :: rop_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(IN   ) :: u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(IN   ) :: v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(IN   ) :: w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(INOUT) :: mu_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(OUT  ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), INTENT(OUT  ) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)


      integer(c_int), intent(in   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)

      real(c_double), intent(in   ) :: pvol(max_pip)
      real(c_double), intent(in   ) :: des_radius(max_pip)

      real(c_double), intent(in   ) :: des_pos_new(max_pip,3)
      real(c_double), intent(in   ) :: des_vel_new(max_pip,3)

      integer(c_int)         , intent(in   ) :: particle_state(max_pip)
      integer(c_int)         , intent(in   ) :: particle_phase(max_pip)


!-----------------------------------------------------------------------

      ! Calculate all physical properties, transport properties,
      ! and exchange rates.
      CALL CALC_COEFF(flag, 2, ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, mu_g, &
         f_gds, drag_bm,  particle_phase, particle_state, pvol, &
         des_pos_new, des_vel_new, des_radius)

      IF (DES_EXPLICITLY_COUPLED) CALL CALC_DRAG_DES_EXPLICIT(flag, ep_g, &
         u_g, v_g, w_g, ro_g, mu_g, f_gds, drag_bm,  &
         particle_phase,  particle_state, &
         pvol, des_pos_new, des_vel_new, des_radius)

      END SUBROUTINE CALC_COEFF_ALL

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
      SUBROUTINE CALC_COEFF(flag, pLevel, ro_g, p_g, ep_g, rop_g, u_g, v_g, &
         w_g, mu_g, f_gds, drag_bm,  particle_phase, particle_state, &
         pvol, des_pos_new, des_vel_new, des_radius)

      use fld_const, only: ro_g0
      use compar   , only: istart3,iend3,jstart3,jend3,kstart3,kend3
      use discretelement, only: max_pip
      use discretelement, only: DES_EXPLICITLY_COUPLED
      use discretelement, only: DES_CONTINUUM_COUPLED

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
! Level to calculate physical properties.
! 0) Only density
! 1) Everything but density
! 2) All physical properties
      INTEGER, intent(in) :: pLevel
      double precision, intent(inout) :: ro_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) ::  p_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) :: ep_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) :: u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) :: v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) :: w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: mu_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(out  ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(out  ) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
      integer         , intent(in   ) :: flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

      double precision, intent(in   ) :: pvol(max_pip)
      double precision, intent(in   ) :: des_radius(max_pip)
      double precision, intent(in   ) :: des_pos_new(max_pip,3)
      double precision, intent(in   ) :: des_vel_new(max_pip,3)

      integer         , intent(in   ) :: particle_state(max_pip)
      integer         , intent(in   ) :: particle_phase(max_pip)


!-----------------------------------------------------------------------
      integer IER

! Calculate physical properties: (density, specific heat, diameter)
      CALL PHYSICAL_PROP(IER, pLevel, ro_g, p_g, ep_g, rop_g, ro_g0, flag)

! Calculate interphase coeffs: (momentum and energy)
      IF (DES_CONTINUUM_COUPLED .AND. .NOT.DES_EXPLICITLY_COUPLED)  &
         CALL CALC_DRAG_DES_2FLUID(ep_g, u_g, v_g, w_g, ro_g, mu_g, &
         f_gds, drag_bm,  particle_state, &
         particle_phase, pvol, des_pos_new, des_vel_new, des_radius)

      END SUBROUTINE CALC_COEFF

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
      SUBROUTINE CALC_TRD_AND_TAU(tau_u_g,tau_v_g,tau_w_g,trd_g,&
         ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag) &
        bind(C, name="mfix_calc_trd_and_tau")

      use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3

      implicit none

      ! Stress tensor cross terms.
      real(c_double), INTENT(INOUT) :: tau_u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(INOUT) :: tau_v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(INOUT) :: tau_w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(INOUT) :: trd_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)

      real(c_double), INTENT(IN   ) :: ep_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(IN   ) :: u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(IN   ) :: v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(IN   ) :: w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(IN   ) :: lambda_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), INTENT(IN   ) :: mu_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      integer(c_int), INTENT(IN   ) :: flag&
            (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

      ! Calculate the trace of the stress tensor (gas phase; m=0)
      CALL CALC_TRD_G(trd_g,u_g,v_g,w_g,flag)

      ! Calculate the cross terms of the stress tensor (gas phase; m=0)
      CALL CALC_TAU_U_G (TAU_U_G,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag)
      CALL CALC_TAU_V_G (TAU_V_G,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag)
      CALL CALC_TAU_W_G (TAU_W_G,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g,flag)

      END SUBROUTINE CALC_TRD_AND_TAU

end module calc_coeff_module
