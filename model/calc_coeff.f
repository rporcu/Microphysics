module calc_coeff_module

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
      SUBROUTINE CALC_COEFF_ALL(ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, mu_g, FLAG)

! Global variables:
!-----------------------------------------------------------------------
      use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3

      ! Flag for explcit coupling between the fluid and particles.
      use discretelement, only: DES_EXPLICITLY_COUPLED
      use calc_drag_des_module

      implicit none

      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) ::  p_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)

!-----------------------------------------------------------------------
! FLAG = 0, overwrite the coeff arrays, (e.g. start of a time step)
! FLAG = 1, do not overwrite
      INTEGER, intent(in) :: FLAG
!-----------------------------------------------

      ! Calculate all physical properties, transport properties,
      ! and exchange rates.
      CALL CALC_COEFF(ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, mu_g, 2)

      IF (DES_EXPLICITLY_COUPLED) CALL CALC_DRAG_DES_EXPLICIT(ep_g,u_g,v_g,w_g,ro_g,rop_g,mu_g)

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
      SUBROUTINE CALC_COEFF(ro_g, p_g, ep_g, rop_g, u_g, v_g, w_g, mu_g, pLevel)

      use fld_const, only: ro_g0
      use compar   , only: istart3,iend3,jstart3,jend3,kstart3,kend3
      use discretelement, only: DES_EXPLICITLY_COUPLED
      use discretelement, only: DES_CONTINUUM_COUPLED

      use calc_drag_des_module

      implicit none

      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) ::  p_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: rop_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)

! Dummy arguments
!-----------------------------------------------------------------------
! Level to calculate physical properties.
! 0) Only density
! 1) Everything but density
! 2) All physical properties
      INTEGER, intent(in) :: pLevel
!-----------------------------------------------------------------------
      integer IER

! Calculate physical properties: (density, specific heat, diameter)
      CALL PHYSICAL_PROP(IER, pLevel, ro_g, p_g, ep_g, rop_g, ro_g0)

! Calculate interphase coeffs: (momentum and energy)
      IF (DES_CONTINUUM_COUPLED .AND. .NOT.DES_EXPLICITLY_COUPLED) &
         CALL CALC_DRAG_DES_2FLUID(ep_g,u_g,v_g,w_g,ro_g,mu_g)

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
                                  ep_g,u_g,v_g,w_g,lambda_g,mu_g)

      use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3

      implicit none

      ! Stress tensor cross terms.
      DOUBLE PRECISION, INTENT(INOUT) :: tau_u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: tau_v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: tau_w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: trd_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)

      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: lambda_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: mu_g&
            (istart3:iend3,jstart3:jend3,kstart3:kend3)

      ! Calculate the trace of the stress tensor (gas phase; m=0)
      CALL CALC_TRD_G(trd_g,u_g,v_g,w_g)

      ! Calculate the cross terms of the stress tensor (gas phase; m=0)
      CALL CALC_TAU_U_G (TAU_U_G,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g)
      CALL CALC_TAU_V_G (TAU_V_G,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g)
      CALL CALC_TAU_W_G (TAU_W_G,trd_g,ep_g,u_g,v_g,w_g,lambda_g,mu_g)

      END SUBROUTINE CALC_TRD_AND_TAU

end module calc_coeff_module
