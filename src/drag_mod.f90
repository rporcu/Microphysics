module drag

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   use param1, only: one, half, zero, small_number, large_number

! Drag model options (see drag_gs for full details)
! default is syam_obrien (may enforce a corrected Umf by defining
! drag_c1 and drag_d1 accordingly)
      character(64) :: drag_type
      integer :: drag_type_enum
      integer,parameter :: syam_obrien=0
      integer,parameter :: gidaspow=1
      integer,parameter :: gidaspow_pcf=2
      integer,parameter :: gidaspow_blend=3
      integer,parameter :: gidaspow_blend_pcf=4
      integer,parameter :: wen_yu=5
      integer,parameter :: wen_yu_pcf=6
      integer,parameter :: koch_hill=7
      integer,parameter :: koch_hill_pcf=8
      integer,parameter :: bvk=9
      integer,parameter :: user_drag=11

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function(s): C_DsxRe                                                C
!  Purpose:                                                            C
!     Calculate single sphere drag correlation multiplied by           C
!     the Reynolds number or                                           C
!     Calculate the single sphere drag correlation                     C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

! Dalla Valle (1948)
!----------------------------------------------------------------->>>
  real(c_real) FUNCTION C_DSXRE_DV(RE)

    IMPLICIT NONE
    real(c_real), intent(IN) :: RE ! Reynolds number

    C_DSXRE_DV = (0.63D0*SQRT(RE) + 4.8D0)**2
    RETURN
  END FUNCTION C_DSXRE_DV

! Schiller and Naumann (1933)
!----------------------------------------------------------------->>>
  real(c_real) FUNCTION C_DS_SN(RE)

    IMPLICIT NONE
    real(c_real), intent(IN) :: RE ! Reynolds number

    C_DS_SN = 24.D0*(1.D0 + 0.15D0*RE**0.687D0)/(RE+SMALL_NUMBER)
    RETURN
  END FUNCTION C_DS_SN
!-----------------------------------------------------------------<<<


!-----------------------------------------------------------------<<<

! Turton and Levenspiel (1986)
!----------------------------------------------------------------->>>
      real(c_real) FUNCTION C_DSXRE_TL(RE)
      IMPLICIT NONE
      real(c_real), intent(IN) :: RE ! Reynolds number

      C_DSXRE_TL = 24.D0*(1.D0 + 0.173D0*RE**0.657D0) + &
         0.413D0*RE**2.09D0/(RE**1.09D0 + 16300.D0)
      RETURN
      END FUNCTION C_DSXRE_TL
!-----------------------------------------------------------------<<<

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_SYAM_OBRIEN                                        C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Syamlal M, O'Brien TJ (1988). International Journal of           C
!        Multiphase Flow 14: 473-481.                                  C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_SYAM_OBRIEN(lDgA,EPg,Mug,ROg,VREL,&
                 DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use constant, only : drag_c1, drag_d1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      real(c_real), intent(OUT) :: ldGA
! gas volume fraction
      real(c_real), intent(IN) :: EPg
! gas laminar viscosity
      real(c_real), intent(IN) :: Mug
! gas density
      real(c_real), intent(IN) :: ROg
! Magnitude of gas-solids relative velocity
      real(c_real), intent(IN) :: VREL
! particle diameter of solids phase M
      real(c_real), intent(IN) :: DPM
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
!     Parameters in the Cluster-effect model
!     PARAMETER (a1 = 250.)   !for G_s = 98 kg/m^2.s
!     PARAMETER (a1 = 1500.)  !for G_s = 147 kg/m^2.s
!     a1 depends upon solids flux.  It has been represented by C(1)
!     defined in the data file.
!     real(c_real), PARAMETER :: A2 = 0.005D0
!     real(c_real), PARAMETER :: A3 = 90.0D0
!     real(c_real), PARAMETER :: RE_C = 5.D0
!     real(c_real), PARAMETER :: EP_C = 0.92D0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Variables which are function of EP_g
      real(c_real) :: A, B
! Ratio of settling velocity of a multiparticle system to
! that of a single particle
      real(c_real) :: V_rm
! Reynolds number
      real(c_real) :: RE
!-----------------------------------------------

      IF(Mug > ZERO) THEN
         RE = DPM*VREL*ROg/Mug
      ELSE
         RE = LARGE_NUMBER
      ENDIF

! Calculate V_rm
      A = EPg**4.14D0
      IF (EPg <= 0.85D0) THEN
         B = drag_c1*EPg**1.28D0
      ELSE
        B = EPg**drag_d1
      ENDIF

      V_RM=HALF*(A-0.06D0*RE+&
           SQRT((3.6D-3)*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A) )

!------------------Begin cluster correction --------------------------
! uncomment the following four lines ...
!       V_RM = V_RM * (ONE + C(1)*&
!                      EXP(-A2*(RE-RE_C)**2 - A3*(EPg-EP_C)**2)* &
!                      RE*(1. - EPg))
!------------------End cluster correction ----------------------------

      lDgA = 0.75D0*Mug*EPg*C_DSXRE_DV(RE/V_RM)/(V_RM*DPM*DPM)

      IF (RE < EPSILON(RE)) lDgA = ZERO

      RETURN

      END SUBROUTINE DRAG_SYAM_OBRIEN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_GIDASPOW                                           C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Ding J, Gidaspow D (1990). AIChE Journal 36: 523-538.            C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_GIDASPOW(lDgA,EPg,Mug,ROg,ROPg,VREL, DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      real(c_real), intent(OUT) :: lDgA
! gas volume fraction
      real(c_real), intent(IN) :: EPg
! gas laminar viscosity
      real(c_real), intent(IN) :: Mug
! gas density
      real(c_real), intent(IN) :: ROg
! gas density*EP_g
      real(c_real), intent(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      real(c_real), intent(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      real(c_real), intent(IN) :: DPM

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      real(c_real) :: RE
! Single sphere drag coefficient
      real(c_real) :: C_d
!-----------------------------------------------

! Note the presence of gas volume fraction in ROPG
      RE = merge(DPM*VREL*ROPg/Mug, LARGE_NUMBER, MUg > ZERO)

! Dense phase
      IF(EPg <= 0.8D0) THEN
         lDgA = 150D0*(ONE-EPg)*Mug / (EPg*DPM**2) + &
                1.75D0*ROg*VREL/DPM
      ELSE
! Dilute phase - EP_g >= 0.8
         IF(RE <= 1000D0)THEN
! this could be replaced with the function C_DS_SN
            C_d = C_DS_SN(RE)
         ELSE
            C_d = 0.44D0
         ENDIF
         lDgA = 0.75D0*C_d*VREL*ROPg*EPg**(-2.65D0) / DPM
      ENDIF

      IF (RE < EPSILON(RE)) lDgA = ZERO

      RETURN

      END SUBROUTINE DRAG_GIDASPOW

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_GIDASPOW_BLEND                                     C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Author: Charles E.A. Finney                        Date: 23-Mar-06  C
!  Reviewer: Sreekanth Pannala                        Date: 24-Mar-06  C
!                                                                      C
!  Literature/Document References:                                     C
!     original source unknown:                                         C
!     Lathouwers D, Bellan J (2000). Proceedings of the 2000 U.S. DOE  C
!        Hydrogen Program Review NREL/CP-570-28890. Available from     C
!     http://www.eere.energy.gov/hydrogenandfuelcells/pdfs/28890k.pdf. C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_GIDASPOW_BLEND(lDgA,EPg,Mug,ROg,ROPg,VREL,&
                 DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use constant, only : PI
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      real(c_real), intent(OUT) :: lDgA
! gas volume fraction
      real(c_real), intent(IN) :: EPg
! gas laminar viscosity
      real(c_real), intent(IN) :: Mug
! gas density
      real(c_real), intent(IN) :: ROg
! gas density*EP_g
      real(c_real), intent(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      real(c_real), intent(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      real(c_real), intent(IN) :: DPM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      real(c_real) :: RE
! Single sphere drag coefficient
      real(c_real) :: C_d
! Gidaspow switch function variables
      real(c_real) :: Ergun
      real(c_real) :: WenYu
      real(c_real) :: PHI_gs
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPM*VREL*ROPg/Mug
      ELSE
         RE = LARGE_NUMBER
      ENDIF

! Dense phase - EP_g <= 0.8
      Ergun = 150D0*(ONE-EPg)*Mug / (EPg*DPM**2) + &
               1.75D0*ROg*VREL/DPM
! Dilute phase - EP_g >= 0.8
      IF(RE <= 1000D0)THEN
         C_d = (24.D0/(RE+SMALL_NUMBER)) * &
           (ONE + 0.15D0 * RE**0.687D0)
      ELSE
         C_d = 0.44D0
      ENDIF
      WenYu = 0.75D0*C_d*VREL*ROPg*EPg**(-2.65D0) / DPM

! Switch function
      PHI_gs = ATAN(150.D0*1.75D0*(EPg - 0.8D0))/PI + 0.5D0

! Blend the models
      lDgA = (1.D0-PHI_gs)*Ergun + PHI_gs*WenYu
      IF (RE < EPSILON(RE)) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_GIDASPOW_BLEND



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_WEN_YU                                             C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Wen CY, Yu YH (1966). Chemical Engineering Progress Symposium    C
!        Series 62: 100-111.                                           C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_WEN_YU(lDgA,EPg,Mug,ROPg,VREL,DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      real(c_real), intent(OUT) :: lDgA
! gas volume fraction
      real(c_real), intent(IN) :: EPg
! gas laminar viscosity
      real(c_real), intent(IN) :: Mug
! gas density*EP_g
      real(c_real), intent(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      real(c_real), intent(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      real(c_real), intent(IN) :: DPM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      real(c_real) :: RE
! Single sphere drag coefficient
      real(c_real) :: C_d
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPM*VREL*ROPg/Mug
      ELSE
         RE = LARGE_NUMBER
      ENDIF

      IF(RE <= 1000.0D0)THEN
         C_d = (24.D0/(RE+SMALL_NUMBER)) * (ONE + 0.15D0*RE**0.687D0)
      ELSE
         C_d = 0.44D0
      ENDIF

      lDgA = 0.75D0 * C_d * VREL * ROPg * EPg**(-2.65D0) / DPM
      IF (RE < EPSILON(RE)) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_WEN_YU


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_KOCH_HILL                                          C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Author: Clay Sutton (Lehigh University)            Date: 14-Jul-04  C
!                                                                      C
!  Revision: 1                                                         C
!  Author: Sofiane Benyahia                           Date: 21-Jan-05  C
!                                                                      C
!  Literature/Document References:                                     C
!     Benyahia S, Syamlal M, O'Brien TJ (2006). Powder Technology      C
!        162: 166-174.                                                 C
!     Hill RJ, Koch DL, Ladd JC (2001). Journal of Fluid Mechanics     C
!        448: 213-241.                                                 C
!     Hill RJ, Koch DL, Ladd JC (2001). Journal of Fluid Mechanics     C
!        448: 243-278.                                                 C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_KOCH_HILL(lDgA,EPg,Mug,ROPg,VREL,&
                 DPM,DPA,PHIS)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      real(c_real), intent(OUT) :: lDgA
! gas volume fraction
      real(c_real), intent(IN) :: EPg
! gas laminar viscosity
      real(c_real), intent(IN) :: Mug
! gas density*EP_g
      real(c_real), intent(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      real(c_real), intent(IN) :: VREL
! particle diameter of solids phase M
      real(c_real), intent(IN) :: DPM
! average particle diameter if pcf otherwise DPM again
      real(c_real), intent(IN) :: DPA
! total solids volume fraction of solids phases
      real(c_real), intent(IN) :: PHIS
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      real(c_real) :: RE
! transition Reynolds numbers
      real(c_real) :: Re_Trans_1, Re_Trans_2
! Stokes Drag Force
      real(c_real) :: F_STOKES
! zero Re function for low Reynolds number
      real(c_real) :: F_0
! inertial function for low Reynolds number
      real(c_real) :: F_1
! zero Re function for high Reynolds number
      real(c_real) :: F_2
! inertial function for high Reynolds number
      real(c_real) :: F_3
! dimensionless drag force F
      real(c_real) :: F
! weighting factor to compute F_0 and F_2
      real(c_real) :: w
!-----------------------------------------------


      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG and factor of 1/2
         RE = 0.5D0*DPA*VREL*ROPG/Mug        ! if pcf DPA otherwise DPM
      ELSE
         RE = LARGE_NUMBER
      ENDIF

      F_STOKES = 18.D0*Mug*EPg*EPg/DPM**2    ! use DPM
      w = EXP(-10.0D0*(0.4D0-phis)/phis)

      IF(phis > 0.01D0 .AND. phis < 0.4D0) THEN
         F_0 = (1.0D0-w) * (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + &
            135.0D0/64.0D0*phis*LOG(phis) + 17.14D0*phis) / &
            (1.0D0 + 0.681D0*phis - 8.48D0*phis*phis + &
            8.16D0*phis**3) + w*10.0D0*phis/(1.0D0-phis)**3
      ELSEIF(phis >= 0.4D0) THEN
         F_0 = 10.0D0*phis/(1.0D0-phis)**3
      ENDIF

      IF(phis > 0.01D0 .AND. phis <= 0.1D0) THEN
        F_1 = dsqrt(2.0D0/phis) / 40.0D0
      ELSE IF(phis > 0.1D0) THEN
        F_1 = 0.11D0 + 5.1D-04 * exp(11.6D0*phis)
      ENDIF

      IF(phis < 0.4D0) THEN
        F_2 = (1.0D0-w) * (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + &
           135.0D0/64.0D0*phis*LOG(phis) + 17.89D0*phis) / &
           (1.0D0 + 0.681D0*phis - 11.03D0*phis*phis + &
           15.41D0*phis**3)+ w*10.0D0*phis/(1.0D0-phis)**3
      ELSE
         F_2 = 10.0D0*phis/(1.0D0-phis)**3
      ENDIF

      IF(phis < 0.0953D0) THEN
         F_3 = 0.9351D0*phis + 0.03667D0
      ELSE
         F_3 = 0.0673D0 + 0.212D0*phis +0.0232D0/(1.0-phis)**5
      ENDIF

      Re_Trans_1 = (F_2 - 1.0D0)/(3.0D0/8.0D0 - F_3)
      Re_Trans_2 = (F_3 + dsqrt(F_3*F_3 - 4.0D0*F_1 &
           *(F_0-F_2))) / (2.0D0*F_1)

      IF(phis <= 0.01D0 .AND. RE <= Re_Trans_1) THEN
         F = 1.0D0 + 3.0D0/8.0D0*RE
      ELSEIF(phis > 0.01D0 .AND. RE <= Re_Trans_2) THEN
         F = F_0 + F_1*RE*RE
      ELSEIF(phis <= 0.01D0 .AND. RE > Re_Trans_1 .OR.   &
         phis >  0.01D0 .AND. RE > Re_Trans_2) THEN
         F = F_2 + F_3*RE
      ELSE
         F = zero
      ENDIF

      lDgA = F * F_STOKES
      IF (RE < EPSILON(RE)) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_KOCH_HILL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_BVK                                                C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Beetstra, van der Hoef, Kuipers, Chem. Eng. Science 62           C
!     (Jan 2007)                                                       C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_BVK(lDgA,EPg,Mug,ROPg,VREL,&
                 DPM,DPA,PHIS)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      real(c_real), intent(OUT) :: lDgA
! gas volume fraction
      real(c_real), intent(IN) :: EPg
! gas laminar viscosity
      real(c_real), intent(IN) :: Mug
! gas density*EP_g
      real(c_real), intent(IN) :: ROPg
! magnitude of gas-solids relative velocity
      real(c_real), intent(IN) :: VREL
! particle diameter of solids phase M or
      real(c_real), intent(IN) :: DPM
! average particle diameter
      real(c_real), intent(IN) :: DPA
! total solids volume fraction of solids phases
      real(c_real), intent(IN) :: PHIS
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      real(c_real) :: RE
! Stokes Drag Force
      real(c_real) :: F_STOKES
! dimensionless drag force F
      real(c_real) :: F
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPA*VREL*ROPg/Mug        ! use DPA
      ELSE
         RE = LARGE_NUMBER
      ENDIF

! eq(9) BVK J. fluid. Mech. 528, 2005
! (this F_Stokes is /= of Koch_Hill by a factor of ep_g)
      F_STOKES = 18D0*Mug*EPg/DPM**2   ! use DPM

      F = 10d0*phis/EPg**2 + EPg**2*(ONE+1.5d0*DSQRT(phis))
      F = F + 0.413d0*RE/(24.d0*EPg**2) * &
             (ONE/EPg + 3d0*EPg*phis + 8.4d0/RE**0.343) / &
             (ONE+10.d0**(3d0*phis)/RE**(0.5+2.d0*phis))

      lDgA = F*F_STOKES
      IF (RE < EPSILON(RE)) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_BVK

END MODULE drag
