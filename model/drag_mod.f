!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: drag                                                   C
!  Purpose: Common block containing drag arrays                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE drag

! Gas-solids drag
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  F_gs

! Solids-solids drag
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  F_ss

! Off diagonal friction coefficient in HYS drag relation
  DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  beta_ij

CONTAINS

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
  DOUBLE PRECISION FUNCTION C_DSXRE_DV(RE)
    USE param
    USE param1
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RE ! Reynolds number

    C_DSXRE_DV = (0.63D0*SQRT(RE) + 4.8D0)**2
    RETURN
  END FUNCTION C_DSXRE_DV

! Schiller and Naumann (1933)
!----------------------------------------------------------------->>>
  DOUBLE PRECISION FUNCTION C_DS_SN(RE)
    USE param
    USE param1
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RE ! Reynolds number

    C_DS_SN = 24.D0*(1.D0 + 0.15D0*RE**0.687D0)/(RE+SMALL_NUMBER)
    RETURN
  END FUNCTION C_DS_SN
!-----------------------------------------------------------------<<<


!-----------------------------------------------------------------<<<

! Turton and Levenspiel (1986)
!----------------------------------------------------------------->>>
      DOUBLE PRECISION FUNCTION C_DSXRE_TL(RE)
      USE param
      USE param1
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: RE ! Reynolds number

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
      USE constant, only : drag_c1, drag_d1
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: ldGA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: DPM
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
!     Parameters in the Cluster-effect model
!     PARAMETER (a1 = 250.)   !for G_s = 98 kg/m^2.s
!     PARAMETER (a1 = 1500.)  !for G_s = 147 kg/m^2.s
!     a1 depends upon solids flux.  It has been represented by C(1)
!     defined in the data file.
!     DOUBLE PRECISION, PARAMETER :: A2 = 0.005D0
!     DOUBLE PRECISION, PARAMETER :: A3 = 90.0D0
!     DOUBLE PRECISION, PARAMETER :: RE_C = 5.D0
!     DOUBLE PRECISION, PARAMETER :: EP_C = 0.92D0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Variables which are function of EP_g
      DOUBLE PRECISION :: A, B
! Ratio of settling velocity of a multiparticle system to
! that of a single particle
      DOUBLE PRECISION :: V_rm
! Reynolds number
      DOUBLE PRECISION :: RE
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

      IF (RE == ZERO) lDgA = ZERO

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
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient
      DOUBLE PRECISION :: C_d
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

      IF (RE == ZERO) lDgA = ZERO

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
      USE param
      USE param1
      USE constant, only : PI
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient
      DOUBLE PRECISION :: C_d
! Gidaspow switch function variables
      DOUBLE PRECISION :: Ergun
      DOUBLE PRECISION :: WenYu
      DOUBLE PRECISION :: PHI_gs
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
      IF (RE == ZERO) lDgA = ZERO

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
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient
      DOUBLE PRECISION :: C_d
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
      IF (RE == ZERO) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_WEN_YU



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_DRAG_IGCI                                       C
!  Purpose: Calculate subgrid correction to the gas-solids drag        C
!           coefficient developed by Wen-Yu                            C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Revision: 1                                                         C
!  Purpose: Minor changes, e.g., fix inconsistenty with analogous      C
!     calls in des_drag_gp and with new variable density feature       C
!  Author: Janine Carney, June 2013                                    C
!                                                                      C
!  Literature/Document References:                                     C
!     Igci, Y., Pannala, S., Benyahia, S., & Sundaresan S.,            C
!        Validation studies on filtered model equations for gas-       C
!        particle flows in risers, Industrial & Engineering Chemistry  C
!        Research, 2012, 51(4), 2094-2103                              C
!                                                                      C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE SUBGRID_DRAG_IGCI(lDgA,EPg,Mug,ROg,DPM,ROs,IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE run, only : filter_size_ratio, SUBGRID_WALL
      USE constant, only : GRAVITY
      USE geometry, only : VOL,AXY,DO_K
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(INOUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! particle diameter of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: DPM
! particle density of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: ROs
! current cell index
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! factor to correct the drag for subgrid domain effects
      DOUBLE PRECISION :: F_Subgrid
! factor to correct the drag for subgrid domain effects arising from
! wall
      DOUBLE PRECISION :: F_SubGridWall
! particle terminal settling velocity from stokes' formulation
      DOUBLE PRECISION :: vt
! filter size which is a function of each grid cell volume
      DOUBLE PRECISION :: filtersize
! inverse Froude number, or dimensionless filter size
      DOUBLE PRECISION :: Inv_Froude
! total solids volume fraction
      DOUBLE PRECISION :: EPs
! Variables for Igci model
      DOUBLE PRECISION :: GG_phip, h_phip, h_phip2, c_function,&
                          f_filter
!-----------------------------------------------

! initialize
      F_Subgrid = ONE
      F_SubgridWall = ONE

! particle terminal settling velocity: vt = g*d^2*(Rho_s - Rho_g) / 18 * Mu_g
      vt = GRAVITY*DPM*DPM*(ROs - ROg) / (18.0d0*Mug)
! filter size calculation for each specific gridcell volume
      IF(DO_K) THEN
         filtersize = filter_size_ratio * (VOL(IJK)**(ONE/3.0d0))
      ELSE
         filtersize = filter_size_ratio * DSQRT(AXY(IJK))
      ENDIF

! dimensionless inverse of Froude number
      IF(ABS(vt) > SMALL_NUMBER) THEN
         Inv_Froude =  filtersize * GRAVITY / vt**2
      ELSE
         Inv_Froude =  LARGE_NUMBER
      ENDIF

! total solids volume fraction
      EPs = ONE - EPg

      IF (EPs .LT. 0.0012d0) THEN
         h_phip = 2.7d0*(EPs**0.234)
      ELSEIF (EPs .LT. 0.014d0) THEN
         h_phip = -0.019d0/(EPs**0.455) + 0.963d0
      ELSEIF (EPs .LT. 0.25d0) THEN
         h_phip = 0.868d0*EXP((-0.38*EPs)) - &
            0.176d0*EXP((-119.2*EPs))
      ELSEIF (EPs .LT. 0.455d0) THEN
         h_phip = -4.59d-5*EXP((19.75*EPs)) + &
            0.852d0*EXP((-0.268*EPs))
      ELSEIF (EPs .LE. 0.59d0) THEN
         h_phip = (EPs - 0.59d0) * (-1501.d0*(EPs**3) + &
            2203.d0*(EPs**2) - 1054.d0*EPs + 162.d0)
      ELSE
         h_phip=ZERO
      ENDIF

      IF (EPs .LT. 0.18d0) THEN
          GG_phip = (EPs**0.24)*(1.48d0 + EXP(-18.0*EPs))
      ELSE
          GG_phip = ONE
      ENDIF

! a filter function needed in Igci Filtered/subgrid Model [dimensionless]
      f_filter = (Inv_Froude**1.6) / ((Inv_Froude**1.6)+0.4d0)
      h_phip2=h_phip*GG_phip
      c_function=-h_phip2*f_filter
      F_Subgrid =(ONE + c_function)

      IF (SUBGRID_WALL) THEN
         CALL SUBGRID_DRAG_WALL(F_SubgridWall,vt,IJK)
      ENDIF

      lDgA = F_SubgridWall*F_Subgrid * lDgA

      RETURN
      END SUBROUTINE SUBGRID_DRAG_IGCI



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_DRAG_MILIOLI                                    C
!  Purpose: Calculate subgrid correction to the gas-solids drag        C
!           coefficient developed by Wen-Yu                            C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Revision: 1                                                         C
!  Purpose: Minor changes, e.g., fix inconsistenty with analogous      C
!     calls in des_drag_gp and with new variable density feature       C
!  Author: Janine Carney, June 2013                                    C
!                                                                      C
!  Literature/Document References:                                     C
!     Milioli, C. C., et al., Filtered two-fluid models of fluidized   C
!        gas-particle flows: new constitutive relations, AICHE J,      C
!        doi: 10.1002/aic.14130                                        C
!                                                                      C
!  Comments:                                                           C
!     Still needs to be reviewed for accuracy with source material     C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE SUBGRID_DRAG_MILIOLI(lDgA,EPg,Mug,ROg,VREL,DPM,ROs,&
         IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE run, only : filter_size_ratio, SUBGRID_WALL
      USE constant, only : GRAVITY
      USE geometry, only : VOL,AXY,DO_K
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(INOUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: DPM
! particle density of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: ROs
! current cell index
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! factor to correct the drag for subgrid domain effects
      DOUBLE PRECISION :: F_Subgrid
! factor to correct the drag for subgrid domain effects arising from
! wall
      DOUBLE PRECISION :: F_SubGridWall
! particle terminal settling velocity from stokes' formulation
      DOUBLE PRECISION :: vt
! filter size which is a function of each grid cell volume
      DOUBLE PRECISION :: filtersize
! inverse Froude number, or dimensionless filter size
      DOUBLE PRECISION :: Inv_Froude
! dimensionless slip velocity = VREL/vt
      DOUBLE PRECISION :: vslip
! total solids volume fraction
      DOUBLE PRECISION :: EPs
! Variables for Milioli model
      DOUBLE PRECISION :: h1, henv, hlin
!-----------------------------------------------

! initialize
      F_Subgrid = ONE
      F_SubgridWall = ONE

! particle terminal settling velocity: vt = g*d^2*(Rho_s - Rho_g) / 18 * Mu_g
      vt = GRAVITY*DPM*DPM*(ROs - ROg) / (18.0d0*Mug)
! filter size calculation for each specific gridcell volume
      IF(DO_K) THEN
         filtersize = filter_size_ratio * (VOL(IJK)**(ONE/3.0d0))
      ELSE
         filtersize = filter_size_ratio * DSQRT(AXY(IJK))
      ENDIF
! dimensionless inverse of Froude number
      IF(ABS(vt) > SMALL_NUMBER) THEN
         Inv_Froude =  filtersize * GRAVITY / vt**2
      ELSE
         Inv_Froude =  LARGE_NUMBER
      ENDIF
! total solids volume fractionn
      EPs = ONE - EPg
! dimensionless slip velocity between gas and solids phase M
      Vslip = VREL / vt

      IF (Inv_Froude .LE. 1.028d0) THEN
         h1 = (1.076d0 + 0.12d0*Vslip - (0.02d0/(Vslip+0.01d0)))*EPs + &
            (0.084d0 + 0.09d0*Vslip - (0.01d0/(0.1d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.53d0) THEN
            henv = (6.8d0*(ONE+EPs)*(EPs**0.3)) / &
               (10.d0*(EPs**1.5) + 5.d0)
         ELSEIF (EPs .GT. 0.53d0 .AND. EPs .LE. 0.65d0) THEN
            henv = (2.23d0*((0.65d0-EPs)**(0.45))) / &
               ((ONE/EPs)-ONE)
         ELSEIF (EPs .GT. 0.65d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 1.028d0 .AND. &
              Inv_Froude .LE. 2.056d0) THEN
         h1 = (1.268d0 - (0.2d0*Vslip) + (0.14d0/(Vslip+0.01d0)))*EPs + &
            (0.385d0 + 0.09d0*Vslip - (0.05d0/(0.2d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.53d0) THEN
            henv = (8.6d0*(ONE+EPs)*(EPs**0.2)) / (10.d0*EPs + 6.3d0)
         ELSEIF (EPs .GT. 0.53d0 .AND. EPs .LE. 0.65d0) THEN
            henv = (0.423d0*((0.65d0-EPs)**0.3)) / (ONE-(EPs**0.4))
         ELSEIF (EPs .GT. 0.65d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 2.056d0 .AND. &
              Inv_Froude .LE. 4.112d0) THEN
         h1 = ((0.018d0*Vslip + 0.1d0)/(0.14d0*Vslip + 0.01d0))*EPs + &
            (0.9454d0 - (0.09d0/(0.2d0*Vslip + 0.01d0)))
         IF (EPs .LE. 0.5d0) THEN
            henv = (7.9d0*(ONE+EPs)*(EPs**0.2)) / &
               (10.d0*(EPs**0.9) + 5.d0)
         ELSEIF (EPs .GT. 0.5d0 .AND. EPs .LE. 0.63d0) THEN
            henv = (0.705d0*((0.63d0-EPs)**0.3)) / (ONE-(EPs**0.7))
         ELSEIF (EPs .GT. 0.63d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 4.112d0 .AND. &
              Inv_Froude .LE. 8.224d0) THEN
         h1 = ((0.05d0*Vslip+0.3d0)/(0.4d0*Vslip+0.06d0))*EPs + &
            (0.9466d0 - (0.05d0/(0.11d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.45d0) THEN
            henv = (7.9d0*(ONE+EPs)*(EPs**0.2)) / &
               ((10.d0*(EPs**0.6)) + 3.6d0)
         ELSEIF (EPs .GT. 0.45d0 .AND. EPs .LE. 0.57d0) THEN
            henv = (0.78d0*((0.57d0-EPs)**0.2)) / (ONE-(EPs**0.9))
         ELSEIF (EPs .GT. 0.57d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 8.224d0 .AND. &
              Inv_Froude .LE. 12.336d0) THEN
         h1 = ((1.3d0*Vslip+2.2d0)/(5.2d0*Vslip+0.07d0))*EPs + &
            (0.9363d0-(0.11d0/(0.3d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.35d0) THEN
            henv = (7.6d0*(ONE+EPs)*(EPs**0.2)) / &
               ((10.d0*(EPs**0.6)) + 3.3d0)
         ELSEIF (EPs .GT. 0.35d0 .AND. EPs .LE. 0.55d0) THEN
            henv = (0.81d0*((0.55d0-EPs)**0.3)) / (ONE-(EPs**0.7))
         ELSEIF (EPs .GT. 0.55d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 12.336d0 .AND. &
              Inv_Froude .LE. 16.448d0) THEN
         h1 = ((2.6d0*Vslip+4.d0)/(10.d0*Vslip+0.08d0))*EPs + &
            (0.926d0-(0.17d0/(0.5d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.25d0) THEN
            henv = (8.4d0*(ONE+EPs)*(EPs**0.2)) / &
               ((10.d0*(EPs**0.5)) + 3.3d0)
         ELSEIF (EPs .GT. 0.25d0 .AND. EPs .LE. 0.52d0) THEN
            henv = (1.01d0*((0.52d0-EPs)**0.03))/(ONE-(EPs**0.9))
         ELSEIF (EPs .GT. 0.52d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 16.448d0 .AND. &
              Inv_Froude .LE. 20.56d0) THEN
         h1 = ((2.5d0*Vslip+4.d0)/(10.d0*Vslip+0.08d0))*EPs + &
            (0.9261d0-(0.17d0/(0.5d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.25d0) THEN
            henv = (8.4d0*(ONE+EPs)*(EPs**0.2)) / &
               ((10.d0*(EPs**0.5)) + 3.3d0)
         ELSEIF (EPs .GT. 0.25d0 .AND. EPs .LE. 0.52d0) THEN
            henv = (1.065d0*((0.52d0-EPs)**0.3))/(ONE-EPs)
         ELSEIF (EPs .GT.0.52d0) THEN
            henv=ZERO
         ENDIF

      ELSEIF (Inv_Froude .GT. 20.56d0) THEN
         h1 = ((1.6d0*Vslip+4.d0)/(7.9d0*Vslip+0.08d0))*EPs + &
            (0.9394d0 - (0.22d0/(0.6d0*Vslip+0.01d0)))
         IF (EPs .LE. 0.25d0) THEN
            henv = (9.d0*(ONE+EPs)*(EPs**0.15)) / &
               (10.d0*(EPs**0.45) + 4.2d0)
         ELSEIF (EPs .GT. 0.25d0 .AND. EPs .LE. 0.52d0) THEN
            henv = (0.91d0*((0.52d0-EPs)**0.4))/(ONE-(EPs**0.6))
         ELSEIF (EPs .GT. 0.52d0) THEN
            henv=ZERO
         ENDIF
      ENDIF

      IF (h1 .GT. ZERO) THEN
         hlin=h1
      ELSE
         hlin=ZERO
      ENDIF

      IF (Inv_Froude .LT. 1.028d0) THEN
! for very small filtered size, the drag wont be changed:
! F_Subgrid = 1.0 - H where H = 0.0
         F_Subgrid = ONE
      ELSE
! MIN(henv,hlin) is H in Milioli paper, 2013
         F_Subgrid = ONE - MIN(henv,hlin)
      ENDIF

! Filtered drag = (1 - H)*Microscopic_drag; it is strange Milioli takes EPs
!     F_Subgrid = EPs*(ONE-hmili)

      IF (SUBGRID_WALL) THEN
         CALL SUBGRID_DRAG_WALL(F_SubgridWall,vt,IJK)
      ENDIF

      lDgA = F_SubgridWall*F_Subgrid * lDgA


      RETURN
      END SUBROUTINE SUBGRID_DRAG_MILIOLI



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_DRAG_WALL                                       C
!  Purpose: Calculate subgrid correction arising from wall to the      C
!     gas-solids drag coefficient                                      C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Revision: 1                                                         C
!  Author: Janine Carney, June 2013                                    C
!                                                                      C
!  Literature/Document References:                                     C
!     Igci, Y., and Sundaresan, S., Verification of filtered two-      C
!        fluid models for gas-particle flows in risers, AICHE J.,      C
!        2011, 57 (10), 2691-2707.                                     C
!                                                                      C
!  Comments: Currently only valid for free-slip wall but no checks     C
!     are made to ensure user has selected free-slip wall when this    C
!     option is invoked                                                C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE SUBGRID_DRAG_WALL(lSubgridWall,vt,IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE constant, only : GRAVITY
      USE cutcell, only : DWALL
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! factor to correct the drag for subgrid domain effects arising from
! wall
      DOUBLE PRECISION, INTENT(OUT) :: lSubGridWall
! particle terminal settling velocity from stokes' formulation
      DOUBLE PRECISION, INTENT(IN) :: vt
! current cell index
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! values are only correct for FREE-Slip walls
      DOUBLE PRECISION, PARAMETER :: a22=6.0d0, b22=0.295d0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! dimensionless distance to wall
      DOUBLE PRECISION :: x_d
!-----------------------------------------------
! initialize
      lSubgridWall = ONE

! dimensionless distance to the Wall
      x_d = DWALL(IJK) * GRAVITY / vt**2

! decrease exponentionally away from the wall
! more complex model could be implemented with JJ wall model
      lSubgridWall = ONE / ( ONE + a22 * (EXP(-b22*x_d)) )

      RETURN
      END SUBROUTINE SUBGRID_DRAG_WALL



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
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: DPM
! average particle diameter if pcf otherwise DPM again
      DOUBLE PRECISION, INTENT(IN) :: DPA
! total solids volume fraction of solids phases
      DOUBLE PRECISION, INTENT(IN) :: PHIS
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! transition Reynolds numbers
      DOUBLE PRECISION :: Re_Trans_1, Re_Trans_2
! Stokes Drag Force
      DOUBLE PRECISION :: F_STOKES
! zero Re function for low Reynolds number
      DOUBLE PRECISION :: F_0
! inertial function for low Reynolds number
      DOUBLE PRECISION :: F_1
! zero Re function for high Reynolds number
      DOUBLE PRECISION :: F_2
! inertial function for high Reynolds number
      DOUBLE PRECISION :: F_3
! dimensionless drag force F
      DOUBLE PRECISION :: F
! weighting factor to compute F_0 and F_2
      DOUBLE PRECISION :: w
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
      IF (RE == ZERO) lDgA = ZERO

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
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density*EP_g
      DOUBLE PRECISION, INTENT(IN) :: ROPg
! magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
      DOUBLE PRECISION, INTENT(IN) :: DPM
! average particle diameter
      DOUBLE PRECISION, INTENT(IN) :: DPA
! total solids volume fraction of solids phases
      DOUBLE PRECISION, INTENT(IN) :: PHIS
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Stokes Drag Force
      DOUBLE PRECISION :: F_STOKES
! dimensionless drag force F
      DOUBLE PRECISION :: F
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
      IF (RE == ZERO) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_BVK

END MODULE drag
