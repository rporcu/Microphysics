module drag

   use amrex_fort_module, only : rt => amrex_real
   use iso_c_binding , only: c_int

   use param, only: one, half, zero, small_number, large_number

   integer,parameter :: invalid_drag=-1
   integer,parameter :: user_drag=0
   integer,parameter :: wen_yu=1
   integer,parameter :: gidaspow=2
   integer,parameter :: bvk2=3
   integer,parameter :: wen_yu_pcf=21
   integer,parameter :: gidaspow_pcf=22

   character(64) :: drag_type
   integer :: drag_type_enum = invalid_drag

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
  real(rt) FUNCTION C_DSXRE_DV(RE)

    IMPLICIT NONE
    real(rt), intent(IN) :: RE ! Reynolds number

    C_DSXRE_DV = (0.63D0*SQRT(RE) + 4.8D0)**2
    RETURN
  END FUNCTION C_DSXRE_DV

! Schiller and Naumann (1933)
!----------------------------------------------------------------->>>
  real(rt) FUNCTION C_DS_SN(RE)

    IMPLICIT NONE
    real(rt), intent(IN) :: RE ! Reynolds number

    C_DS_SN = 24.D0*(1.D0 + 0.15D0*RE**0.687D0)/(RE+SMALL_NUMBER)
    RETURN
  END FUNCTION C_DS_SN
!-----------------------------------------------------------------<<<


!-----------------------------------------------------------------<<<

! Turton and Levenspiel (1986)
!----------------------------------------------------------------->>>
      real(rt) FUNCTION C_DSXRE_TL(RE)
      IMPLICIT NONE
      real(rt), intent(IN) :: RE ! Reynolds number

      C_DSXRE_TL = 24.D0*(1.D0 + 0.173D0*RE**0.657D0) + &
         0.413D0*RE**2.09D0/(RE**1.09D0 + 16300.D0)
      RETURN
      END FUNCTION C_DSXRE_TL
!-----------------------------------------------------------------<<<

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
      real(rt), intent(OUT) :: lDgA
! gas volume fraction
      real(rt), intent(IN) :: EPg
! gas laminar viscosity
      real(rt), intent(IN) :: Mug
! gas density*EP_g
      real(rt), intent(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      real(rt), intent(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      real(rt), intent(IN) :: DPM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      real(rt) :: RE
! Single sphere drag coefficient
      real(rt) :: C_d
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
!  Subroutine: DRAG_GIDASPOW                                           C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Author: Charles E.A. Finney                        Date: 23-Mar-06  C
!  Reviewer: Sreekanth Pannala                        Date: 24-Mar-06  C
!                                                                      C
!  Literature/Document References:                                     C
!     Ding J, Gidaspow D (1990). AIChE Journal 36: 523-538.            C
!     original source unknown:                                         C
!     Lathouwers D, Bellan J (2000). Proceedings of the 2000 U.S. DOE  C
!        Hydrogen Program Review NREL/CP-570-28890. Available from     C
!     http://www.eere.energy.gov/hydrogenandfuelcells/pdfs/28890k.pdf. C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_GIDASPOW(lDgA,EPg,Mug,ROg,ROPg,VREL,DPM)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use amrex_constants_module, only : M_PI
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      real(rt), intent(OUT) :: lDgA
! gas volume fraction
      real(rt), intent(IN) :: EPg
! gas laminar viscosity
      real(rt), intent(IN) :: Mug
! gas density
      real(rt), intent(IN) :: ROg
! gas density*EP_g
      real(rt), intent(IN) :: ROPg
! Magnitude of gas-solids relative velocity
      real(rt), intent(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      real(rt), intent(IN) :: DPM
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      real(rt) :: RE
! Single sphere drag coefficient
      real(rt) :: C_d
! Gidaspow switch function variables
      real(rt) :: Ergun
      real(rt) :: WenYu
      real(rt) :: PHI_gs
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
      PHI_gs = ATAN(150.D0*1.75D0*(EPg - 0.8D0))/M_PI + 0.5D0

! Blend the models
      lDgA = (1.D0-PHI_gs)*Ergun + PHI_gs*WenYu
      IF (RE < EPSILON(RE)) lDgA = ZERO

      RETURN
      END SUBROUTINE DRAG_GIDASPOW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_BVK2                                               C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Literature/Document References:                                     C
!     Beetstra, van der Hoef, Kuipers, Chem. Eng. Science 62           C
!     (Jan 2007)                                                       C
!                                                                      C
!     Tang, Peters, Kuipers, Kriebitzsch, van der Hoef, AIChEJ,        C
!     61(2) (Feb 2015)                                                 C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_BVK2(lDgA,EPg,Mug,ROPg,VREL,&
                 DPM,DPA,PHIS)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      real(rt), intent(OUT) :: lDgA
! gas volume fraction
      real(rt), intent(IN) :: EPg
! gas laminar viscosity
      real(rt), intent(IN) :: Mug
! gas density*EP_g
      real(rt), intent(IN) :: ROPg
! magnitude of gas-solids relative velocity
      real(rt), intent(IN) :: VREL
! particle diameter of solids phase M or
      real(rt), intent(IN) :: DPM
! average particle diameter
      real(rt), intent(IN) :: DPA
! total solids volume fraction of solids phases
      real(rt), intent(IN) :: PHIS
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      real(rt) :: RE
! Stokes Drag Force
      real(rt) :: F_STOKES
! dimensionless drag force F
      real(rt) :: F
! 1.0d0/EPg**4
      real(rt) :: oEPgfour
!-----------------------------------------------

      IF(Mug > ZERO) THEN
! Note the presence of gas volume fraction in ROPG
         RE = DPA*VREL*ROPg/Mug        ! use DPA
      ELSE
         RE = LARGE_NUMBER
      ENDIF

      if(RE > epsilon(0.0_rt)) then

         oEPgfour = 1.0d0/EPg**4
      
! eq(9) BVK J. fluid. Mech. 528, 2005
! (this F_Stokes is /= of Koch_Hill by a factor of ep_g)
         F_STOKES = 18D0*Mug*EPg/DPM**2   ! use DPM

         F = 10d0*phis/EPg**2 + EPg**2*(ONE+1.5d0*DSQRT(phis))

         F = F + RE*(0.11d0*phis*(1.0d0+phis) - 4.56d-3*oEPgfour + &
                RE**-0.343d0*(0.169d0*EPg + 6.44d-2*oEPgfour))       

!rm         F = F + 0.413d0*RE/(24.d0*EPg**2) * &
!rm              (ONE/EPg + 3d0*EPg*phis + 8.4d0/RE**0.343) / &
!rm              (ONE+10.d0**(3d0*phis)/RE**(0.5+2.d0*phis))

         lDgA = F*F_STOKES

      else

         lDga = ZERO
      endif


      RETURN
      END SUBROUTINE DRAG_BVK2

END MODULE drag
