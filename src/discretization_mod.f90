!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DISCRETIZATION                                         C
!  Purpose: A collection of functions for higher-order discretization  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-MAR-97  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE discretization

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   IMPLICIT NONE

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      real(c_real) FUNCTION SUPERBEE (PHI_C)
      use param, only: one, half, zero
      IMPLICIT NONE
      real(c_real), INTENT(IN) :: PHI_C

      real(c_real) :: TH

      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region
         TH = PHI_C/(ONE - PHI_C)
         SUPERBEE = HALF*MAX(ZERO,MIN(ONE,2.d0*TH),MIN(2.d0,TH))
      ELSE IF (abs(PHI_C-ONE) < epsilon(ONE)) THEN
         SUPERBEE = ONE
      ELSE                                       !first order upwinding
         SUPERBEE = ZERO
      ENDIF
      RETURN
      END FUNCTION SUPERBEE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      real(c_real) FUNCTION SMART (PHI_C)
      use param, only: zero, half, one
      IMPLICIT NONE
      real(c_real), INTENT(IN) :: PHI_C

      real(c_real) :: TH

      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region
         TH = PHI_C/(ONE - PHI_C)
         SMART = HALF*MAX(ZERO,MIN(4.0D0*TH,0.75D0 + 0.25D0*TH,2.0D0))
      ELSE IF (abs(PHI_C-ONE) < epsilon(ONE)) THEN
         SMART = ONE
      ELSE                                       !first order upwinding
         SMART = ZERO
      ENDIF
      RETURN
      END FUNCTION SMART




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      real(c_real) FUNCTION ULTRA_QUICK (PHI_C, CF)
      use param, only: zero, half, one
      IMPLICIT NONE
      real(c_real), INTENT(IN) :: PHI_C
      real(c_real), INTENT(IN) :: CF

      real(c_real), PARAMETER :: FIVEOSIX = 5.D0/6.D0
      real(c_real) :: TH, OCF

      OCF = MAX(ONE,ONE/MAX(1.D-2,CF))
      IF (PHI_C > ONE) THEN
         ULTRA_QUICK = HALF
      ELSE IF (PHI_C > FIVEOSIX) THEN
         ULTRA_QUICK = ONE
      ELSE IF (PHI_C > 3.D0/(8.D0*OCF - 6.D0)) THEN
         TH = PHI_C/(ONE - PHI_C)
         ULTRA_QUICK = HALF - 0.125D0*(ONE - TH)
      ELSE IF (PHI_C > ZERO) THEN
         TH = PHI_C/(ONE - PHI_C)
         ULTRA_QUICK = (OCF - ONE)*TH
      ELSE
         TH = PHI_C/(ONE - PHI_C)
         ULTRA_QUICK = HALF*TH
      ENDIF
      RETURN
      END FUNCTION ULTRA_QUICK


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      real(c_real) FUNCTION QUICKEST (PHI_C, CF, ODXC, ODXUC, ODXCD)
      use param, only: zero, half, one, small_number
      IMPLICIT NONE
      real(c_real), INTENT(IN) :: PHI_C
      real(c_real), INTENT(IN) :: CF
      real(c_real), INTENT(IN) :: ODXC
      real(c_real), INTENT(IN) :: ODXUC
      real(c_real), INTENT(IN) :: ODXCD

      real(c_real) :: FCF, TH

      IF (PHI_C>ZERO .AND. PHI_C<ONE) THEN       !monotonic region
         FCF = -(ONE - CF*CF)/3.D0
         TH = PHI_C/(ONE - PHI_C + SMALL_NUMBER)
         QUICKEST = HALF*(ONE - CF) + FCF*(ODXC/ODXCD - ODXC*ODXUC*TH/ODXCD**2)
         IF(PHI_C<CF)QUICKEST=MIN(QUICKEST,(ONE/CF-ONE)*PHI_C/(ONE-PHI_C))
         QUICKEST = MAX(ZERO,MIN(ONE,QUICKEST))
      ELSE                                       !first order upwinding
         QUICKEST = ZERO
      ENDIF
      RETURN
      END FUNCTION QUICKEST


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      real(c_real) FUNCTION MUSCL (PHI_C)
      use param, only: zero, half, one
      IMPLICIT NONE
      real(c_real), INTENT(IN) :: PHI_C

      real(c_real) :: TH

      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region
         TH = PHI_C/(ONE - PHI_C)
         MUSCL = HALF*MAX(ZERO,MIN(2.0D0*TH,(ONE + TH)/2.0D0,2.0D0))
      ELSE IF (abs(PHI_C-ONE) < epsilon(ONE)) THEN
         MUSCL = ONE
      ELSE                                       !first order upwinding
         MUSCL = ZERO
      ENDIF
      RETURN
      END FUNCTION MUSCL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      real(c_real) FUNCTION VANLEER (PHI_C)
      use param, only: zero, one
      IMPLICIT NONE
      real(c_real), INTENT(IN) :: PHI_C

      IF (PHI_C>=ZERO .AND. PHI_C<=ONE) THEN     !monotonic region
         VANLEER = PHI_C
      ELSE                                       !first order upwinding
         VANLEER = ZERO
      ENDIF
      RETURN
      END FUNCTION VANLEER


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      real(c_real) FUNCTION MINMOD (PHI_C)
      use param, only: zero, half, one
      IMPLICIT NONE
      real(c_real), INTENT(IN) :: PHI_C

      IF (PHI_C>=ZERO .AND. PHI_C<=ONE) THEN     !monotonic region
         IF (PHI_C >= HALF) THEN                 !central differencing
            MINMOD = HALF
         ELSE                                    !second order upwinding
            MINMOD = HALF*PHI_C/(ONE - PHI_C)
         ENDIF
      ELSE                                       !first order upwinding
         MINMOD = ZERO
      ENDIF

      RETURN
      END FUNCTION MINMOD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Central scheme.                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      real(c_real) FUNCTION CENTRAL_SCHEME ()
      use param, only: half
      IMPLICIT NONE

      CENTRAL_SCHEME = HALF
      RETURN
      END FUNCTION CENTRAL_SCHEME


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      real(c_real) FUNCTION UMIST (PHI_C)
      use param, only: zero, half, one
      IMPLICIT NONE
      real(c_real), INTENT(IN) :: PHI_C

      real(c_real) :: TH

      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region
         TH = PHI_C/(ONE - PHI_C)
         UMIST = HALF*MAX(ZERO,MIN(2.0D0*TH,0.75D0 + 0.25D0*TH,2.0D0))
      ELSE IF (abs(PHI_C-ONE) < EPSILON(ONE)) THEN
         UMIST = ONE
      ELSE                                       !first order upwinding
         UMIST = ZERO
      ENDIF
      RETURN
      END FUNCTION UMIST


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      real(c_real) FUNCTION XSI (V, DWF)
      use param, only: zero, one
      IMPLICIT NONE
      real(c_real), INTENT(IN) :: V
      real(c_real), INTENT(IN) :: DWF

      IF (V >= ZERO) THEN
         XSI = DWF
      ELSE
         XSI = ONE - DWF
      ENDIF
      RETURN
      END FUNCTION XSI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      real(c_real) FUNCTION PHI_C_OF (PHI_U, PHI_C, PHI_D)
      use param, only: zero
      use param, only: equal
      IMPLICIT NONE
      real(c_real), INTENT(IN) :: PHI_U
      real(c_real), INTENT(IN) :: PHI_C
      real(c_real), INTENT(IN) :: PHI_D

      real(c_real) :: DEN

      IF (equal(PHI_D,PHI_U)) THEN
         PHI_C_OF = ZERO
      ELSE
         DEN = PHI_D - PHI_U
         PHI_C_OF = (PHI_C - PHI_U)/DEN
      ENDIF
      RETURN
      END FUNCTION PHI_C_OF

      END MODULE discretization
