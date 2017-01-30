MODULE ADJUST_LEQ_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_LEQ(RESID, LEQ_IT, LEQ_METHOD, LEQI, LEQM, IER) C
!  Purpose: Adjusts liner equation solver method and iterations        C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 23-MAY-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE ADJUST_LEQ(LEQ_ITL, LEQ_METHODL, LEQI, LEQM)

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LEQ_ITL, LEQ_METHODL, LEQI, LEQM
!
!  The adjustment is disabled, because it was adversely affecting species
!  conservation
!      IF (LEQ_ADJUST .AND. RESID<=TOL_RESID*0.1) THEN
!         LEQM = LEQ_METHOD_CONV
!         LEQI = MIN(LEQ_IT_CONV,LEQ_ITL)
!      ELSE
         LEQM = LEQ_METHODL
         LEQI = LEQ_ITL
!      ENDIF
!
      RETURN
      END SUBROUTINE ADJUST_LEQ
END MODULE ADJUST_LEQ_MODULE
