!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ACCUM_RESID                                            C
!                                                                      C
!  Purpose: Accumulate all the residuals and calculate max_resid       C
!                                                                      C
!                                                                      C
!  Author: S. Pannala                                 Date: 14-Jun-07  C
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
      SUBROUTINE ACCUM_RESID
!
!-----------------------------------------------
!     M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE residual
      USE run
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER          M, NN

!-----------------------------------------------
!     Local variables
!-----------------------------------------------
      INTEGER              LOCAL_INDEX

!
      IF(DEBUG_RESID) Return
!
      LOCAL_INDEX = 0

! Pack the numerators and denominators into one vector for performing single global operation
      DO NN = 2, NRESID
         LOCAL_INDEX = LOCAL_INDEX + 1
         RESID_PACK(LOCAL_INDEX) = NUM_RESID(NN)
         LOCAL_INDEX = LOCAL_INDEX + 1
         RESID_PACK(LOCAL_INDEX) = DEN_RESID(NN)
      ENDDO

      call global_all_sum(RESID_PACK)

! Unpack the numerators and denominators from the global sum vector

      LOCAL_INDEX = 0

      DO NN = 2, NRESID
         LOCAL_INDEX = LOCAL_INDEX + 1
         NUM_RESID(NN) = RESID_PACK(LOCAL_INDEX)
         LOCAL_INDEX = LOCAL_INDEX + 1
         DEN_RESID(NN) = RESID_PACK(LOCAL_INDEX)
      ENDDO

      DO NN = 2, NRESID
         IF (DEN_RESID(NN) > ZERO) THEN
            RESID(NN) = NUM_RESID(NN)/DEN_RESID(NN)
         ELSE IF (NUM_RESID(NN) == ZERO) THEN
            RESID(NN) = ZERO
         ELSE
            RESID(NN) = UNDEFINED
!     WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!     CALL WRITE_ERROR ('ACCUM_RESID', LINE, 1)
         ENDIF
      ENDDO

      RETURN
      END
