!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_CORNER_CELLS(IER)                                  C
!  Author: M. Syamlal                                 Date: 08-JUL-98  C
!                                                                      C
!  Purpose: Identify wall cells with more than one fulid cell as       C
!           a neighbor.  No heat mass, momentum, or energy transfer    C
!           is allowed to such cells to avoid ambiguity.               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_CORNER_CELLS()

      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE functions, only: iminus, iplus, jminus, jplus, kminus, kplus
      USE functions, only: fluid_at, cyclic_at, wall_at
      USE param1   , only: max_ncorn
      USE funits   , only: dmp_log, unit_log
      USE geometry , only: do_k
      use matrix   , only: e, w, s, n, t, b

      IMPLICIT NONE

!                      Loop index
      INTEGER          L

!
!                      indices
      INTEGER          i, j, k
!
!                      number of faces adjacent to a fluid cell
      INTEGER          NUM
!
!                      fluid face location, whether not a corner
      LOGICAL          dir(-3:3), NotCorner
!
!-----------------------------------------------
!
!                      indices of corner cells
      INTEGER          IJK_CORN (MAX_NCORN)
!
!                      Number of corner cells
      INTEGER          NCORN

      NCORN = 0
!
      do k = kstart3, kend3
         do j = jstart3, jend3
           do i = istart3, iend3

         IF (WALL_AT(i,j,k).AND..NOT.CYCLIC_AT(i,j,k)) THEN
!----------------------------------------------------------------
            NUM = 0
!
            IF (fluid_at(iminus(i,j,k),j,k)) then
               NUM = NUM + 1
               DIR(W) = .TRUE.
            ELSE
               DIR(W) = .FALSE.
            ENDIF
!
            IF (fluid_at(iplus(i,j,k),j,k)) then
               NUM = NUM + 1
               DIR(E) = .TRUE.
            ELSE
               DIR(E) = .FALSE.
            ENDIF
!
            IF (fluid_at(i,jminus(i,j,k),k)) then
               NUM = NUM + 1
               DIR(S) = .TRUE.
            ELSE
               DIR(S) = .FALSE.
            ENDIF
!
            IF (fluid_at(i,jplus(i,j,k),k)) then
               NUM = NUM + 1
               DIR(N) = .TRUE.
            ELSE
               DIR(N) = .FALSE.
            ENDIF
!
            IF (fluid_at(i,j,kminus(i,j,k))) then
               NUM = NUM + 1
               DIR(B) = .TRUE.
            ELSE
               DIR(B) = .FALSE.
            ENDIF
!
            IF (fluid_at(i,j,kplus(i,j,k))) then
               NUM = NUM + 1
               DIR(T) = .TRUE.
            ELSE
               DIR(T) = .FALSE.
            ENDIF
!
            IF (NUM > 1) THEN
!
!
               NOTCORNER = .TRUE.
!
!           check for single cell thick internal walls
               IF (DIR(W) .AND. DIR(E) .OR. DIR(S) .AND. DIR(N) .OR. DIR(T)&
                   .AND. DIR(B)) THEN
!
                  IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
               ENDIF
!
!           check for corner cells
!
               IF (DIR(E)) THEN
!
                  IF (DIR(N)) THEN
                     IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                  ENDIF
!
                  IF (DIR(S)) THEN
                     IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                  ENDIF
!
                  IF (DO_K) THEN
                     IF (DIR(T)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
!
                     IF (DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                  ENDIF
!
               ENDIF
!
               IF (DIR(W)) THEN
!
                  IF (DIR(N)) THEN
                     IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                  ENDIF
!
                  IF (DIR(S)) THEN
                     IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                  ENDIF
!
                  IF (DO_K) THEN
                     IF (DIR(T)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
!
                     IF (DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                  ENDIF
!
               ENDIF
!
               IF (DIR(N)) THEN
!
!
                  IF (DO_K) THEN
                     IF (DIR(T)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
!
                     IF (DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                  ENDIF
!
               ENDIF
!
               IF (DIR(S)) THEN
!
!
                  IF (DO_K) THEN
                     IF (DIR(T)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
!
                     IF (DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                  ENDIF
!
               ENDIF
!
            ENDIF
!
         ENDIF
          end do
        end do
      end do

      IF (NCORN > 0) THEN
            IF(DMP_LOG)WRITE (UNIT_LOG, 1000)
!
         DO L = 1, NCORN
            IF(DMP_LOG)WRITE (UNIT_LOG, 1100) IJK_CORN(L)
         END DO
         IF(DMP_LOG)WRITE (UNIT_LOG, 1300)
      ENDIF
!
      RETURN
!
 1000 FORMAT(/1X,70('*')//' From: Get_Corner_Cells',/&
         ' Warning: The following wall-cells are adjacent to two or',/,&
         ' more fluid-cells.  Mass, momentum, and energy transfer ',/,&
         ' to these wall-cells have been set to zero.',/,&
         '     IJK')
 1100 FORMAT(3X,I6)
!
 1300 FORMAT(/1X,70('*')/)
      END SUBROUTINE GET_CORNER_CELLS
!
      SUBROUTINE ADDCORN(NOTCORNER, NCORN)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE compar
      USE exit_mod, only: mfix_exit
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER NCORN
      LOGICAL NOTCORNER
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(LEN=80) :: LINE
!-----------------------------------------------
!
!                      error message
!
      NCORN = NCORN + 1
      IF (NCORN > MAX_NCORN) THEN
         WRITE (LINE, '(A)') 'Error: Increase MAX_NCORN in param1.inc.'
         CALL WRITE_ERROR ('AddCorn', LINE, 1)
         CALL MFIX_EXIT(myPE)
      ENDIF
!
      NOTCORNER = .FALSE.
!
      RETURN
      END SUBROUTINE ADDCORN
