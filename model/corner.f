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

      USE compar
      USE functions
      USE funits
      USE geometry
      use matrix, only: e, w, s, n, t, b
      USE param
      USE param1
      USE physprop

      IMPLICIT NONE

!                      Loop index
      INTEGER          L

!
!                      indices
      INTEGER          IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP
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
!                      IJK indices of corner cells
      INTEGER          IJK_CORN (MAX_NCORN)
!
!                      Number of corner cells
      INTEGER          NCORN

      NCORN = 0
!
      do k = kstart3, kend3
         do j = jstart3, jend3
           do i = istart3, iend3

           ijk = funijk(i,j,k)

         IF (WALL_AT(IJK).AND..NOT.CYCLIC_AT(IJK)) THEN
!
!----------------------------------------------------------------
            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IPJK = FUNIJK(iplus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IJPK = FUNIJK(i,jplus(i,j,k),k)
            IJKM = FUNIJK(i,j,kminus(i,j,k))
            IJKP = FUNIJK(i,j,kplus(i,j,k))
!----------------------------------------------------------------
            NUM = 0
!
            IF (fluid_cell(iminus(i,j,k),j,k)) then
               NUM = NUM + 1
               DIR(W) = .TRUE.
            ELSE
               DIR(W) = .FALSE.
            ENDIF
!
            IF (fluid_cell(iplus(i,j,k),j,k)) then
               NUM = NUM + 1
               DIR(E) = .TRUE.
            ELSE
               DIR(E) = .FALSE.
            ENDIF
!
            IF (fluid_cell(i,jminus(i,j,k),k)) then
               NUM = NUM + 1
               DIR(S) = .TRUE.
            ELSE
               DIR(S) = .FALSE.
            ENDIF
!
            IF (fluid_cell(i,jplus(i,j,k),k)) then
               NUM = NUM + 1
               DIR(N) = .TRUE.
            ELSE
               DIR(N) = .FALSE.
            ENDIF
!
            IF (fluid_cell(i,j,kminus(i,j,k))) then
               NUM = NUM + 1
               DIR(B) = .TRUE.
            ELSE
               DIR(B) = .FALSE.
            ENDIF
!
            IF (fluid_cell(i,j,kplus(i,j,k))) then
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
            IJK = IJK_CORN(L)
            IF(DMP_LOG)WRITE (UNIT_LOG, 1100) IJK
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
         CALL MFIX_EXIT()
      ENDIF
!
      NOTCORNER = .FALSE.
!
      RETURN
      END SUBROUTINE ADDCORN
