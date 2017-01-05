MODULE CORNER_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   use write_error_module, only: write_error

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: get_corner_cells(IER)                                  C
!  Author: M. Syamlal                                 Date: 08-JUL-98  C
!                                                                      C
!  Purpose: Identify wall cells with more than one fulid cell as       C
!           a neighbor.  No heat mass, momentum, or energy transfer    C
!           is allowed to such cells to avoid ambiguity.               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine get_corner_cells(slo,shi,lo,hi,flag)

      use functions    , only: iminus, iplus, jminus, jplus, kminus, kplus
      use funits       , only: dmp_log, unit_log
      use param1       , only: max_ncorn
      use matrix       , only: e, w, s, n, t, b

      IMPLICIT NONE

      integer(c_int), intent(in) :: slo(3), shi(3), lo(3), hi(3)
      integer(c_int), intent(in) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

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
      do k = slo(3),shi(3)
         do j = slo(2),shi(2)
            do i = slo(1),shi(1)

               IF (flag(i,j,k,1)>=100 .AND. flag(i,j,k,1) /= 106 .and. &
                  flag(i,j,k,1) /= 107) then
!----------------------------------------------------------------
                  NUM = 0

                  IF (1.eq.flag(iminus(i,j,k),j,k,1)) then
                     NUM = NUM + 1
                     DIR(W) = .TRUE.
                  ELSE
                     DIR(W) = .FALSE.
                  ENDIF

                  IF (1.eq.flag(iplus(i,j,k),j,k,1)) then
                     NUM = NUM + 1
                     DIR(E) = .TRUE.
                  ELSE
                     DIR(E) = .FALSE.
                  ENDIF

                  IF (1.eq.flag(i,jminus(i,j,k),k,1)) then
                     NUM = NUM + 1
                     DIR(S) = .TRUE.
                  ELSE
                     DIR(S) = .FALSE.
                  ENDIF

                  IF (1.eq.flag(i,jplus(i,j,k),k,1)) then
                     NUM = NUM + 1
                     DIR(N) = .TRUE.
                  ELSE
                     DIR(N) = .FALSE.
                  ENDIF

                  IF (1.eq.flag(i,j,kminus(i,j,k),1)) then
                     NUM = NUM + 1
                     DIR(B) = .TRUE.
                  ELSE
                     DIR(B) = .FALSE.
                  ENDIF

                  IF (1.eq.flag(i,j,kplus(i,j,k),1)) then
                     NUM = NUM + 1
                     DIR(T) = .TRUE.
                  ELSE
                     DIR(T) = .FALSE.
                  ENDIF

                  IF (NUM > 1) THEN

                     NOTCORNER = .TRUE.

!           check for single cell thick internal walls
                     IF (DIR(W) .AND. DIR(E) .OR. DIR(S) .AND. DIR(N) .OR. DIR(T)&
                        .AND. DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF

!           check for corner cells
                     IF (DIR(E)) THEN
                        IF (DIR(N)) THEN
                           IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                        ENDIF
                        IF (DIR(S)) THEN
                           IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                        ENDIF
                        IF (DIR(T)) THEN
                           IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                        ENDIF
                        IF (DIR(B)) THEN
                           IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                        ENDIF
                     ENDIF

                  ENDIF

                  IF (DIR(W)) THEN
                     IF (DIR(N)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                     IF (DIR(S)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                     IF (DIR(T)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                     IF (DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                  ENDIF

                  IF (DIR(N)) THEN
                     IF (DIR(T)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                     IF (DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                  ENDIF

                  IF (DIR(S)) THEN
                     IF (DIR(T)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                     IF (DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                  ENDIF
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
      CONTAINS
!
      subroutine addcorn(NOTCORNER, NCORN)

      USE param1, only: max_ncorn
      USE compar, only: mype
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

      end subroutine addcorn

      end subroutine get_corner_cells
END MODULE CORNER_MODULE
