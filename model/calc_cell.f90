MODULE CALC_CELL_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: CALC_CELL                                          C
!  Purpose: calculate the i, j or k cell index for the corresponding   C
!     x y or z reactor location. the index returned depends on which   C
!     half of the i, j or k cell that the x, y, or z position          C
!     intersects                                                       C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-DEC-91  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_CELL(REACTOR_LOC, D_DIR, N_DIR, CELL_LOC)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only : half
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! the x, y or z location along the axis for which the cell index (i, j
! or k) is to be found
      real(c_real), INTENT(IN) :: REACTOR_LOC
! number of cells in the corresponding direction (IMAX, JMAX, or KMAX)
      INTEGER, INTENT(IN) :: N_DIR
! the cell lengths along the corresponding axis (DX, DY or DZ)
      real(c_real), INTENT(IN) :: D_DIR
! the i, j, or k cell index that corresponds to the x, y or z
! reactor_location (calculated value)
      INTEGER, INTENT(INOUT) :: CELL_LOC
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop counter
      INTEGER :: LC
! start and end coordinate for cell
      real(c_real) :: CELL_START, CELL_END
!-----------------------------------------------

      CELL_LOC = -1
      CELL_START = 0.0d0
      DO LC = 2, N_DIR + 1
         CELL_END = CELL_START + D_DIR
         IF (REACTOR_LOC <= CELL_START + HALF*D_DIR) THEN
            CELL_LOC = LC - 1
            RETURN
         ELSEIF (REACTOR_LOC <= CELL_END + HALF*D_DIR) THEN
            CELL_LOC = LC
            RETURN
         ENDIF
         CELL_START = CELL_END
      ENDDO
      RETURN
      END SUBROUTINE CALC_CELL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: CALC_LOC                                           C
!  Purpose: calculate the x, y, or z position corresponding to the     C
!     given i, j or k cell index. this call returns the x, y or z      C
!     position corresponding to the east, north or top face of the     C
!     cell with the given i, j or k index.                             C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-DEC-91  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_LOC(D_DIR, CELL_LOC, REACTOR_LOC)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! the i, j, or k cell index that corresponds to the x, y or z
! reactor_location to be found
      INTEGER, INTENT(IN) :: CELL_LOC
! the cell lengths along the corresponding axis (DX, DY or DZ)
      real(c_real), INTENT(IN) :: D_DIR
! the x, y or z location along the axis that corresponds to the i, j
! k cell index  (calculated value)
      real(c_real), INTENT(INOUT) :: REACTOR_LOC
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop counter
      INTEGER :: LC
!-----------------------------------------------

      REACTOR_LOC = 0.0d0
      LC = 2
      IF (CELL_LOC - 1 > 0) THEN
         REACTOR_LOC = REACTOR_LOC + D_DIR*(CELL_LOC-2)
         LC = CELL_LOC + 1
      ENDIF
      RETURN
      END SUBROUTINE CALC_LOC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_CELL_INTERSECT                                     !
!  Author: J.Musser                                   Date: 23-APR-14  !
!                                                                      !
!  Purpose: Calculate the cell index that intersects LOC. Unlike the   !
!  base routine (CALC_CELL), this routine does not shift the index     !
!  based on which half of the cell the point intersects.               !
!                                                                      !
!  Comment: This is a brute force approach and should not be called    !
!  from within any critical routines/loops.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_CELL_INTERSECT(LOC, D_DIR, N_DIR, CELL)

      IMPLICIT NONE

! Passed Arguments:
!---------------------------------------------------------------------//
! Point to check for intersection.
      real(c_real), INTENT(in) :: LOC
! Number of cells in this direction (IMAX,JMAX,KMAX)
      INTEGER, INTENT(in) :: N_DIR
! Cell lengths (DX,DY,DZ)
      real(c_real), INTENT(IN) :: D_DIR
! Cell indices corresponding to LOC
      INTEGER, INTENT(out) :: CELL

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: LC
! Start/End coordinates for cell
      real(c_real) :: CELL_START, CELL_END
!......................................................................!

     CELL = -1

      CELL_START = 0.0d0
      DO LC=2, N_DIR+1
         CELL_END = CELL_START + D_DIR
         IF(CELL_START <= LOC .AND. LOC <= CELL_END) THEN
            CELL = LC
            RETURN
         ENDIF
         CELL_START=CELL_END
      ENDDO

      RETURN
      END SUBROUTINE CALC_CELL_INTERSECT
END MODULE CALC_CELL_MODULE
