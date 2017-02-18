MODULE utilities

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   IMPLICIT NONE

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  function: mfix_isnan                                                !
!  Purpose: check whether argument is NAN                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      LOGICAL FUNCTION mfix_isnan(x)

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      real(c_real) :: x
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      CHARACTER(LEN=80) :: notnumber
!-----------------------------------------------

      mfix_isnan = .False.
      write(notnumber,*) x
! To check for NaN's in x, see if x (a real number) contain a letter "N"
! "n" or symbol "?", in which case it is a NaN (Not a Number)

      IF(INDEX(notnumber,'?') > 0 .OR.     &
         INDEX(notnumber,'n') > 0 .OR.     &
         INDEX(notnumber,'N') > 0 ) THEN
        mfix_isnan = .TRUE.
         RETURN
      ENDIF

      RETURN
    END FUNCTION mfix_isnan

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function: CHECK_VEL_BOUND()                                         C
!  Purpose: Check velocities upper bound to be less than speed of      C
!           sound                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: 25-AUG-05  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      logical function check_vel_bound (slo,shi,ulo,uhi,vlo,vhi,wlo,whi,u_g,v_g,w_g,ep_g)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE toleranc, only: max_inlet_vel

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

      real(c_real)  , intent(in   ) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real)  , intent(in   ) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real)  , intent(in   ) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real)  , intent(in   ) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer :: I,J,K

      check_vel_bound = .FALSE.

      do k = wlo(3),whi(3)
      do j = wlo(2),whi(2)
      do i = wlo(1),whi(1)
         if ( abs(w_g(I,J,K)) > max_inlet_vel) THEN
            check_vel_bound = .true.
            write(*,1000) max_inlet_vel, I, J, K, &
                          ep_g(I,J,K), w_g(I,J,K)
         end if
      end do
      end do
      end do

      do k = vlo(3),vhi(3)
      do j = vlo(2),vhi(2)
      do i = vlo(1),vhi(1)
         if (abs(v_g(I,J,K)) > max_inlet_vel) then
            check_vel_bound = .true.
            write(*,1000) max_inlet_vel, I, J, K, &
                          ep_g(I,J,K), v_g(I,J,K)
         end if
      end do
      end do
      end do

      do k = ulo(3),uhi(3)
      do j = ulo(2),uhi(2)
      do i = ulo(1),uhi(1)
         if (abs(u_g(I,J,K)) > max_inlet_vel) then
            check_vel_bound = .true.
            write(*,1000) max_inlet_vel, I, J, K, &
                          ep_g(I,J,K), u_g(I,J,K)
         end if
      end do
      end do
      end do

 1000 FORMAT(1X,'Message from: CHECK_VEL_BOUND',/&
            'WARNING: w velocity higher than maximum allowed velocity: ', &
            G12.5, '(to change this adjust the scale factor max_inlet_vel_FAC)'/&
            'in this cell: ','I = ',I4,2X,' J = ',I4,2X,' K = ',I4, /&
            '  ','Epg = ', G12.5, 'Wg = ', G12.5)
 1001 FORMAT(1X,'Message from: CHECK_VEL_BOUND',/&
            'WARNING: v velocity higher than maximum allowed velocity: ', &
            G12.5, '(to change this adjust the scale factor max_inlet_vel_FAC)'/&
            'in this cell: ','I = ',I4,2X,' J = ',I4,2X,' K = ',I4, /&
            '  ','Epg = ', G12.5, 'Vg = ', G12.5)
 1002 FORMAT(1X,'Message from: CHECK_VEL_BOUND',/&
            'WARNING: u velocity higher than maximum allowed velocity: ', &
            G12.5, '(to change this adjust the scale factor max_inlet_vel_FAC)'/&
            'in this cell: ','I = ',I4,2X,' J = ',I4,2X,' K = ',I4, /&
            '  ','Epg = ', G12.5, 'Ug = ', G12.5)

      end function check_vel_bound

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function name: SEEK_COMMENT (LINE_MAXCOL)                           !
!  Author: P.Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
!  Purpose: Returns the index to where a comment character was found   !
!  in the input data line.  Equals MAXCOL + 1 if no-comment characters !
!  in the line.                                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      PURE INTEGER FUNCTION SEEK_COMMENT (LINE, MAXCOL)

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Input data line
      CHARACTER(len=*), intent(IN) :: LINE
! Maximum column of input data line to search
      INTEGER, intent(IN) :: MAXCOL

! Local Variables
!---------------------------------------------------------------------//
! The number of designated comment characters
      INTEGER, PARAMETER :: DIM_COMMENT = 2
! The comment characters
      CHARACTER, PARAMETER :: COMMENT_CHAR(DIM_COMMENT) = (/'#', '!'/)
! Loop indicies
      INTEGER :: L, L2
!.......................................................................!

      DO L = 1, MAXCOL
         DO L2 = 1, DIM_COMMENT
            IF (LINE(L:L) == COMMENT_CHAR(L2)) THEN
               SEEK_COMMENT = L
               RETURN
            ENDIF
         END DO
      END DO
      SEEK_COMMENT = MAXCOL + 1
!
      RETURN
      END FUNCTION SEEK_COMMENT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function name: SEEK_END (LINE, MAXCOL)                              !
!  Author: P.Nicoletti, M. Syamlal                    Date: 7-AUG-92   !
!                                                                      !
!  Purpose: Return the index to where the last character was found in  !
!  the input data line.  Equals MAXCOL if no trailing blank characters !
!  in the line.                                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      PURE INTEGER FUNCTION SEEK_END (LINE, MAXCOL)

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! input data line
      CHARACTER, intent(IN) :: LINE*(*)
! maximum column of input data line to search
      INTEGER, intent(IN) :: MAXCOL

! Local Variables
!---------------------------------------------------------------------//
      INTEGER :: L
!.......................................................................!

      SEEK_END = 0
      DO L = 1, MAXCOL
         IF (LINE(L:L) /= ' ') SEEK_END = L
      END DO
      RETURN
      END FUNCTION SEEK_END

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Function name: LINE_TOO_BIG (LINE,LINE_LEN,MAXCOL)                  !
!  Author: P.Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
!  Purpose: Return a value greater than 0 to indicate an error         !
!  condition (data passed column MAXCOL in LINE)                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      PURE INTEGER FUNCTION LINE_TOO_BIG (LINE, LINE_LEN, MAXCOL)

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! input data line
      CHARACTER(LEN=*), intent(IN) :: LINE
! length of input data line
      INTEGER, intent(IN) :: LINE_LEN
! maximum column that non-blank charcater are in the input data line
      INTEGER, intent(IN) :: MAXCOL

! Local Variables
!---------------------------------------------------------------------//
      INTEGER :: L
!.......................................................................!

      DO L = MAXCOL + 1, LINE_LEN
         IF (LINE(L:L) /= ' ') THEN
            LINE_TOO_BIG = L
            RETURN
         ENDIF
      END DO
      LINE_TOO_BIG = 0
      RETURN
      END FUNCTION LINE_TOO_BIG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Function: BLANK_LINE                                                !
! Author: P. Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
! Purpose: Return .TRUE. if a line contains no input or only spaces.   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      PURE LOGICAL FUNCTION BLANK_LINE (line)

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
      CHARACTER, intent(IN) :: LINE*(*)

! Local Variables
!---------------------------------------------------------------------//
      INTEGER :: L
!.......................................................................!

      BLANK_LINE = .FALSE.
      DO L=1, len(line)
         IF(line(L:L)/=' ' .and. line(L:L)/='    ')RETURN
      ENDDO

      BLANK_LINE = .TRUE.
      RETURN
      END FUNCTION BLANK_LINE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MAKE_UPPER_CASE (LINE_STRING,MAXCOL)                   C
!  Author: P.Nicoletti                                Date: 26-NOV-91  C
!                                                                      C
!  Purpose: change lowercase characters to uppercase in input line     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE MAKE_UPPER_CASE(LINE_STRING, MAXCOL)

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Input line to change to uppercase
      CHARACTER(len=*), intent(INOUT) :: LINE_STRING
! Number of characters to look at in LINE_STRING
      INTEGER, intent(IN) :: MAXCOL

! Local Variables
!---------------------------------------------------------------------//
! ICHAR value for UPPERCASE A, lowercase a, lowercase z
      INTEGER, PARAMETER :: A_UP = ICHAR('A')
      INTEGER, PARAMETER :: A_LO = ICHAR('a')
      INTEGER, PARAMETER :: Z_LO = ICHAR('z')
! ICHAR differnce between lower and uppercase letters
      INTEGER, PARAMETER :: A_DIFF = A_LO - A_UP
! Holds ICHAR value of current character
      INTEGER :: INT_C
! loop index
      INTEGER :: L
!.......................................................................!

      DO L = 1, MAXCOL
         INT_C = ICHAR(LINE_STRING(L:L))
         IF (A_LO<=INT_C .AND. INT_C<=Z_LO) THEN
            INT_C = INT_C - A_DIFF
            LINE_STRING(L:L) = CHAR(INT_C)
         ENDIF
      END DO
      RETURN
      END SUBROUTINE MAKE_UPPER_CASE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: REPLACE_TAB (LINE_STRING,MAXCOL)                       !
!  Author: M. Syamlal                                 Date: 10-JUL-03  !
!                                                                      !
!  Purpose: replace tab characters with space                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE REPLACE_TAB(LINE_STRING, MAXCOL)

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Input line to change to uppercase
      CHARACTER(len=*), intent(INOUT) :: LINE_STRING
! Number of characters to look at in LINE_STRING
      INTEGER, intent(IN) :: MAXCOL

! Local Variables
!---------------------------------------------------------------------//
      CHARACTER, PARAMETER :: TAB = CHAR(9)
      CHARACTER, PARAMETER :: CRET = CHAR(13)
! Loop index
      INTEGER :: L
!.......................................................................!

      DO L = 1, MAXCOL
        if(LINE_STRING(L:L) .eq. TAB) LINE_STRING(L:L) = ' '
        if(LINE_STRING(L:L) .eq. CRET) LINE_STRING(L:L) = ' '
      END DO
      RETURN
      END SUBROUTINE REPLACE_TAB

END MODULE utilities
