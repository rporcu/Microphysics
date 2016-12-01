MODULE utilities

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
      double precision x
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      CHARACTER(LEN=80) :: notnumber
!-----------------------------------------------

      mfix_isnan = .False.
      WRITE(notnumber,*) x
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
!  function: MAX_VEL_INLET                                             C
!  Purpose: Find maximum velocity at inlets.                           C
!                                                                      C
!  Author: S. Benyahia                                Date: 26-AUG-05  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      DOUBLE PRECISION FUNCTION MAX_VEL_INLET(u_g,v_g,w_g)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_defined, bc_type
      use bc, only: bc_plane
      use bc, only: bc_k_b, bc_k_t, bc_j_s, bc_j_n, bc_i_w, bc_i_e
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      use functions, only: iminus, jminus, kminus
      use param    , only: dimension_bc
      use param1   , only: zero, small_number
      use run      , only: units
      use toleranc , only: max_allowed_vel, max_inlet_vel, max_inlet_vel_fac

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: L, I, J, K
      DOUBLE PRECISION :: maxVEL
!---------------------------------------------------------------------//

      maxVEL = ZERO

      DO L = 1, DIMENSION_BC

         IF(.NOT.(BC_DEFINED(L).AND.(BC_TYPE(L)=='MASS_INFLOW' .OR. &
            BC_TYPE(L) == 'P_INFLOW'))) CYCLE

         DO K = BC_K_B(L), BC_K_T(L)
         DO J = BC_J_S(L), BC_J_N(L)
         DO I = BC_I_W(L), BC_I_E(L)

            SELECT CASE (BC_PLANE(L))
            CASE ('E'); maxVEL = max(maxVEL,abs(U_G(I,J,K)))
            CASE ('N'); maxVEL = max(maxVEL,abs(V_G(I,J,K)))
            CASE ('T'); maxVEL = max(maxVEL,abs(W_G(I,J,K)))
            CASE ('W'); maxVEL = max(maxVEL,abs(U_G(iminus(i,j,k),j,k)))
            CASE ('S'); maxVEL = max(maxVEL,abs(V_G(i,jminus(i,j,k),k)))
            CASE ('B'); maxVEL = max(maxVEL,abs(W_G(i,j,kminus(i,j,k))))
            END SELECT

         ENDDO
         ENDDO
         ENDDO
      ENDDO

      ! CALL GLOBAL_ALL_MAX(maxVEL, MAX_VEL_INLET)

! If no inlet velocity is specified, use an upper limit defined in
! toleranc_mod.f
      IF(MAX_VEL_INLET <= SMALL_NUMBER) THEN
         MAX_VEL_INLET = MAX_ALLOWED_VEL
         IF(UNITS == 'SI') MAX_INLET_VEL = 1D-2*MAX_ALLOWED_VEL
      ELSE
! Scale the value using a user defined scale factor
         MAX_VEL_INLET = 100.0d0*MAX_INLET_VEL_FAC*MAX_VEL_INLET
      ENDIF

      RETURN
      END FUNCTION MAX_VEL_INLET


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

      LOGICAL FUNCTION CHECK_VEL_BOUND (u_g,v_g,w_g,ep_g)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE toleranc, only: max_inlet_vel
      USE functions, only: fluid_at

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I,J,K
      LOGICAL :: ALL_IS_ERROR
!-----------------------------------------------

! initializing
      CHECK_VEL_BOUND = .FALSE.
      ALL_IS_ERROR    = .FALSE.

LOOP_FLUID: DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IF (fluid_at(i,j,k)) THEN
            IF(ABS(U_G(I,J,K)) > MAX_INLET_VEL .OR. &
               ABS(V_G(I,J,K)) > MAX_INLET_VEL .OR. &
               ABS(W_G(I,J,K)) > MAX_INLET_VEL) THEN
               CHECK_VEL_BOUND = .TRUE.
               WRITE(*,1000) MAX_INLET_VEL, I, J, K, &
                             EP_g(I,J,K), U_G(I,J,K), V_G(I,J,K), W_G(I,J,K)
               EXIT LOOP_FLUID
            ENDIF
         ENDIF

      ENDDO
      ENDDO
      ENDDO LOOP_FLUID

      ! CALL GLOBAL_ALL_OR(CHECK_VEL_BOUND, ALL_IS_ERROR)
      IF(ALL_IS_ERROR) CHECK_VEL_BOUND = .TRUE.

      RETURN
 1000 FORMAT(1X,'Message from: CHECK_VEL_BOUND',/&
            'WARNING: velocity higher than maximum allowed velocity: ', &
            G12.5, '(to change this adjust the scale factor MAX_INLET_VEL_FAC)'/&
            'in this cell: ','I = ',I4,2X,' J = ',I4,2X,' K = ',I4, /&
            '  ','Epg = ', G12.5, 'Ug = ', G12.5, 'Vg = ', G12.5, 'Wg = ', G12.5)

      END FUNCTION CHECK_VEL_BOUND

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
      CHARACTER(len=*), INTENT(IN) :: LINE
! Maximum column of input data line to search
      INTEGER, INTENT(IN) :: MAXCOL

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
      CHARACTER, INTENT(IN) :: LINE*(*)
! maximum column of input data line to search
      INTEGER, INTENT(IN) :: MAXCOL

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
      CHARACTER(LEN=*), INTENT(IN) :: LINE
! length of input data line
      INTEGER, INTENT(IN) :: LINE_LEN
! maximum column that non-blank charcater are in the input data line
      INTEGER, INTENT(IN) :: MAXCOL

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
      CHARACTER, INTENT(IN) :: LINE*(*)

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
      CHARACTER(len=*), INTENT(INOUT) :: LINE_STRING
! Number of characters to look at in LINE_STRING
      INTEGER, INTENT(IN) :: MAXCOL

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
      CHARACTER(len=*), INTENT(INOUT) :: LINE_STRING
! Number of characters to look at in LINE_STRING
      INTEGER, INTENT(IN) :: MAXCOL

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
