! This file is included by functions_mod.f

! These function are are defined in this file so it can be included by
! individual subprograms in order to inline the (short) functions
! defined in this file.

! All functions that may benefits from inlining should be defined in this file.

! For code that is not performance-critical, do:

! SUBROUTINE subname
! USE functions
! ...
! END SUBROUTINE subname

! For code that is performance-critical, do:

! SUBROUTINE subname
! ...
! CONTAINS
! INCLUDE 'functions.inc'
! END SUBROUTINE subname

!---------------------------------------------------------------------//
! Functions for generating the LOCAL 3-D array index IJK from the
! 1-D indices I, J, and K.
!//FUNIJK is moved to compar for debugging purposes - Sreekanth-10/26/99
!     FUNIJK (LI, LJ, LK) = c0 + LI + (LJ-jstart3_all(myPE))*c1 + (LK-kstart3_all(myPE))* c2
!      funijk(li,lj,lk) = lj + c0 + li*c1 + lk*c2
      INTEGER FUNCTION funijk(li,lj,lk)
      USE compar
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LI, LJ, LK
      funijk = IJK_ARRAY_OF(li,lj,lk)
      END FUNCTION funijk

      INTEGER FUNCTION funijk_0(li,lj,lk)
      USE compar
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LI, LJ, LK
      funijk_0 = lj + c0 + li*c1 + lk*c2
      END FUNCTION funijk_0

! Function for generating the LOCAL 3-D array index IJK from the
! the 1-D indices I, J, K and IPROC.
!     FUNIJK_PROC(LI, LJ, LK, LIPROC) = 1 + (LI - istart3_all(LIPROC))+ &
!     (LJ-jstart3_all(LIPROC))*(iend3_all(LIPROC)-istart3_all(LIPROC)+1) &
!     + (LK-kstart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1)* &
!     (iend3_all(LIPROC)-istart3_all(LIPROC)+1)
      INTEGER FUNCTION FUNIJK_PROC(LI, LJ, LK, LIPROC)
      USE compar
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LI, LJ, LK, LIPROC
      FUNIJK_PROC = 1 + (LJ - jstart3_all(LIPROC))+ &
         (LI-Istart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1) &
         + (LK-kstart3_all(LIPROC))*(jend3_all(LIPROC)-jstart3_all(LIPROC)+1)* &
         (iend3_all(LIPROC)-istart3_all(LIPROC)+1)
      END FUNCTION FUNIJK_PROC

! Function for generating the GLOBAL 3-D array index IJK from the
! 1-D indices I, J, and K.
!     FUNIJK_GL (LI, LJ, LK) = 1 + (LI - imin3) + (LJ-jmin3)*(imax3-imin3+1) &
!     + (LK-kmin3)*(jmax3-jmin3+1)*(imax3-imin3+1)
      INTEGER FUNCTION FUNIJK_GL (LI, LJ, LK)
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LI, LJ, LK
      FUNIJK_GL = 1 + (LJ - jmin3) + (LI-imin3)*(jmax3-jmin3+1) &
         + (LK-kmin3)*(jmax3-jmin3+1)*(imax3-imin3+1)
      END FUNCTION FUNIJK_GL

! Function for generating the 3-D array index IJK from the 1-D indices
! I, J, and K in IO format
      INTEGER FUNCTION FUNIJK_IO (LI, LJ, LK)
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LI, LJ, LK
      FUNIJK_IO = 1 + (LI - imin2) + (LJ-jmin2)*(imax2-imin2+1) &
         + (LK-kmin2)*(jmax2-jmin2+1)*(imax2-imin2+1)
      END FUNCTION FUNIJK_IO

!----------------------------------------------------------------------!
!  Function: IS_ON_myPE_OWNS                                           !
!                                                                      !
!  Purpose: Returns TRUE if the I,J,K values point to a computational  !
!  cell that is OWNED by the current process.                          !
!                                                                      !
!  o Ownership is defined as belonging to the current PE's domain but  !
!    as a cell in any of the PE's ghost layers.                        !
!                                                                      !
!  o Each computational cell is owned by one -and only one- PE.        !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION IS_ON_myPE_OWNS(LI, LJ, LK)
        USE compar
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: LI, LJ, LK

      IS_ON_MYPE_OWNS = &
         LI >= ISTART .AND. LI <= IEND .AND. &
         LJ >= JSTART .AND. LJ <= JEND .AND. &
         LK >= KSTART .AND. LK <= KEND

      RETURN
      END FUNCTION IS_ON_MYPE_OWNS

!----------------------------------------------------------------------!
!  Function: IS_ON_myPE_WOBND                                          !
!                                                                      !
!  Purpose: Returns TRUE if the I,J,K values point to a computational  !
!  cell that is OWNED by the current process and not a exterior ghost  !
!  cell.                                                               !
!                                                                      !
!  o This is a subset of IS_ON_myPE_OWNS.                              !
!                                                                      !
!  o Exterior ghost cells are those in cells surrounding the domain.   !
!    These are cells created to fully define boundary conditions       !
!    (e.g., I == 1 where X_E(1) == ZERO).                              !
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION IS_ON_myPE_wobnd (LI, LJ, LK)
        USE compar
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: LI, LJ, LK

      IS_ON_MYPE_WOBND = &
         LI >= ISTART1 .AND. LI <= IEND1 .AND. &
         LJ >= JSTART1 .AND. LJ <= JEND1 .AND. &
         LK >= KSTART1 .AND. LK <= KEND1

      RETURN
      END FUNCTION IS_ON_myPE_wobnd

!----------------------------------------------------------------------!
!  Function: IS_ON_myPE_Plus1Layer                                     !
!                                                                      !
!  Purpose: Returns TRUE if the I,J,K values point to a computational  !
!  cell that is OWNED by the current process or contained in the fisrt !
!  layer of ghost cells seen by the current PE.                        !
!                                                                      !
!  o This is a superset of IS_ON_myPE_OWNS.                            !
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION IS_ON_myPE_plus1layer (LI, LJ, LK)
        USE compar
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: LI, LJ, LK

      IS_ON_MYPE_PLUS1LAYER = &
         LI >= ISTART2 .AND. LI <= IEND2 .AND. &
         LJ >= JSTART2 .AND. LJ <= JEND2 .AND. &
         LK >= KSTART2 .AND. LK <= KEND2

      RETURN
      END FUNCTION IS_ON_myPE_plus1layer

!----------------------------------------------------------------------!
!  Function: IS_ON_myPE_Plus2Layer                                     !
!                                                                      !
!  Purpose: Returns TRUE if the I,J,K values point to a computational  !
!  cell that is OWNED by the current process or contained in the fisrt !
!  two layers of ghost cells seen by the current PE.                   !
!                                                                      !
!  o This is a superset of IS_ON_Plus1Layer.                           !
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION IS_ON_myPE_plus2layers (LI, LJ, LK)
        USE compar
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: LI, LJ, LK

      IS_ON_MYPE_PLUS2LAYERS = &
         LI >= ISTART3 .AND. LI <= IEND3 .AND. &
         LJ >= JSTART3 .AND. LJ <= JEND3 .AND. &
         LK >= KSTART3 .AND. LK <= KEND3

      RETURN
      END FUNCTION IS_ON_myPE_plus2layers

!---------------------------------------------------------------------//
! logical function to identify a fluid cell
      LOGICAL FUNCTION FLUID_AT(IJK)
      USE geometry, only: flag
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      FLUID_AT    = FLAG(IJK) .EQ. 1
      END FUNCTION FLUID_AT

      logical function fluid_cell(i,j,k)
         use geometry, only: flag
         implicit none
         integer, intent(in) :: i,j,k
         fluid_cell = (flag(funijk(i,j,k)).eq.1)
      end function fluid_cell

! logical function to identify a specified pressure inflow cell
      LOGICAL FUNCTION P_FLOW_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      P_FLOW_AT = FLAG(IJK) .EQ. 10 .OR. &
         FLAG(IJK) .EQ. 11
      END FUNCTION P_FLOW_AT

! logical function to identify a specified pressure outflow cell
      LOGICAL FUNCTION P_OUTFLOW_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      P_OUTFLOW_AT= FLAG(IJK) .EQ. 11
      END FUNCTION P_OUTFLOW_AT

! logical function to identify either a specified pressure inflow
! or outflow cell or a fluid cell (simplified check)
! FLUID_AT or P_FLOW_AT (simplified check)
      logical function FLUIDorP_FLOW_AT(i,j,k)
      USE geometry
      implicit none
      integer, intent(in) :: i,j,k
      FLUIDorP_FLOW_AT = (flag(FUNIJK(i,j,k)) .LE. 11)
      end function FLUIDorP_FLOW_AT

! logical function to identify a specified mass outflow cell
      LOGICAL FUNCTION MASS_OUTFLOW_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      MASS_OUTFLOW_AT= FLAG(IJK) .EQ. 21
      END FUNCTION MASS_OUTFLOW_AT

! logical function to identify a specified outflow cell
      LOGICAL FUNCTION OUTFLOW_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      OUTFLOW_AT  = FLAG(IJK) .EQ. 31
      END FUNCTION OUTFLOW_AT

! logical function to identify any type of flow in/out at cell
! pressure inflow/outflow, mass inflow/outflow or outflow
      LOGICAL FUNCTION FLOW_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      FLOW_AT     = FLAG(IJK) .GE. 10 .AND. FLAG(IJK) .LE. 31
      END FUNCTION FLOW_AT

! Logical function to identify default walls
      LOGICAL FUNCTION WALL_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      WALL_AT     = FLAG(IJK) .GE. 100
      END FUNCTION WALL_AT

      logical function wall_cell(i,j,k)
         use geometry, only: flag
         implicit none
         integer, intent(in) :: i,j,k
         wall_cell = (flag(funijk(i,j,k)).ge.100)
      end function wall_cell

! Logical function to identify a No-slip wall cell
      LOGICAL FUNCTION NS_WALL_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      NS_WALL_AT  = FLAG(IJK) .EQ. 100
      END FUNCTION NS_WALL_AT

      ! Logical function to identify a no-slip wall cell
      logical function ns_wall_cell(i,j,k)
         use geometry, only: flag
         implicit none
         integer, intent(in) :: i,j,k
         ns_wall_cell = (flag(funijk(i,j,k)).eq.100)
      end function ns_wall_cell

! Logical function to identify a Free-slip wall cell
      LOGICAL FUNCTION FS_WALL_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      FS_WALL_AT  = FLAG(IJK) .EQ. 101
      END FUNCTION FS_WALL_AT

      logical function fs_wall_cell(i,j,k)
         use geometry, only: flag
         implicit none
         integer, intent(in) :: i,j,k
         fs_wall_cell = (flag(funijk(i,j,k)).eq.101)
      end function fs_wall_cell

! Logical function to identify a Partial-slip wall cell
      LOGICAL FUNCTION PS_WALL_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      PS_WALL_AT  = FLAG(IJK) .EQ. 102
      END FUNCTION PS_WALL_AT

! Logical function to identify wall ICBC_FLAG
      LOGICAL FUNCTION WALL_ICBC_FLAG(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      WALL_ICBC_FLAG = ICBC_FLAG(IJK)(1:1) .EQ. 'W' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 'S' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 's' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 'c' .OR. &
         ICBC_FLAG(IJK)(1:1) .EQ. 'C'
      END FUNCTION WALL_ICBC_FLAG

      LOGICAL FUNCTION DEFAULT_WALL_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      DEFAULT_WALL_AT = ICBC_FLAG(IJK)(2:3) .EQ. '--' .AND. &
         (ICBC_FLAG(IJK)(1:1) .NE. 'c'  .AND. &
         ICBC_FLAG(IJK)(1:1) .NE. 'C')
      END FUNCTION DEFAULT_WALL_AT


! Cyclic
!---------------------------------------------------------------------//
      LOGICAL FUNCTION CYCLIC_AT(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      CYCLIC_AT = FLAG(IJK) .EQ. 106 .OR. &
         FLAG(IJK) .EQ. 107
      END FUNCTION CYCLIC_AT

! logical function to identify cyclic condition at east boundary
      LOGICAL FUNCTION CYCLIC_AT_E(IJK)
      USE geometry, only: flag_e
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      CYCLIC_AT_E   = FLAG_E(IJK) .EQ. 2000
      END FUNCTION CYCLIC_AT_E

! logical function to identify cyclic condition at north boundary
      LOGICAL FUNCTION CYCLIC_AT_N(IJK)
      USE geometry, only: flag_n
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      CYCLIC_AT_N   = FLAG_N(IJK) .EQ. 2000
      END FUNCTION CYCLIC_AT_N

! logical function to identify cyclic condition at top boundary
      LOGICAL FUNCTION CYCLIC_AT_T(IJK)
      USE geometry, only: flag_t
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      CYCLIC_AT_T   = FLAG_T(IJK) .EQ. 2000
      END FUNCTION CYCLIC_AT_T


! Flow boundaries
!---------------------------------------------------------------------//
! identify flow at east boundary
      LOGICAL FUNCTION FLOW_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      FLOW_AT_E   = FLAG_E(IJK) .GE. 2000 .AND.&
         FLAG_E(IJK) .LE. 2011
      END FUNCTION FLOW_AT_E

! identify specified flow north boundary
      LOGICAL FUNCTION FLOW_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      FLOW_AT_N   = FLAG_N(IJK) .GE. 2000 .AND.&
         FLAG_N(IJK) .LE. 2011
      END FUNCTION FLOW_AT_N

! identify specified flow top boundary
      LOGICAL FUNCTION FLOW_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJK
      FLOW_AT_T   = FLAG_T(IJK) .GE. 2000 .AND.&
         FLAG_T(IJK) .LE. 2011
      END FUNCTION FLOW_AT_T

! identify const. pressure flow top boundary
      LOGICAL FUNCTION PFLOW_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      PFLOW_AT_E  = FLAG_E(IJK) .EQ. 2010 .OR.&
         FLAG_E(IJK) .EQ. 2011
      END FUNCTION PFLOW_AT_E

! identify const. pressure flow north boundary
      LOGICAL FUNCTION PFLOW_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      PFLOW_AT_N  = FLAG_N(IJK) .EQ. 2010 .OR.&
         FLAG_N(IJK) .EQ. 2011
      END FUNCTION PFLOW_AT_N

! identify const. pressure flow east boundary
      LOGICAL FUNCTION PFLOW_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      PFLOW_AT_T  = FLAG_T(IJK) .EQ. 2010 .OR.&
         FLAG_T(IJK) .EQ. 2011
      END FUNCTION PFLOW_AT_T

! identify specified flow east boundary
      LOGICAL FUNCTION MFLOW_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      MFLOW_AT_E  = FLAG_E(IJK) .EQ. 2020 .OR. &
         FLAG_E(IJK) .EQ. 2021 .OR. &
         FLAG_E(IJK) .EQ. 2031
      END FUNCTION MFLOW_AT_E

! identify specified flow north boundary
      LOGICAL FUNCTION MFLOW_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      MFLOW_AT_N  = FLAG_N(IJK) .EQ. 2020 .OR. &
         FLAG_N(IJK) .EQ. 2021 .OR. &
         FLAG_N(IJK) .EQ. 2031
      END FUNCTION MFLOW_AT_N

! identify specified flow top boundary
      LOGICAL FUNCTION MFLOW_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      MFLOW_AT_T  = FLAG_T(IJK) .EQ. 2020 .OR. &
         FLAG_T(IJK) .EQ. 2021 .OR. &
         FLAG_T(IJK) .EQ. 2031
      END FUNCTION MFLOW_AT_T


! Functions to identify a impermeable and/or semi-permeable surface at
! indicated boundary (specific type of internal surface)
!---------------------------------------------------------------------//
! Logical function to identify IP (impermeable surface) East of the cell
      logical function ip_at_e(i,j,k)
      use geometry
      implicit none
      integer, intent(in) :: i,j,k
      ip_at_e = flag_e(funijk(i,j,k)) .LT. 1000
      end function ip_at_e

! Logical function to identify IP (impermeable surface) North of the cell
      logical function ip_at_n(i,j,k)
      use geometry
      implicit none
      integer, intent(in) :: i,j,k
      ip_at_n = flag_n(funijk(i,j,k)) .LT. 1000
      end function ip_at_n

! Logical function to identify IP (impermeable surface) Top of the cell
      logical function ip_at_t(i,j,k)
      use geometry
      implicit none
      integer, intent(in) :: i,j,k
      ip_at_t = flag_t(funijk(i,j,k)) .LT. 1000
      end function ip_at_t

! Logical function to identify SP or IP (semi or impermeable surface)
! at east of the cell
      LOGICAL FUNCTION SIP_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SIP_AT_E    = (FLAG_E(IJK) .LT. 2000)
      END FUNCTION SIP_AT_E

! Logical function to identify SP or IP (semi or impermeable surface)
! at north of the cell
      LOGICAL FUNCTION SIP_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SIP_AT_N    = (FLAG_N(IJK) .LT. 2000)
      END FUNCTION SIP_AT_N

! Logical function to identify SP or IP (semi or impermeable surface)
! at top of the cell
      LOGICAL FUNCTION SIP_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SIP_AT_T    = (FLAG_T(IJK) .LT. 2000)
      END FUNCTION SIP_AT_T

! Logical function to identify SP (semi-permeable surface) at east
! of cell
      LOGICAL FUNCTION SP_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SP_AT_E     = (FLAG_E(IJK) .LT. 2000) .AND. &
         (FLAG_E(IJK) .GE. 1000)
      END FUNCTION SP_AT_E

! Logical function to identify SP (semi-permeable surface) at north
! of cell
      LOGICAL FUNCTION SP_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SP_AT_N     = (FLAG_N(IJK) .LT. 2000) .AND. &
         (FLAG_N(IJK) .GE. 1000)
      END FUNCTION SP_AT_N

! Logical function to identify SP (semi-permeable surface) at top
! of cell
      LOGICAL FUNCTION SP_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      SP_AT_T     = (FLAG_T(IJK) .LT. 2000) .AND. &
         (FLAG_T(IJK) .GE. 1000)
      END FUNCTION SP_AT_T


! Logical functions concerning general internal surfaces
! Integer functions to return internal surface ID
!---------------------------------------------------------------------//
! Internal surface ID for east face
      INTEGER FUNCTION IS_ID_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_ID_AT_E = FLAG_E(IJK) - 1000
      END FUNCTION IS_ID_AT_E

! Internal surface ID for north face
      INTEGER FUNCTION IS_ID_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_ID_AT_N  = FLAG_N(IJK) - 1000
      END FUNCTION IS_ID_AT_N

! Internal surface ID for top face
      INTEGER FUNCTION IS_ID_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_ID_AT_T  = FLAG_T(IJK) - 1000
      END FUNCTION IS_ID_AT_T

! Logical function to identify IS at East of the cell
      LOGICAL FUNCTION IS_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_AT_E     = FLAG_E(IJK) .LT. 2000
      END FUNCTION IS_AT_E

! Logical function to identify IS at North of the cell
      LOGICAL FUNCTION IS_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_AT_N     = FLAG_N(IJK) .LT. 2000
      END FUNCTION IS_AT_N

! Logical function to identify IS at Top of the cell
      LOGICAL FUNCTION IS_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      IS_AT_T     = FLAG_T(IJK) .LT. 2000
      END FUNCTION IS_AT_T

! Logical function to identify No IS at East of the cell
      LOGICAL FUNCTION NO_IS_AT_E(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      NO_IS_AT_E  = FLAG_E(IJK) .GE. 2000
      END FUNCTION NO_IS_AT_E

! Logical function to identify No IS at North of the cell
      LOGICAL FUNCTION NO_IS_AT_N(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      NO_IS_AT_N  = FLAG_N(IJK) .GE. 2000
      END FUNCTION NO_IS_AT_N

! Logical function to identify No IS at Top of the cell
      LOGICAL FUNCTION NO_IS_AT_T(IJK)
      USE geometry
      IMPLICIT NONE
      INTEGER IJK
      NO_IS_AT_T  = FLAG_T(IJK) .GE. 2000
      END FUNCTION NO_IS_AT_T

! Misc
!---------------------------------------------------------------------//

! Function that returns the maximum of zero or input
      DOUBLE PRECISION FUNCTION ZMAX(XXX)
      USE param1, only: zero
      IMPLICIT NONE
      DOUBLE PRECISION XXX
      ZMAX       = MAX(XXX, ZERO)
      END FUNCTION ZMAX

      LOGICAL FUNCTION IS_NONEXISTENT(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NONEXISTENT
        INTEGER, INTENT(IN) :: PP
        IS_NONEXISTENT = (PARTICLE_STATE(PP)==NONEXISTENT)
      END FUNCTION IS_NONEXISTENT

      LOGICAL FUNCTION IS_NORMAL(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NORMAL_PARTICLE
        INTEGER, INTENT(IN) :: PP
        IS_NORMAL = (PARTICLE_STATE(PP)==NORMAL_PARTICLE)
      END FUNCTION IS_NORMAL

      LOGICAL FUNCTION IS_ENTERING(PP)
        USE discretelement, ONLY: PARTICLE_STATE, ENTERING_PARTICLE
        INTEGER, INTENT(IN) :: PP
        IS_ENTERING = (PARTICLE_STATE(PP)==ENTERING_PARTICLE)
      END FUNCTION IS_ENTERING

      LOGICAL FUNCTION IS_EXITING(PP)
        USE discretelement, ONLY: PARTICLE_STATE, EXITING_PARTICLE
        INTEGER, INTENT(IN) :: PP
        IS_EXITING = (PARTICLE_STATE(PP)==EXITING_PARTICLE)
      END FUNCTION IS_EXITING

      LOGICAL FUNCTION IS_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NORMAL_GHOST
        INTEGER, INTENT(IN) :: PP
        IS_GHOST = (PARTICLE_STATE(PP)==NORMAL_GHOST)
      END FUNCTION IS_GHOST

      LOGICAL FUNCTION IS_ENTERING_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, ENTERING_GHOST
        INTEGER, INTENT(IN) :: PP
        IS_ENTERING_GHOST = (PARTICLE_STATE(PP)==ENTERING_GHOST)
      END FUNCTION IS_ENTERING_GHOST

      LOGICAL FUNCTION IS_EXITING_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, EXITING_GHOST
        INTEGER, INTENT(IN) :: PP
        IS_EXITING_GHOST = (PARTICLE_STATE(PP)==EXITING_GHOST)
      END FUNCTION IS_EXITING_GHOST

      LOGICAL FUNCTION IS_ANY_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NORMAL_GHOST
        USE discretelement, ONLY: ENTERING_GHOST, EXITING_GHOST
        INTEGER, INTENT(IN) :: PP
        IS_ANY_GHOST = ((PARTICLE_STATE(PP)==NORMAL_GHOST) .OR.        &
           (PARTICLE_STATE(PP)==ENTERING_GHOST) .OR.                   &
           (PARTICLE_STATE(PP)==EXITING_GHOST))
      END FUNCTION IS_ANY_GHOST

      SUBROUTINE SET_NONEXISTENT(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NONEXISTENT
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=NONEXISTENT
      END SUBROUTINE SET_NONEXISTENT

      SUBROUTINE SET_NORMAL(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NORMAL_PARTICLE
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=NORMAL_PARTICLE
      END SUBROUTINE SET_NORMAL

      SUBROUTINE SET_ENTERING(PP)
        USE discretelement, ONLY: PARTICLE_STATE, ENTERING_PARTICLE
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=ENTERING_PARTICLE
      END SUBROUTINE SET_ENTERING

      SUBROUTINE SET_EXITING(PP)
        USE discretelement, ONLY: PARTICLE_STATE, EXITING_PARTICLE
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=EXITING_PARTICLE
      END SUBROUTINE SET_EXITING

      SUBROUTINE SET_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, NORMAL_GHOST
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=NORMAL_GHOST
      END SUBROUTINE SET_GHOST

      SUBROUTINE SET_ENTERING_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, ENTERING_GHOST
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=ENTERING_GHOST
      END SUBROUTINE SET_ENTERING_GHOST

      SUBROUTINE SET_EXITING_GHOST(PP)
        USE discretelement, ONLY: PARTICLE_STATE, EXITING_GHOST
        INTEGER, INTENT(IN) :: PP
        PARTICLE_STATE(PP)=EXITING_GHOST
      END SUBROUTINE SET_EXITING_GHOST

      ! ************************************************

      integer function iminus(i,j,k)
      use compar, only: istart3
      use geometry, only  : imin1,imax1,cyclic_x
      implicit none
      integer, intent(in) :: i,j,k
      integer             :: ijk

      ijk = funijk(i,j,k)

      if (i.eq.imin1 .and. (wall_at(ijk) .or. flow_at(ijk)) ) then
         iminus = i
      else
         iminus = max(istart3, i-1)
      end if

      if (cyclic_x .and. (iminus < imin1) .and. .not.wall_at(ijk)) iminus = imax1

      end function iminus

      ! ************************************************

      integer function iplus(i,j,k)
      use compar, only: iend3
      use geometry, only  : imin1,imax1,cyclic_x
      implicit none
      integer, intent(in) :: i,j,k
      integer             :: ijk

      ijk = funijk(i,j,k)

      if (i.eq.imax1 .and. (wall_at(ijk) .or. flow_at(ijk)) ) then
         iplus = i
      else
         iplus = min(iend3, i+1)
      end if

      if (cyclic_x .and. (iplus > imax1) .and. .not.wall_at(ijk)) iplus = imin1

      end function iplus

      ! ************************************************

      integer function jminus(i,j,k)
      use compar, only: jstart3
      use geometry, only  : jmin1,jmax1,cyclic_y
      implicit none
      integer, intent(in) :: i,j,k
      integer             :: ijk

      ijk = funijk(i,j,k)

      if (j.eq.jmin1 .and. (wall_at(ijk) .or. flow_at(ijk)) ) then
         jminus = j
      else
         jminus = max(jstart3, j-1)
      end if

      if (cyclic_y .and. (jminus < jmin1) .and. .not.wall_at(ijk)) jminus = jmax1

      end function jminus

      ! ************************************************

      integer function jplus(i,j,k)
      use compar, only: jend3
      use geometry, only  : jmin1,jmax1,cyclic_y
      implicit none
      integer, intent(in) :: i,j,k
      integer             :: ijk

      ijk = funijk(i,j,k)

      if (j.eq.jmax1 .and. (wall_at(ijk) .or. flow_at(ijk)) ) then
         jplus = j
      else
         jplus = min(jend3, j+1)
      end if

      if (cyclic_y .and. (jplus > jmax1) .and. .not.wall_at(ijk)) jplus = jmin1

      end function jplus

      ! ************************************************

      integer function kminus(i,j,k)

      use compar, only: kstart3
      use geometry, only  : kmin1,kmax1,cyclic_z
      implicit none
      integer, intent(in) :: i,j,k
      integer             :: ijk

      ijk = funijk(i,j,k)

      if (k.eq.kmin1 .and. (wall_at(ijk) .or. flow_at(ijk)) ) then
         kminus = k
      else
         kminus = max(kstart3, k-1)
      end if

      if (cyclic_z .and. (kminus < kmin1) .and. .not.wall_at(ijk)) kminus = kmax1

      end function kminus

      ! ************************************************

      integer function kplus(i,j,k)

      use compar, only: kend3
      use geometry, only  : kmin1,kmax1,cyclic_z
      implicit none
      integer, intent(in) :: i,j,k
      integer             :: ijk

      ijk = funijk(i,j,k)

      if (k.eq.kmax1 .and. (wall_at(ijk) .or. flow_at(ijk)) ) then
         kplus = k
      else
         kplus = min(kend3, k+1)
      end if

      if (cyclic_z .and. (kplus > kmax1) .and. .not.wall_at(ijk)) kplus = kmin1

      end function kplus

      ! ************************************************

      integer function ieast(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k
      integer             :: ipjk

      ipjk = funijk(iplus(i,j,k),j,k)

      if (wall_at(ipjk)) then
         ieast = i
      else
         ieast = iplus(i,j,k)
      end if

      end function ieast

      ! ************************************************

      integer function iwest(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k
      integer             :: imjk

      imjk = funijk(iminus(i,j,k),j,k)

      if (wall_at(imjk)) then
         iwest = i
      else
         iwest = iminus(i,j,k)
      end if

      end function iwest

      ! ************************************************

      integer function jnorth(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k
      integer             :: ijpk

      ijpk = funijk(i,jplus(i,j,k),k)

      if (wall_at(ijpk)) then
         jnorth = j
      else
         jnorth = jplus(i,j,k)
      end if

      end function jnorth

      ! ************************************************

      integer function jsouth(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k
      integer             :: ijmk

      ijmk = funijk(i,jminus(i,j,k),k)

      if (wall_at(ijmk)) then
         jsouth = j
      else
         jsouth = jminus(i,j,k)
      end if

      end function jsouth

      ! ************************************************

      integer function ktop(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k
      integer             :: kp,ijkp

      kp = kplus(i,j,k)
      ijkp = funijk(i,j,kp)

      if (wall_at(ijkp)) then
         ktop = k
      else
         ktop = kp
      end if

      end function ktop

      ! ************************************************

      integer function kbot(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k
      integer             :: km,ijkm

      km = kminus(i,j,k)
      ijkm = funijk(i,j,km)

      if (wall_at(ijkm)) then
         kbot = k
      else
         kbot = km
      end if

      end function kbot
