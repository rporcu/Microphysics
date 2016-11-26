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
      INTEGER FUNCTION funijk(li,lj,lk)
      USE compar, only: c0, c1, c2
      implicit none
      INTEGER, INTENT(IN) :: LI, LJ, LK
      funijk = c0 + li + lj*c1 + lk*c2
      END FUNCTION funijk

! Function for generating the LOCAL 3-D array index IJK from the
! the 1-D indices I, J, K and IPROC.
      INTEGER FUNCTION FUNIJK_PROC(LI, LJ, LK, lproc)
      USE compar
      implicit none
      INTEGER, INTENT(IN) :: LI, LJ, LK, lproc
      FUNIJK_PROC = 1 + (li - istart3_all(lproc))+ &
         (lj-jstart3_all(lproc))*(iend3_all(lproc)-istart3_all(lproc)+1) &
         + (LK-kstart3_all(lproc))*(iend3_all(lproc)-istart3_all(lproc)+1)* &
         (jend3_all(lproc)-jstart3_all(lproc)+1)
      END FUNCTION FUNIJK_PROC

! Function for generating the GLOBAL 3-D array index IJK from the
! 1-D indices I, J, and K.
      INTEGER FUNCTION FUNIJK_GL (LI, LJ, LK)
      USE geometry
      implicit none
      INTEGER, INTENT(IN) :: LI, LJ, LK
      FUNIJK_GL = 1 + (LJ - jmin3) + (LI-imin3)*(jmax3-jmin3+1) &
         + (LK-kmin3)*(jmax3-jmin3+1)*(imax3-imin3+1)
      END FUNCTION FUNIJK_GL

! Function for generating the 3-D array index IJK from the 1-D indices
! I, J, and K in IO format
      INTEGER FUNCTION FUNIJK_IO (LI, LJ, LK)
      USE geometry
      implicit none
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
        implicit none

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
        implicit none

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
        implicit none

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
        implicit none

        INTEGER, INTENT(IN) :: LI, LJ, LK

      IS_ON_MYPE_PLUS2LAYERS = &
         LI >= ISTART3 .AND. LI <= IEND3 .AND. &
         LJ >= JSTART3 .AND. LJ <= JEND3 .AND. &
         LK >= KSTART3 .AND. LK <= KEND3

      RETURN
      END FUNCTION IS_ON_myPE_plus2layers

!---------------------------------------------------------------------//
      ! logical function to identify a fluid cell
      logical function fluid_at(i,j,k)
         use geometry, only: flag
         implicit none
         integer, intent(in) :: i,j,k
         fluid_at = (flag(i,j,k).eq.1)
      end function fluid_at

      ! logical function to identify a specified pressure inflow cell
      logical function p_flow_at(i,j,k)
      USE geometry
      implicit none
      integer i,j,k
      P_FLOW_AT = flag(i,j,k) .EQ. 10 .OR. &
                  flag(i,j,k) .EQ. 11
      END FUNCTION P_FLOW_AT

      ! logical function to identify a specified pressure outflow cell
      logical function p_outflow_at(i,j,k)
      USE geometry
      implicit none
      integer i,j,k
      P_OUTFLOW_AT= flag(i,j,k) .EQ. 11
      END FUNCTION P_OUTFLOW_AT

! logical function to identify either a specified pressure inflow
! or outflow cell or a fluid cell (simplified check)
! FLUID_AT or P_FLOW_AT (simplified check)
      logical function FLUIDorP_FLOW_AT(i,j,k)
      USE geometry
      implicit none
      integer, intent(in) :: i,j,k
      FLUIDorP_FLOW_AT = (flag(i,j,k) .LE. 11)
      end function FLUIDorP_FLOW_AT

! logical function to identify a specified mass outflow cell
      logical function mass_outflow_at(i,j,k)
      USE geometry
      implicit none
      integer i,j,k
      MASS_OUTFLOW_AT= flag(i,j,k) .EQ. 21
      END FUNCTION MASS_OUTFLOW_AT

! logical function to identify a specified outflow cell
      logical function outflow_at(i,j,k)
      USE geometry
      implicit none
      integer i,j,k
      OUTFLOW_AT  = flag(i,j,k) .EQ. 31
      END FUNCTION OUTFLOW_AT

! logical function to identify any type of flow in/out at cell
! pressure inflow/outflow, mass inflow/outflow or outflow
      logical function flow_at(i,j,k)
      USE geometry
      implicit none
      integer i,j,k
      FLOW_AT     = flag(i,j,k) .GE. 10 .AND. flag(i,j,k) .LE. 31
      END FUNCTION FLOW_AT

! Logical function to identify default walls
      logical function wall_at(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      WALL_AT     = flag(i,j,k) .GE. 100
      END FUNCTION WALL_AT

! Logical function to identify a No-slip wall cell
      logical function ns_wall_at(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      NS_WALL_AT  = flag(i,j,k) .EQ. 100
      end function ns_wall_at

! Logical function to identify a Free-slip wall cell
      logical function fs_wall_at(i,j,k)
         use geometry, only: flag
         implicit none
         integer, intent(in) :: i,j,k
         fs_wall_at = (flag(i,j,k).eq.101)
      end function fs_wall_at

! Logical function to identify a Partial-slip wall cell
      logical function ps_wall_at(i,j,k)
      USE geometry
      implicit none
      integer, intent(in) :: i,j,k
      ps_wall_at  = flag(i,j,k) .EQ. 102
      end function ps_wall_at

      ! Logical function to identify wall ICBC_FLAG
      logical function WALL_ICBC_FLAG(i,j,k)
      use geometry
      use ic
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i,j,k
      WALL_ICBC_FLAG = ICBC_FLAG(i,j,k) .EQ. icbc_no_s  .OR. &
                       ICBC_FLAG(i,j,k) .EQ. icbc_free  .OR. &
                       ICBC_FLAG(i,j,k) .EQ. icbc_pslip .OR. &
                       ICBC_FLAG(i,j,k) .EQ. icbc_cycl  .OR. &
                       ICBC_FLAG(i,j,k) .EQ. icbc_cyclp
      END FUNCTION WALL_ICBC_FLAG

      logical function DEFAULT_WALL_AT(i,j,k)
      use geometry
      use ic
      implicit none
      integer, intent(in) :: i,j,k
      DEFAULT_WALL_AT = (ICBC_FLAG(i,j,k) .NE. icbc_cycl .and. &
                         ICBC_FLAG(i,j,k) .NE. icbc_cyclp)
      END FUNCTION DEFAULT_WALL_AT

! Cyclic
!---------------------------------------------------------------------//
      LOGICAL FUNCTION CYCLIC_AT(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      CYCLIC_AT = (FLAG(i,j,k) .EQ. 106 .OR. FLAG(i,j,k) .EQ. 107)
      END FUNCTION CYCLIC_AT

! logical function to identify cyclic condition at east boundary
      LOGICAL FUNCTION CYCLIC_AT_E(i,j,k)
      USE geometry, only: flag_e
      implicit none
      INTEGER, INTENT(IN) :: i,j,k
      CYCLIC_AT_E   = FLAG_E(i,j,k) .EQ. 2000
      END FUNCTION CYCLIC_AT_E

! logical function to identify cyclic condition at north boundary
      LOGICAL FUNCTION CYCLIC_AT_N(i,j,k)
      USE geometry, only: flag_n
      implicit none
      INTEGER, INTENT(IN) :: i,j,k
      CYCLIC_AT_N   = FLAG_N(i,j,k) .EQ. 2000
      END FUNCTION CYCLIC_AT_N

! logical function to identify cyclic condition at top boundary
      LOGICAL FUNCTION CYCLIC_AT_T(i,j,k)
      USE geometry, only: flag_t
      implicit none
      INTEGER, INTENT(IN) :: i,j,k
      CYCLIC_AT_T   = FLAG_T(i,j,k) .EQ. 2000
      END FUNCTION CYCLIC_AT_T


! Flow boundaries
!---------------------------------------------------------------------//
! identify flow at east boundary
      LOGICAL FUNCTION FLOW_AT_E(i,j,k)
      USE geometry
      implicit none
      INTEGER, INTENT(IN) :: i,j,k
      FLOW_AT_E   = FLAG_E(i,j,k) .GE. 2000 .AND.&
                    FLAG_E(i,j,k) .LE. 2011
      END FUNCTION FLOW_AT_E

! identify specified flow north boundary
      LOGICAL FUNCTION FLOW_AT_N(i,j,k)
      USE geometry
      implicit none
      INTEGER, INTENT(IN) :: i,j,k
      FLOW_AT_N   = FLAG_N(i,j,k) .GE. 2000 .AND.&
                    FLAG_N(i,j,k) .LE. 2011
      END FUNCTION FLOW_AT_N

! identify specified flow top boundary
      LOGICAL FUNCTION FLOW_AT_T(i,j,k)
      USE geometry
      implicit none
      INTEGER, INTENT(IN) :: i,j,k
      FLOW_AT_T   = FLAG_T(i,j,k) .GE. 2000 .AND.&
                    FLAG_T(i,j,k) .LE. 2011
      END FUNCTION FLOW_AT_T

! identify const. pressure flow top boundary
      LOGICAL FUNCTION PFLOW_AT_E(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      PFLOW_AT_E  = FLAG_E(i,j,k) .EQ. 2010 .OR.&
                    FLAG_E(i,j,k) .EQ. 2011
      END FUNCTION PFLOW_AT_E

! identify const. pressure flow north boundary
      LOGICAL FUNCTION PFLOW_AT_N(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      PFLOW_AT_N  = FLAG_N(i,j,k) .EQ. 2010 .OR.&
                    FLAG_N(i,j,k) .EQ. 2011
      END FUNCTION PFLOW_AT_N

! identify const. pressure flow east boundary
      LOGICAL FUNCTION PFLOW_AT_T(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      PFLOW_AT_T  = FLAG_T(i,j,k) .EQ. 2010 .OR.&
                    FLAG_T(i,j,k) .EQ. 2011
      END FUNCTION PFLOW_AT_T

! identify specified flow east boundary
      LOGICAL FUNCTION MFLOW_AT_E(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      MFLOW_AT_E  = FLAG_E(i,j,k) .EQ. 2020 .OR. &
                    FLAG_E(i,j,k) .EQ. 2021 .OR. &
                    FLAG_E(i,j,k) .EQ. 2031
      END FUNCTION MFLOW_AT_E

! identify specified flow north boundary
      LOGICAL FUNCTION MFLOW_AT_N(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      MFLOW_AT_N  = FLAG_N(i,j,k) .EQ. 2020 .OR. &
                    FLAG_N(i,j,k) .EQ. 2021 .OR. &
                    FLAG_N(i,j,k) .EQ. 2031
      END FUNCTION MFLOW_AT_N

! identify specified flow top boundary
      LOGICAL FUNCTION MFLOW_AT_T(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      MFLOW_AT_T  = FLAG_T(i,j,k) .EQ. 2020 .OR. &
                    FLAG_T(i,j,k) .EQ. 2021 .OR. &
                    FLAG_T(i,j,k) .EQ. 2031
      END FUNCTION MFLOW_AT_T


! Functions to identify a impermeable and/or semi-permeable surface at
! indicated boundary (specific type of internal surface)
!---------------------------------------------------------------------//
! Logical function to identify IP (impermeable surface) East of the cell
      logical function ip_at_e(i,j,k)
      use geometry
      implicit none
      integer, intent(in) :: i,j,k
      ip_at_e = flag_e(i,j,k) .LT. 1000
      end function ip_at_e

! Logical function to identify IP (impermeable surface) North of the cell
      logical function ip_at_n(i,j,k)
      use geometry
      implicit none
      integer, intent(in) :: i,j,k
      ip_at_n = flag_n(i,j,k) .LT. 1000
      end function ip_at_n

! Logical function to identify IP (impermeable surface) Top of the cell
      logical function ip_at_t(i,j,k)
      use geometry
      implicit none
      integer, intent(in) :: i,j,k
      ip_at_t = flag_t(i,j,k) .LT. 1000
      end function ip_at_t

! Logical function to identify SP or IP (semi or impermeable surface)
! at east of the cell
      LOGICAL FUNCTION SIP_AT_E(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      SIP_AT_E    = (FLAG_E(i,j,k) .LT. 2000)
      END FUNCTION SIP_AT_E

! Logical function to identify SP or IP (semi or impermeable surface)
! at north of the cell
      LOGICAL FUNCTION SIP_AT_N(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      SIP_AT_N    = (FLAG_N(i,j,k) .LT. 2000)
      END FUNCTION SIP_AT_N

! Logical function to identify SP or IP (semi or impermeable surface)
! at top of the cell
      LOGICAL FUNCTION SIP_AT_T(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      SIP_AT_T    = (FLAG_T(i,j,k) .LT. 2000)
      END FUNCTION SIP_AT_T

! Logical function to identify SP (semi-permeable surface) at east
! of cell
      LOGICAL FUNCTION SP_AT_E(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      SP_AT_E     = (FLAG_E(i,j,k) .LT. 2000) .AND. &
         (FLAG_E(i,j,k) .GE. 1000)
      END FUNCTION SP_AT_E

! Logical function to identify SP (semi-permeable surface) at north
! of cell
      LOGICAL FUNCTION SP_AT_N(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      SP_AT_N     = (FLAG_N(i,j,k) .LT. 2000) .AND. &
         (FLAG_N(i,j,k) .GE. 1000)
      END FUNCTION SP_AT_N

! Logical function to identify SP (semi-permeable surface) at top
! of cell
      LOGICAL FUNCTION SP_AT_T(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      SP_AT_T     = (FLAG_T(i,j,k) .LT. 2000) .AND. &
         (FLAG_T(i,j,k) .GE. 1000)
      END FUNCTION SP_AT_T


! Logical functions concerning general internal surfaces
! Integer functions to return internal surface ID
!---------------------------------------------------------------------//
! Internal surface ID for east face
      INTEGER FUNCTION IS_ID_AT_E(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      IS_ID_AT_E = FLAG_E(i,j,k) - 1000
      END FUNCTION IS_ID_AT_E

! Internal surface ID for north face
      INTEGER FUNCTION IS_ID_AT_N(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      IS_ID_AT_N  = FLAG_N(i,j,k) - 1000
      END FUNCTION IS_ID_AT_N

! Internal surface ID for top face
      INTEGER FUNCTION IS_ID_AT_T(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      IS_ID_AT_T  = FLAG_T(i,j,k) - 1000
      END FUNCTION IS_ID_AT_T

! Logical function to identify IS at East of the cell
      LOGICAL FUNCTION IS_AT_E(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      IS_AT_E     = FLAG_E(i,j,k) .LT. 2000
      END FUNCTION IS_AT_E

! Logical function to identify IS at North of the cell
      LOGICAL FUNCTION IS_AT_N(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      IS_AT_N     = FLAG_N(i,j,k) .LT. 2000
      END FUNCTION IS_AT_N

! Logical function to identify IS at Top of the cell
      LOGICAL FUNCTION IS_AT_T(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      IS_AT_T     = FLAG_T(i,j,k) .LT. 2000
      END FUNCTION IS_AT_T

! Logical function to identify No IS at East of the cell
      LOGICAL FUNCTION NO_IS_AT_E(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      NO_IS_AT_E  = FLAG_E(i,j,k) .GE. 2000
      END FUNCTION NO_IS_AT_E

! Logical function to identify No IS at North of the cell
      LOGICAL FUNCTION NO_IS_AT_N(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      NO_IS_AT_N  = FLAG_N(i,j,k) .GE. 2000
      END FUNCTION NO_IS_AT_N

! Logical function to identify No IS at Top of the cell
      LOGICAL FUNCTION NO_IS_AT_T(i,j,k)
      USE geometry
      implicit none
      INTEGER i,j,k
      NO_IS_AT_T  = FLAG_T(i,j,k) .GE. 2000
      END FUNCTION NO_IS_AT_T

! Misc
!---------------------------------------------------------------------//

! Function that returns the maximum of zero or input
      DOUBLE PRECISION FUNCTION ZMAX(XXX)
      USE param1, only: zero
      implicit none
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

      if (i.eq.imin1 .and. (wall_at(i,j,k) .or. flow_at(i,j,k)) ) then
         iminus = i
      else
         iminus = max(istart3, i-1)
      end if

      if (cyclic_x .and. (iminus < imin1) .and. .not.wall_at(i,j,k)) iminus = imax1

      end function iminus


      ! ************************************************

      integer function iplus(i,j,k)
      use compar, only: iend3
      use geometry, only  : imin1,imax1,cyclic_x
      implicit none
      integer, intent(in) :: i,j,k

      if (i.eq.imax1 .and. (wall_at(i,j,k) .or. flow_at(i,j,k)) ) then
         iplus = i
      else
         iplus = min(iend3, i+1)
      end if

      if (cyclic_x .and. (iplus > imax1) .and. .not.wall_at(i,j,k)) iplus = imin1

      end function iplus

      ! ************************************************

      integer function jminus(i,j,k)
      use compar, only: jstart3
      use geometry, only  : jmin1,jmax1,cyclic_y
      implicit none
      integer, intent(in) :: i,j,k

      if (j.eq.jmin1 .and. (wall_at(i,j,k) .or. flow_at(i,j,k)) ) then
         jminus = j
      else
         jminus = max(jstart3, j-1)
      end if

      if (cyclic_y .and. (jminus < jmin1) .and. .not.wall_at(i,j,k)) jminus = jmax1

      end function jminus

      ! ************************************************

      integer function jplus(i,j,k)
      use compar, only: jend3
      use geometry, only  : jmin1,jmax1,cyclic_y
      implicit none
      integer, intent(in) :: i,j,k

      if (j.eq.jmax1 .and. (wall_at(i,j,k) .or. flow_at(i,j,k)) ) then
         jplus = j
      else
         jplus = min(jend3, j+1)
      end if

      if (cyclic_y .and. (jplus > jmax1) .and. .not.wall_at(i,j,k)) jplus = jmin1

      end function jplus

      ! ************************************************

      integer function kminus(i,j,k)

      use compar, only: kstart3
      use geometry, only  : kmin1,kmax1,cyclic_z
      implicit none
      integer, intent(in) :: i,j,k

      if (k.eq.kmin1 .and. (wall_at(i,j,k) .or. flow_at(i,j,k)) ) then
         kminus = k
      else
         kminus = max(kstart3, k-1)
      end if

      if (cyclic_z .and. (kminus < kmin1) .and. .not.wall_at(i,j,k)) kminus = kmax1

      end function kminus

      ! ************************************************

      integer function kplus(i,j,k)

      use compar, only: kend3
      use geometry, only  : kmin1,kmax1,cyclic_z
      implicit none
      integer, intent(in) :: i,j,k

      if (k.eq.kmax1 .and. (wall_at(i,j,k) .or. flow_at(i,j,k)) ) then
         kplus = k
      else
         kplus = min(kend3, k+1)
      end if

      if (cyclic_z .and. (kplus > kmax1) .and. .not.wall_at(i,j,k)) kplus = kmin1

      end function kplus

      ! ************************************************

      integer function ieast(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k

      if (wall_at(iplus(i,j,k),j,k)) then
         ieast = i
      else
         ieast = iplus(i,j,k)
      end if

      end function ieast

      ! ************************************************

      integer function iwest(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k

      if (wall_at(iminus(i,j,k),j,k)) then
         iwest = i
      else
         iwest = iminus(i,j,k)
      end if

      end function iwest

      ! ************************************************

      integer function jnorth(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k

      if (wall_at(i,jplus(i,j,k),k)) then
         jnorth = j
      else
         jnorth = jplus(i,j,k)
      end if

      end function jnorth

      ! ************************************************

      integer function jsouth(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k

      if (wall_at(i,jminus(i,j,k),k)) then
         jsouth = j
      else
         jsouth = jminus(i,j,k)
      end if

      end function jsouth

      ! ************************************************

      integer function ktop(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k
      integer             :: kp

      kp = kplus(i,j,k)

      if (wall_at(i,j,kp)) then
         ktop = k
      else
         ktop = kp
      end if

      end function ktop

      ! ************************************************

      integer function kbot(i,j,k)

      implicit none
      integer, intent(in) :: i,j,k
      integer             :: km

      km = kminus(i,j,k)

      if (wall_at(i,j,km)) then
         kbot = k
      else
         kbot = km
      end if

      end function kbot

! ************************************************

      integer function im1(i)

      use compar, only: nodesi, istart3
      use geometry, only  : imin1, imin2, imax1, cyclic_x
      implicit none
      integer, intent(in) :: i

      if(i==imin2) then
         im1 = imin2
      elseif(cyclic_x .and. nodesi == 1 .and. i==imin1 ) then
         im1 = imax1
      else
         im1 = max(istart3, i-1)
      endif

      end function im1

! ************************************************

      integer function ip1(i)

      use compar, only: nodesi, iend3
      use geometry, only  : imin1, imax2, imax1, cyclic_x
      implicit none
      integer, intent(in) :: i

      if(i==imax2) then
         ip1 = imax2
      elseif(cyclic_x .and. nodesi == 1 .and. i==imax1 ) then
         ip1 = imin1
      else
         ip1 = min(iend3, i+1)
      endif

      end function ip1

! ************************************************

      integer function jm1(j)

      use compar, only: nodesj, jstart3
      use geometry, only  : jmin1, jmin2, jmax1, cyclic_y
      implicit none
      integer, intent(in) :: j

      if(j==jmin2) then
         jm1 = jmin2
      elseif(cyclic_y .and. nodesj == 1 .and. j==jmin1 ) then
         jm1 = jmax1
      else
         jm1 = max(jstart3, j-1)
      endif

      end function jm1

! ************************************************

      integer function jp1(j)

      use compar, only: nodesj, jend3
      use geometry, only  : jmin1, jmax2, jmax1, cyclic_y
      implicit none
      integer, intent(in) :: j

      if(j==jmax2) then
         jp1 = jmax2
      elseif(cyclic_y .and. nodesj == 1 .and. j==jmax1 ) then
         jp1 = jmin1
      else
         jp1 = min(jend3, j+1)
      endif

      end function jp1

! ************************************************

      integer function km1(k)

      use compar, only: nodesk, kstart3
      use geometry, only: no_k, kmin1, kmin2, kmax1, cyclic_z
      implicit none
      integer, intent(in) :: k

      if(no_k) then
         km1 = kmax1
      elseif(k==kmin2) then
         km1 = kmin2
      elseif(cyclic_z .and. nodesk == 1 .and. k==kmin1 ) then
         km1 = kmax1
      else
         km1 = max(kstart3, k-1)
      endif

      end function km1

! ************************************************

      integer function kp1(k)

      use compar, only: nodesk, kend3
      use geometry, only: no_k, kmin1, kmax2, kmax1, cyclic_z
      implicit none
      integer, intent(in) :: k

      if(no_k)then
         kp1 = kmin1
      elseif(k==kmax2) then
         kp1 = kmax2
      elseif(cyclic_z .and. nodesk == 1 .and. k==kmax1 ) then
         kp1 = kmin1
      else
         kp1 = min(kend3, k+1)
      endif

      end function kp1
