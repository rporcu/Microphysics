module set_flags_module

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_flagS                                               C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!  Purpose: This module assigns a flag to a cell to identify its type. C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      subroutine set_flags(slo,shi,lo,hi,flag)

      use ic, only: PINF_, POUT_, MINF_, MOUT_, OUTF_, NSW_, FSW_
      use ic, only: NSW_
      use param1, only: undefined_i

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      integer, intent(inout) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),1)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: i, j, k, ib, jb, kb
!-----------------------------------------------

!  Cell flag definitions
!  flag1   flag0   BC_TYPE        Cell type
!  ----- --------- -------        ---------
!   1       .        -            Cell containing gas or solids or both
!  10       p      P_INFLOW       Specified pressure inflow cell
!  11       P      P_OUTFLOW      Specified pressure outflow cell
!  20       I      MASS_INFLOW    Specified mass flux inflow cell
!  21       O      MASS_OUTFLOW   Specified mass flux outflow cell
!  31       o      OUTFLOW        outflow cell
! 100       W      NO_SLIP_WALL   Internal/external wall with no-slip b.c.
! 101       S      FREE_SLIP_WALL Internal/external wall with free-slip
! 102       s      PAR_SLIP_WALL  Internal/external wall with partial-slip b.c.
! 106       c      CYCLIC         Cyclic b.c.
! 107       C      CYCLIC_PD      Cyclic b.c. with pressure drop
! Flag values greater than 100 are considered to be wall cells
! (see function.inc).

! make the wall cells adjacent to flow boundaries free-slip wall to
! avoid unphysical strain rates in fluid cells adjacent to the flow
! boundary
! ---------------------------------------------------------------->>>
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),shi(1)

              select case (flag(i,j,k,1))
                case (PINF_, POUT_, MINF_, MOUT_, OUTF_)

                ib = min( shi(1), max (slo(1), i+1) )
                jb = min( shi(2), max (slo(2), j) )
                kb = min( shi(3), max (slo(3), k) )
                if (flag(ib,jb,kb,1) == NSW_) &
                    flag(ib,jb,kb,1) = FSW_

                ib = min( shi(1), max (slo(1), i-1) )
                jb = min( shi(2), max (slo(2), j) )
                kb = min( shi(3), max (slo(3), k) )
                if (flag(ib,jb,kb,1) == NSW_) &
                    flag(ib,jb,kb,1) = FSW_

                ib = min( shi(1), max (slo(1), i) )
                jb = min( shi(2), max (slo(2), j+1) )
                kb = min( shi(3), max (slo(3), k) )
                if (flag(ib,jb,kb,1) == NSW_) &
                    flag(ib,jb,kb,1) = FSW_

                ib = min( shi(1), max (slo(1), i) )
                jb = min( shi(2), max (slo(2), j-1) )
                kb = min( shi(3), max (slo(3), k) )
                if (flag(ib,jb,kb,1) == NSW_) &
                   flag(ib,jb,kb,1) = FSW_

                ib = min( shi(1), max (slo(1), i) )
                jb = min( shi(2), max (slo(2), j) )
                kb = min( shi(3), max (slo(3), k+1) )
                if (flag(ib,jb,kb,1) == NSW_) &
                    flag(ib,jb,kb,1) = FSW_

                ib = min( shi(1), max (slo(1), i) )
                jb = min( shi(2), max (slo(2), j) )
                kb = min( shi(3), max (slo(3), k-1) )
                if (flag(ib,jb,kb,1) == NSW_) &
                    flag(ib,jb,kb,1) = FSW_

                if (ib .lt. slo(1) .or. ib .gt. shi(1)) print *,'BAD IB ',ib
                if (jb .lt. slo(2) .or. jb .gt. shi(2)) print *,'BAD IB ',jb
                if (kb .lt. slo(3) .or. kb .gt. shi(3)) print *,'BAD IB ',kb

              end select

            ENDDO
          ENDDO
      ENDDO

      end subroutine set_flags

end module set_flags_module
