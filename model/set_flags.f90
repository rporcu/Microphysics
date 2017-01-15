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

      subroutine set_flags(slo,shi,flag)

      use ic, only: PINF_, POUT_, MINF_, MOUT_, OUTF_, NSW_, FSW_
      use ic, only: NSW_
      use compar, only: mype, pe_io
      use param1, only: undefined_i

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      integer, intent(inout) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

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
      do k = slo(3),shi(3)
         do j = slo(2),shi(2)
            do i = slo(1),shi(1)

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
! ----------------------------------------------------------------<<<
! Define the numerical value of the variable flag for all cells based
! on the corresponding character value of flag0.  By this point the
! flag0 has been defined in all cells
! ---------------------------------------------------------------->>>
      flag(:,:,:,2:4) = UNDEFINED_I

      end subroutine set_flags

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_flagS1                                             C
!  Purpose: Assign IP flag to the faces of wall cells                  C
!                                                                      C
!  Notes: This routine may still leave flag_e, flag_n and flag_t       C
!         undefined in some boundary cells                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified: flag_E, flag_N, flag_T                          C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      subroutine set_flags1(slo,shi,flag)

      USE param1   , only: is_undefined
      USE geometry , only: domlo,domhi
      USE functions, only: iminus, iplus, jminus, jplus, kminus, kplus

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)

      integer, intent(inout) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      integer :: i,j,k

      do k = slo(3),shi(3)
         do j = slo(2),shi(2)
           do i = slo(1),shi(1)

! If the flag is greater than or equal to 2000, there is no
! internal surface.
         IF (flag(i,j,k,1)>=100) then

! ---------------------------------------------------------------->>>
! the default is equivalent to an impermeable surface and these cells
! will be treated as such in the momentum routines
            flag(i,j,k,2) = 0
            flag(i,j,k,3) = 0
            flag(i,j,k,4) = 0
            flag(iminus(i,j,k),j,k,2) = 0
            flag(i,jminus(i,j,k),k,3) = 0
            flag(i,j,kminus(i,j,k),4) = 0

            if(flag(i,j,k,1) == 106 .or. flag(i,j,k,1) == 107) then

               ! make the upper (E, N, T) boundary permeable
               if (I == domhi(1)+1) then
                  if ( (J>=(domlo(2)-1)) .and. (J<=(domhi(2)+1)) .and. &
                       (K>=(domlo(3)-1)) .and. (K<=(domhi(3)+1)) .and. &
                       (flag(iminus(i,j,k),j,k,1) < 100) ) then
                           flag(iminus(i,j,k),j,k,2) = 2000
                  end if
               end if

               if (J == domhi(2)+1) then
                  if ( (I>=(domlo(1)-1)) .and. (I<=(domhi(1)+1)) .and. &
                       (K>=(domlo(3)-1)) .and. (K<=(domhi(3)+1)) .and. &
                       (flag(i,jminus(i,j,k),k,1) < 100) ) then
                           flag(i,jminus(i,j,k),k,3) = 2000
                  end if
               end if

               if (K == domhi(3)+1) then
                  if ( (J>=(domlo(2)-1)) .and. (J<=(domhi(2)+1)) .and. &
                       (I>=(domlo(1)-1)) .and. (I<=(domhi(1)+1)) .and. &
                       (flag(i,j,kminus(i,j,k),1) < 100) ) then
                        flag(i,j,kminus(i,j,k),4) = 2000
                  end if
               end if

            ENDIF

! ----------------------------------------------------------------<<<
         ELSEIF (flag(i,j,k,1)==1) then
! ---------------------------------------------------------------->>>

            IF ( flag(iminus(i,j,k),j,k,1) < 100 .and. &
               IS_UNDEFINED(flag(iminus(i,j,k),j,k,2))) &
               flag(iminus(i,j,k),j,k,2) = 2000 + flag(iminus(i,j,k),j,k,1)

            IF ( flag(iplus(i,j,k),j,k,1) < 100 .and. &
               IS_UNDEFINED(flag(i,j,k,2))) &
               flag(i,j,k,2) = 2000 + flag(iplus(i,j,k),j,k,1)

            IF ( flag(i,jminus(i,j,k),k,1) < 100 .and. &
               IS_UNDEFINED(flag(i,jminus(i,j,k),k,3))) &
               flag(i,jminus(i,j,k),k,3) = 2000 + flag(i,jminus(i,j,k),k,1)

            IF ( flag(i,jplus(i,j,k),k,1) < 100 .and. &
               IS_UNDEFINED(flag(i,j,k,3))) &
               flag(i,j,k,3) = 2000 + flag(i,jplus(i,j,k),k,1)

            IF ( flag(i,j,kminus(i,j,k),1) < 100 .and. &
               IS_UNDEFINED(flag(i,j,kminus(i,j,k),4))) &
               flag(i,j,kminus(i,j,k),4) = 2000 + flag(i,j,kminus(i,j,k),1)

            IF ( flag(i,j,kplus(i,j,k),1) < 100 .and. &
               IS_UNDEFINED(flag(i,j,k,4))) &
               flag(i,j,k,4) = 2000 + flag(i,j,kplus(i,j,k),1)

         ENDIF

          end do
        end do
     end do

! Clean up edge cases
      do k = slo(3),shi(3)
         do j = slo(2),shi(2)
            do i = slo(1),shi(1)
               if(IS_UNDEFINED(flag(i,j,k,2))) flag(i,j,k,2) = 0
               if(IS_UNDEFINED(flag(i,j,k,3))) flag(i,j,k,3) = 0
               if(IS_UNDEFINED(flag(i,j,k,4))) flag(i,j,k,4) = 0
            end do
         end do
      end do
! ----------------------------------------------------------------<<<

! Fill the ghost layers using gather and scatter
      ! call gather( flag_e, flag_temp )
      ! flag_e = flag_temp
      ! call scatter( flag_e, flag_temp )
      ! call gather( flag_n, flag_temp )
      ! flag_n = flag_temp
      ! call scatter( flag_n, flag_temp )
      ! call gather( flag_t, flag_temp )
      ! flag_t = flag_temp
      ! call scatter( flag_t, flag_temp )

! deallocate storage of temporary flag arrays
      ! deallocate( flag_temp )
      ! call send_recv(flag_t,1)
      ! call send_recv(flag_n,1)
      ! call send_recv(flag_e,1)

      end subroutine set_flags1

end module set_flags_module
