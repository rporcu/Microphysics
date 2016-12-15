MODULE SET_FLAGS_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_FLAGS                                               C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!  Purpose: This module assigns a flag to a cell to identify its type. C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_FLAGS(flag)

      use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3
      use ic, only: PINF_, POUT_, MINF_, MOUT_, OUTF_, NSW_, FSW_
      use ic, only: NSW_
      use geometry, only: ijkmax3
      use compar, only: iend3, mype, pe_io
      use param1, only: undefined_i

      IMPLICIT NONE

      integer, intent(inout) :: flag(istart3:iend3,jstart3:jend3,kstart3:kend3,4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: i, j, k, ib, jb, kb
      integer, allocatable :: arr1(:)
!-----------------------------------------------

!  Cell flag definitions
!  FLAG1   FLAG0   BC_TYPE        Cell type
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
      DO i = istart3, iend3
         DO j = jstart3, jend3
            DO k = kstart3, kend3

              select case (flag(i,j,k,1))
                case (PINF_, POUT_, MINF_, MOUT_, OUTF_)

                ib = min( iend3, max (istart3, i+1) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k) )
                if (flag(ib,jb,kb,1) == NSW_) &
                   flag(ib,jb,kb,1) = FSW_

                ib = min( iend3, max (istart3, i-1) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k) )
                if (flag(ib,jb,kb,1) == NSW_) &
                   flag(ib,jb,kb,1) = FSW_

                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j+1) )
                kb = min( kend3, max (kstart3, k) )
                if (flag(ib,jb,kb,1) == NSW_) &
                   flag(ib,jb,kb,1) = FSW_

                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j-1) )
                kb = min( kend3, max (kstart3, k) )
                if (flag(ib,jb,kb,1) == NSW_) &
                   flag(ib,jb,kb,1) = FSW_

                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k+1) )
                if (flag(ib,jb,kb,1) == NSW_) &
                   flag(ib,jb,kb,1) = FSW_

                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k-1) )
                if (flag(ib,jb,kb,1) == NSW_) &
                   flag(ib,jb,kb,1) = FSW_

              end select

            ENDDO
          ENDDO
      ENDDO
! ----------------------------------------------------------------<<<
! Define the numerical value of the variable flag for all cells based
! on the corresponding character value of flag0.  By this point the
! flag0 has been defined in all cells
! ---------------------------------------------------------------->>>
      do k = kstart3, kend3
         do j = jstart3, jend3
           do i = istart3, iend3

! Initialize cell face flags.  UNDEFINED_I should be a large +ve value.
              flag(i,j,k,2) = UNDEFINED_I
              flag(i,j,k,3) = UNDEFINED_I
              flag(i,j,k,4) = UNDEFINED_I

          end do
        end do
      end do

      IF (MYPE.EQ.PE_IO) THEN
         ALLOCATE (ARR1(IJKMAX3))
      ELSE
         ALLOCATE (ARR1(1))
      ENDIF

      ! CALL GATHER(FLAG,ARR1,ROOT)
      ! CALL SCATTER(FLAG,ARR1,ROOT)

      DEALLOCATE (ARR1)


      RETURN

      END SUBROUTINE SET_FLAGS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_FLAGS1                                             C
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
!  Variables modified: FLAG_E, FLAG_N, FLAG_T                          C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_FLAGS1(flag)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1   , only: undefined_i
      USE geometry , only: imax2,jmax2,kmax2,imax3,jmax3,kmax3
      USE compar   , only: istart3,iend3,jstart3,jend3,kstart3,kend3
      USE functions, only: iminus, iplus, jminus, jplus, kminus, kplus

      implicit none

      integer, intent(inout) :: flag(istart3:iend3,jstart3:jend3,kstart3:kend3,4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: i,j,k
!-----------------------------------------------

      do k = kstart3, kend3
         do j = jstart3, jend3
           do i = istart3, iend3
                print *,'IJK ',i,j,k,flag(i,j,k,1)

! If the flag is greater than or equal to 2000, there is no
! internal surface.
         IF (flag(i,j,k,1)>=100) THEN

! ---------------------------------------------------------------->>>
! the default is equivalent to an impermeable surface and these cells
! will be treated as such in the momentum routines
                write(6,*) 'setting 2'
                flush(6)
            FLAG(i,j,k,2) = 0
                write(6,*) 'set     2'
                flush(6)
            FLAG(i,j,k,3) = 0
            FLAG(i,j,k,4) = 0
                print *,'MINUS ',i,j,k,iminus(i,j,k)
                flush(6)
            FLAG(iminus(i,j,k),j,k,2) = 0
                write(6,*) jminus(i,j,k)
                flush(6)
            FLAG(i,jminus(i,j,k),k,3) = 0
                write(6,*) kminus(i,j,k)
                flush(6)
            FLAG(i,j,kminus(i,j,k),4) = 0

                print *,'here now '

            if(flag(i,j,k,1) == 106 .or. flag(i,j,k,1) == 107) then
! make the upper (E, N, T) boundary permeable
               IF (I == IMAX2) THEN
                  IF ((J/=1.AND.J/=0) .AND. (J/=JMAX2.AND.J/=JMAX3)) THEN
                     IF ((K/=1.AND.K/=0) .AND. (K/=KMAX2.AND.K/=KMAX3)) THEN
                        IF(flag(iminus(i,j,k),j,k,1) < 100) &
                           FLAG(iminus(i,j,k),j,k,2) = 2000
                     ENDIF
                  ENDIF
               ENDIF
               IF (J == JMAX2) THEN
                  IF ((I/=1.AND.I/=0) .AND. (I/=IMAX2.AND.I/=IMAX3)) THEN
                     IF ((K/=1.AND.K/=0) .AND. (K/=KMAX2.AND.K/=KMAX3)) THEN
                        IF(flag(i,jminus(i,j,k),k,1)<100) &
                           FLAG(i,jminus(i,j,k),k,3) = 2000
                     ENDIF
                  ENDIF
                ENDIF
               IF (K == KMAX2) THEN
                  IF ((J/=1.AND.J/=0) .AND. (J/=JMAX2.AND.J/=JMAX3)) THEN
                     IF ((I/=1.AND.I/=0) .AND. (I/=IMAX2.AND.I/=IMAX3) .AND. &
                        flag(i,j,kminus(i,j,k),1)<100) &
                        FLAG(i,j,kminus(i,j,k),4) = 2000
                  ENDIF
               ENDIF

                print *,'here end '

            ENDIF

! ----------------------------------------------------------------<<<
         ELSEIF (flag(i,j,k,1)==1) THEN
! ---------------------------------------------------------------->>>

            IF ( flag(iminus(i,j,k),j,k,1) < 100 .AND. &
               FLAG(iminus(i,j,k),j,k,2)==UNDEFINED_I) &
               FLAG(iminus(i,j,k),j,k,2) = 2000 + FLAG(iminus(i,j,k),j,k,1)

            IF ( flag(iplus(i,j,k),j,k,1) < 100 .AND. &
               FLAG(i,j,k,2)==UNDEFINED_I) &
               FLAG(i,j,k,2) = 2000 + FLAG(iplus(i,j,k),j,k,1)

            IF ( flag(i,jminus(i,j,k),k,1) < 100 .AND. &
               FLAG(i,jminus(i,j,k),k,3)==UNDEFINED_I) &
               FLAG(i,jminus(i,j,k),k,3) = 2000 + FLAG(i,jminus(i,j,k),k,1)

            IF ( flag(i,jplus(i,j,k),k,1) < 100 .AND. &
               FLAG(i,j,k,3)==UNDEFINED_I) &
               FLAG(i,j,k,3) = 2000 + FLAG(i,jplus(i,j,k),k,1)

            IF ( flag(i,j,kminus(i,j,k),1) < 100 .AND. &
               FLAG(i,j,kminus(i,j,k),4)==UNDEFINED_I) &
               FLAG(i,j,kminus(i,j,k),4) = 2000 + FLAG(i,j,kminus(i,j,k),1)

            IF ( flag(i,j,kplus(i,j,k),1) < 100 .AND. &
               FLAG(i,j,k,4)==UNDEFINED_I) &
               FLAG(i,j,k,4) = 2000 + FLAG(i,j,kplus(i,j,k),1)


         ENDIF

          end do
        end do
     end do

! Clean up edge cases
     do k = kstart3, kend3
         do j = jstart3, jend3
            do i = istart3, iend3
               if(flag(i,j,k,2) == UNDEFINED_I) flag(i,j,k,2) = 0
               if(flag(i,j,k,3) == UNDEFINED_I) flag(i,j,k,3) = 0
               if(flag(i,j,k,4) == UNDEFINED_I) flag(i,j,k,4) = 0
               write(1200,"(4(i4),3(I11))") i,j,k,flag(i,j,k,:)
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

      RETURN
      END SUBROUTINE SET_FLAGS1
END MODULE SET_FLAGS_MODULE
