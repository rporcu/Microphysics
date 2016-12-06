!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_FLAGS                                               C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!  Purpose: This module assigns a flag to a cell to identify its type. C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_FLAGS

      use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3
      use ic, only: icbc_p_inf, icbc_p_out, icbc_m_inf, icbc_m_out, icbc_outfl, icbc_fluid, icbc_no_s, icbc_free, icbc_pslip
      use ic, only: icbc_cyclp, icbc_no_s, icbc_cycl
      use geometry, only: flag, ijkmax3
      use compar, only: iend3, mype, pe_io
      use param1, only: undefined_i

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: i, j, k, ib, jb, kb, inc
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

              select case (mod(flag(i,j,k,0),1000))
                case (icbc_p_inf, icbc_p_out, icbc_m_inf, icbc_m_out, icbc_outfl)

                ib = min( iend3, max (istart3, i+1) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k) )
                inc = flag(ib,jb,kb,0) - mod(flag(ib,jb,kb,0),1000)
                if (flag(ib,jb,kb,0) == icbc_no_s) &
                   flag(ib,jb,kb,0) = icbc_free+inc

                ib = min( iend3, max (istart3, i-1) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k) )
                inc = flag(ib,jb,kb,0) - mod(flag(ib,jb,kb,0),1000)
                if (flag(ib,jb,kb,0) == icbc_no_s) &
                   flag(ib,jb,kb,0) = icbc_free+inc

                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j+1) )
                kb = min( kend3, max (kstart3, k) )
                inc = flag(ib,jb,kb,0) - mod(flag(ib,jb,kb,0),1000)
                if (flag(ib,jb,kb,0) == icbc_no_s) &
                   flag(ib,jb,kb,0) = icbc_free+inc

                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j-1) )
                kb = min( kend3, max (kstart3, k) )
                inc = flag(ib,jb,kb,0) - mod(flag(ib,jb,kb,0),1000)
                if (flag(ib,jb,kb,0) == icbc_no_s) &
                   flag(ib,jb,kb,0) = icbc_free+inc

                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k+1) )
                inc = flag(ib,jb,kb,0) - mod(flag(ib,jb,kb,0),1000)
                if (flag(ib,jb,kb,0) == icbc_no_s) &
                   flag(ib,jb,kb,0) = icbc_free+inc

                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k-1) )
                inc = flag(ib,jb,kb,0) - mod(flag(ib,jb,kb,0),1000)
                if (flag(ib,jb,kb,0) == icbc_no_s) &
                   flag(ib,jb,kb,0) = icbc_free+inc

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

              select case (mod(flag(i,j,k,0),1000))

                CASE (icbc_fluid);  FLAG(i,j,k,1) = 1
                CASE (icbc_p_inf);  FLAG(i,j,k,1) = 10
                CASE (icbc_p_out);  FLAG(i,j,k,1) = 11
                CASE (icbc_m_inf);  FLAG(i,j,k,1) = 20
                CASE (icbc_m_out);  FLAG(i,j,k,1) = 21
                CASE (icbc_outfl);  FLAG(i,j,k,1) = 31
                CASE (icbc_no_s);   FLAG(i,j,k,1) = 100
                CASE (icbc_free);   FLAG(i,j,k,1) = 101
                CASE (icbc_pslip);  FLAG(i,j,k,1) = 102
                CASE (icbc_cycl);   FLAG(i,j,k,1) = 106
                CASE (icbc_cyclp);  FLAG(i,j,k,1) = 107
                CASE DEFAULT

              end select

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

      SUBROUTINE SET_FLAGS1

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1   , only: undefined_i
      USE geometry , only: imax2,jmax2,kmax2,imax3,jmax3,kmax3,flag
      USE compar   , only: istart3,iend3,jstart3,jend3,kstart3,kend3
      USE functions, only: iminus, iplus, jminus, jplus, kminus, kplus
      USE functions, only: fluid_at

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: i,j,k
!-----------------------------------------------

      do k = kstart3, kend3
         do j = jstart3, jend3
           do i = istart3, iend3
! If the flag is greater than or equal to 2000, there is no
! internal surface.
         IF (flag(i,j,k,1)>=100) THEN
! ---------------------------------------------------------------->>>
! the default is equivalent to an impermeable surface and these cells
! will be treated as such in the momentum routines
            FLAG(i,j,k,2) = 0
            FLAG(i,j,k,3) = 0
            FLAG(i,j,k,4) = 0
            FLAG(iminus(i,j,k),j,k,2) = 0
            FLAG(i,jminus(i,j,k),k,3) = 0
            FLAG(i,j,kminus(i,j,k),4) = 0

            if(flag(i,j,k,1) == 106 .or. flag(i,j,k,1) == 107) then
! make the upper (E, N, T) boundary permeable
               IF (I == IMAX2) THEN
                  IF ((J/=1.AND.J/=0.) .AND. (J/=JMAX2.AND.J/=JMAX3)) THEN
                     IF ((K/=1.AND.K/=0) .AND. (K/=KMAX2.AND.K/=KMAX3)) THEN
                        IF(flag(iminus(i,j,k),j,k,1) <100) &
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
                  IF ((J/=1.AND.J/=0.) .AND. (J/=JMAX2.AND.J/=JMAX3)) THEN
                     IF ((I/=1.AND.I/=0) .AND. (I/=IMAX2.AND.I/=IMAX3) .AND. &
                        flag(i,j,kminus(i,j,k),1)<100) &
                        FLAG(i,j,kminus(i,j,k),4) = 2000
                  ENDIF
               ENDIF

            ENDIF

! ----------------------------------------------------------------<<<
         ELSEIF (fluid_at(i,j,k)) THEN
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
