!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_FLAGS                                               C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!  Purpose: This module assigns a flag to a cell to identify its type. C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_FLAGS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE ic
      USE bc
      USE physprop
      USE funits
      USE compar
      USE boundfunijk
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: i, j, k, ib, jb, kb
      integer, allocatable :: arr1(:)
!-----------------------------------------------

!  Cell flag definitions
!  FLAG  ICBC_FLAG BC_TYPE        Cell type
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
      DO i = istart4, iend4
         DO j = jstart4, jend4
            DO k = kstart4, kend4
              SELECT CASE (icbc_flag(i,j,k))
                CASE (icbc_p_inf, icbc_p_out, icbc_m_inf, icbc_m_out, icbc_outfl)
 
                ib = min( iend3, max (istart3, i+1) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k) )

                if (icbc_flag(ib,jb,kb) == icbc_no_s) icbc_flag(ib,jb,kb) = icbc_free
 
                ib = min( iend3, max (istart3, i-1) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k) )
                if (icbc_flag(ib,jb,kb) == icbc_no_s) icbc_flag(ib,jb,kb) = icbc_free
 
                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j+1) )
                kb = min( kend3, max (kstart3, k) )
                if (icbc_flag(ib,jb,kb) == icbc_no_s) icbc_flag(ib,jb,kb) = icbc_free
 
                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j-1) )
                kb = min( kend3, max (kstart3, k) )
                if (icbc_flag(ib,jb,kb) == icbc_no_s) icbc_flag(ib,jb,kb) = icbc_free
 
                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k+1) )
                if (icbc_flag(ib,jb,kb) == icbc_no_s) icbc_flag(ib,jb,kb) = icbc_free
 
                ib = min( iend3, max (istart3, i) )
                jb = min( jend3, max (jstart3, j) )
                kb = min( kend3, max (kstart3, k-1) )
                if (icbc_flag(ib,jb,kb) == icbc_no_s) icbc_flag(ib,jb,kb) = icbc_free
              END SELECT
            ENDDO
          ENDDO
      ENDDO
! ----------------------------------------------------------------<<<
! Define the numerical value of the variable flag for all cells based
! on the corresponding character value of icbc_flag.  By this point the
! icbc_flag has been defined in all cells
! ---------------------------------------------------------------->>>
      do k = kstart3, kend3
         do j = jstart3, jend3
           do i = istart3, iend3

         SELECT CASE (ICBC_FLAG(i,j,k))
         CASE (icbc_fluid)
            FLAG(i,j,k) = 1
         CASE (icbc_p_inf)
            FLAG(i,j,k) = 10
         CASE (icbc_p_out)
            FLAG(i,j,k) = 11
         CASE (icbc_m_inf)
            FLAG(i,j,k) = 20
         CASE (icbc_m_out)
            FLAG(i,j,k) = 21
         CASE (icbc_outfl)
            FLAG(i,j,k) = 31
         CASE (icbc_no_s)
            FLAG(i,j,k) = 100
         CASE (icbc_free)
            FLAG(i,j,k) = 101
         CASE (icbc_pslip)
            FLAG(i,j,k) = 102
         CASE (icbc_cycl)
            FLAG(i,j,k) = 106
         CASE (icbc_cyclp)
            FLAG(i,j,k) = 107
         CASE DEFAULT

! Access to only one thread at a time
!           IF (DMP_LOG)WRITE (UNIT_LOG, 1000) i,j,k, ICBC_FLAG(i,j,k)
!           call mfix_exit(myPE)
         END SELECT
! ----------------------------------------------------------------<<<

! Initialize cell face flags.  UNDEFINED_I should be a large +ve value.
         FLAG_E(i,j,k) = UNDEFINED_I
         FLAG_N(i,j,k) = UNDEFINED_I
         FLAG_T(i,j,k) = UNDEFINED_I
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
 1000 FORMAT(/1X,70('*')//' From: SET_FLAGS',/&
         ' Message: ICBC_FLAG(',I3,') = ',&
         A3,' is illegal',/1X,70('*')/)

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
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE bc
      USE physprop
      USE funits
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP
      INTEGER :: I, J, K
!
      INTEGER, DIMENSION(:), allocatable :: FLAG_TEMP
      INTEGER :: flag_size
!-----------------------------------------------


! Allocate storage for temporary flag arrays
      flag_size = ijkmax3
      if (myPE.eq.root) then
          flag_size = ijkmax3
      endif
      allocate( flag_temp(flag_size) )


      do k = kstart3, kend3
         do j = jstart3, jend3
           do i = istart3, iend3

           ijk = funijk(i,j,k)
           IMJK = funijk(iminus(i,j,k),j,k)
           IJMK = funijk(i,jminus(i,j,k),k)
           IJKM = funijk(i,j,kminus(i,j,k))
           IPJK = funijk(iplus(i,j,k),j,k)
           IJPK = funijk(i,jplus(i,j,k),k)
           IJKP = funijk(i,j,kplus(i,j,k))

           IF(.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

! If the flag is greater than or equal to 2000, there is no
! internal surface.
         IF (wall_at(i,j,k)) THEN
! ---------------------------------------------------------------->>>
! the default is equivalent to an impermeable surface and these cells
! will be treated as such in the momentum routines
            FLAG_E(i,j,k) = 0
            FLAG_N(i,j,k) = 0
            FLAG_T(i,j,k) = 0
            FLAG_E(iminus(i,j,k),j,k) = 0
            FLAG_N(i,jminus(i,j,k),k) = 0
            FLAG_T(i,j,kminus(i,j,k)) = 0

            IF (CYCLIC_AT(i,j,k)) THEN
! make the upper (E, N, T) boundary permeable
               IF (I == IMAX2) THEN
                  IF ((J/=1.AND.J/=0.) .AND. (J/=JMAX2.AND.J/=JMAX3)) THEN
                     IF (NO_K) THEN
                        IF(.NOT.WALL_AT(iminus(i,j,k),j,k)) FLAG_E(iminus(i,j,k),j,k) = 2000
                     ELSEIF ((K/=1.AND.K/=0) .AND. (K/=KMAX2.AND.K/=KMAX3)) THEN
                        IF(.NOT.WALL_AT(iminus(i,j,k),j,k)) FLAG_E(iminus(i,j,k),j,k) = 2000
                     ENDIF
                  ENDIF
               ENDIF
               IF (J == JMAX2) THEN
                  IF ((I/=1.AND.I/=0) .AND. (I/=IMAX2.AND.I/=IMAX3)) THEN
                     IF (NO_K) THEN
                        IF(.NOT.WALL_AT(i,jminus(i,j,k),k)) FLAG_N(i,jminus(i,j,k),k) = 2000
                     ELSE IF ((K/=1.AND.K/=0) .AND. (K/=KMAX2.AND.K/=KMAX3)) THEN
                        IF(.NOT.WALL_AT(i,jminus(i,j,k),k)) FLAG_N(i,jminus(i,j,k),k) = 2000
                     ENDIF
                  ENDIF
                ENDIF
               IF (K == KMAX2) THEN
                  IF ((J/=1.AND.J/=0.) .AND. (J/=JMAX2.AND.J/=JMAX3)) THEN
                     IF ((I/=1.AND.I/=0) .AND. (I/=IMAX2.AND.I/=IMAX3) .AND. &
                       .NOT.WALL_AT(i,j,kminus(i,j,k))) FLAG_T(i,j,kminus(i,j,k)) = 2000
                  ENDIF
               ENDIF

            ENDIF   ! end if cyclic_at(i,j,k)

! ----------------------------------------------------------------<<<
         ELSEIF (fluid_at(i,j,k)) THEN
! ---------------------------------------------------------------->>>

            IF ( .NOT.WALL_AT(iminus(i,j,k),j,k) .AND. FLAG_E(iminus(i,j,k),j,k)==UNDEFINED_I) &
               FLAG_E(iminus(i,j,k),j,k) = 2000 + FLAG(iminus(i,j,k),j,k)

            IF ( .NOT.WALL_AT(i,jminus(i,j,k),k) .AND. FLAG_N(i,jminus(i,j,k),k)==UNDEFINED_I) &
               FLAG_N(i,jminus(i,j,k),k) = 2000 + FLAG(i,jminus(i,j,k),k)

            IF ( .NOT.WALL_AT(i,j,kminus(i,j,k)) .AND. FLAG_T(i,j,kminus(i,j,k))==UNDEFINED_I) &
               FLAG_T(i,j,kminus(i,j,k)) = 2000 + FLAG(i,j,kminus(i,j,k))

            IF ( .NOT.WALL_AT(iplus(i,j,k),j,k) .AND. FLAG_E(i,j,k)==UNDEFINED_I) &
               FLAG_E(i,j,k) = 2000 + FLAG(iplus(i,j,k),j,k)

            IF ( .NOT.WALL_AT(i,jplus(i,j,k),k) .AND. FLAG_N(i,j,k)==UNDEFINED_I) &
               FLAG_N(i,j,k) = 2000 + FLAG(i,jplus(i,j,k),k)

            IF ( .NOT.WALL_AT(i,j,kplus(i,j,k)) .AND. FLAG_T(i,j,k)==UNDEFINED_I) &
               FLAG_T(i,j,k) = 2000 + FLAG(i,j,kplus(i,j,k))


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
      deallocate( flag_temp )
      ! call send_recv(flag_t,1)
      ! call send_recv(flag_n,1)
      ! call send_recv(flag_e,1)

      RETURN
      END SUBROUTINE SET_FLAGS1
