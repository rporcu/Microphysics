!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Write_Ab_m_var(A_m, b_m, var, IER)                     C
!  Purpose: Write the sparse matrix coefficients and the               C
!           source vector.                                             C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_AB_M_VAR(A_M, B_M, VAR)

         USE compar, only: mype, pe_io
         USE compar, only: istart3, jstart3, kstart3
         USE compar, only: iend3, jend3, kend3
         USE geometry, only: imin2, jmin2, kmin2
         USE geometry, only: imax2, jmax2, kmax2
         USE functions, only: funijk, imap_c, jmap_c, kmap_c
         USE funits, only: dmp_log, unit_log
         USE param, only: DIMENSION_3

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Local index
      INTEGER          L
!
!                      cell index
      INTEGER          IJK
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
!
!                      Source vector
      DOUBLE PRECISION :: b_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

!                      Source vector
      DOUBLE PRECISION var(DIMENSION_3)

!-----------------------------------------------
!
      integer i, j, k


      if (myPE == PE_IO) then
         IF(DMP_LOG)WRITE (UNIT_LOG,*) ' Note : write_am_m is VERY inefficient '
         WRITE (*,*) ' Note : write_am_m is VERY inefficient '
         IF(DMP_LOG)WRITE (UNIT_LOG,*) '  '
         IF(DMP_LOG)WRITE (UNIT_LOG,*) ' A_m and B_m arrays below are in the '
         IF(DMP_LOG)WRITE (UNIT_LOG,*) ' mfix INTERNAL order'
         IF(DMP_LOG)WRITE (UNIT_LOG,*) ' '
         IF(DMP_LOG)WRITE (UNIT_LOG, '(A,A)') &
           '  IJK  I  J  K   b         s         w         p         e       ', &
           '  n         t        Source     Variable'
      end if




      DO K = Kmin2, Kmax2
      DO I = Imin2, Imax2
      DO J = Jmin2, Jmax2


      IJK = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K))

      if (myPE == PE_IO .AND. DMP_LOG)&
         WRITE (UNIT_LOG, '(I5, 3(I3), 9(1X,G9.2))') &
         FUNIJK(I,J,K), I, J, K,&
         (A_M(i,j,k,L),L=-3,3), b_m(I,J,K), var(IJK)

      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE WRITE_AB_M_VAR
