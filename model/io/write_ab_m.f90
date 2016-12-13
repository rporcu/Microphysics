!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Write_Ab_m(A_m, b_m, IJKMAX2, M, IER)                  C                     C
!  Author: M. Syamlal                                 Date: 16-MAY-96  C
!                                                                      C
!  Purpose: Write the sparse matrix coefficients and the               C
!           source vector.                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_AB_M(A_M, B_M, M)    ! pnicol

      USE functions, only: funijk
      USE compar, only: iend3, jend3, kend3
      USE compar, only: istart3, jstart3, kstart3
      USE compar, only: mype, pe_io
      USE funits, only: dmp_log, unit_log
      USE geometry, only: imax2, jmax2, kmax2
      USE geometry, only: imin2, jmin2, kmin2

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!                      Phase index
      INTEGER          M
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! Local index
      INTEGER          L
!-----------------------------------------------
!
      integer i, j, k


      if (myPE == PE_IO) then
         IF(DMP_LOG)WRITE (UNIT_LOG,*) ' Note : write_am_m is VERY inefficient '
         IF(DMP_LOG)WRITE (UNIT_LOG,*) '  '
         IF(DMP_LOG)WRITE (UNIT_LOG,*) ' A_m and B_m arrays below are in the '
         IF(DMP_LOG)WRITE (UNIT_LOG,*) ' mfix INTERNAL order'
         IF(DMP_LOG)WRITE (UNIT_LOG,*) ' '
         IF(DMP_LOG)WRITE (UNIT_LOG, '(A,A)') &
           '  IJK  I  J  K   b         s         w         p         e       ', &
           '  n         t        Source'
      end if



      DO K = Kmin2, Kmax2
      DO I = Imin2, Imax2
      DO J = Jmin2, Jmax2


      if (myPE == PE_IO .AND. DMP_LOG) &
         WRITE (UNIT_LOG, '(I5, 3(I3), 8(1X,G9.2))') &
         FUNIJK(I,J,K), I, J, K,(A_M(i,j,k,L),L=-3,3), B_M(i,j,k)

      END DO
      END DO
      END DO


      RETURN
      END SUBROUTINE WRITE_AB_M
