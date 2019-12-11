!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: write_out0                                              !
!  Author: P. Nicoletti, M. Syamlal                   Date: 04-DEC-91  !
!                                                                      !
!  Purpose: Echo user input.                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine write_out0 (time, dx, dy, dz, xlength, ylength, zlength, domlo, domhi) &
         bind(C, name="write_out0")

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use discretelement, only: kn, kt, kn_w, kt_w
      use discretelement, only: des_etan, des_etat, des_etat_wall, des_etan_wall

      use param, only: dim_ic
      use param, only: half, undefined, zero, is_defined
      use constant, only: mmax

      use ic, only: write_out_ic

      implicit none

      integer(c_int), intent(in   ) :: domlo(3), domhi(3)
      real(rt)  , intent(in   ) :: time, dx, dy, dz
      real(rt)  , intent(in   ) :: xlength, ylength, zlength

      integer :: M, N
      integer :: mmax_tot

      character(LEN= 3), DIMENSION(  3) :: legend

      integer, PARAMETER :: unit_out = 52
      integer :: ier

      !-----------------------------------------------
      !
      mmax_tot = mmax

      open(unit=unit_out, file='MFIX-EXA.OUT', status='unknown', &
         access='sequential', form='formatted', position='append', iostat=ier)

!  Write Headers for .OUT file
!
      write(unit_out,1000)
!
      RETURN
 1000 FORMAT(17X,'MM      MM  FFFFFFFFFF    IIIIII    XX      XX',/17X,&
         'MM      MM  FFFFFFFFFF    IIIIII    XX      XX',/17X,&
         'MMMM  MMMM  FF              II      XX      XX',/17X,&
         'MMMM  MMMM  FF              II      XX      XX',/17X,&
         'MM  MM  MM  FF              II        XX  XX  ',/17X,&
         'MM  MM  MM  FF              II        XX  XX  ',/17X,&
         'MM      MM  FFFFFFFF        II          XX    ',/17X,&
         'MM      MM  FFFFFFFF        II          XX    ',/17X,&
         'MM      MM  FF              II        XX  XX  ',/17X,&
         'MM      MM  FF              II        XX  XX  ',/17X,&
         'MM      MM  FF              II      XX      XX',/17X,&
         'MM      MM  FF              II      XX      XX',/17X,&
         'MM      MM  FF            IIIIII    XX      XX',/17X,&
         'MM      MM  FF            IIIIII    XX      XX',2/20X,&
         'Multiphase Flow with Interphase eXchanges')

!
    contains


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: write_table (legend, ARRAY, DIST_MIN, LSTART, LEND)    C
!  Purpose: To write a table of DX, DY, DZ, and cell wall locations    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine write_table (legend, SCALAR, DIST_MIN, LSTART, LEND)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

!                      Legend
      CHARACTER(LEN=*)    legend(3)

!                      DX, DY, or DZ Array to be written


!                      Starting array index
      integer          LSTART

!                      Ending array index
      integer          LEND

      real(rt) SCALAR

!                      Starting value of distance
      real(rt) DIST_MIN

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------

!                      Number of columns in the table.  When this is changed
!                      remember to change the FORMAT statement also.
      integer, PARAMETER :: NCOL = 5
!
!                      Some dimension large enough for I, J, and K.
      integer, PARAMETER :: DIMENSION_1 = 5000

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      integer          ARRAY1(DIMENSION_1)
!
!                      Array3 to be written
      real(rt) ARRAY3(DIMENSION_1)
!
!                      Number of rows
      integer          NROW
!
!                      Temporary storage for distance calculation
      real(rt) DIST
!
!                      Local array indices
      integer          L, L1, L2, L3
!-----------------------------------------------
!
!
!  Fill arrays 1 and 3
!
      DIST = DIST_MIN
      do L = LSTART, LEND
         ARRAY1(L) = L
         ARRAY3(L) = DIST
         IF (L < LEND) DIST = DIST + SCALAR
      end do
      NROW = (LEND - LSTART + 1)/NCOL
!
      L2 = LSTART - 1
      do L = 1, NROW
         L1 = L2 + 1
         L2 = L1 + NCOL - 1
         write (unit_out, 1010) legend(1), (ARRAY1(L3),L3=L1,L2)
         write (unit_out, 1020) legend(2), (SCALAR,L3=L1,L2)
         write (unit_out, 1030) legend(3), (ARRAY3(L3),L3=L1,L2)
      end do
      if (NROW*NCOL < LEND - LSTART + 1) THEN
         L1 = L2 + 1
         L2 = LEND
         write (unit_out, 1010) legend(1), (ARRAY1(L3),L3=L1,L2)
         write (unit_out, 1020) legend(2), (SCALAR,L3=L1,L2)
         write (unit_out, 1030) legend(3), (ARRAY3(L3),L3=L1,L2)
      end if

      return
!
 1010 FORMAT(7X,A3,2X,5(4X,I3,5X,1X))
 1020 FORMAT(7X,A3,2X,5(G12.5,1X))
 1030 FORMAT(7X,A3,2X,5(G12.5,1X),/)
      end subroutine write_table

      end subroutine write_out0
