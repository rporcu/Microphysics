!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_ARRAY_C (ARRAY,MESSAGE)                            C
!  Purpose: print out a 3D array to standard output (character)        C
!                                                                      C
!  Author: P.Nicoletti                                Date: 10-JAN-92  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: KMAX2                                         C
!  Variables modified: K                                               C
!                                                                      C
!  Local variables: IJK                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OUT_ARRAY_C(ARRAY, MESSAGE)

      USE functions, only: funijk
      USE funits, only: unit_out
      USE geometry, only: kmax2, ijkmax2
      USE in_binary_512i, only: convert_to_io_c

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

! Array to print out
      CHARACTER(LEN=4) :: ARRAY(*)

! Message to print out
      CHARACTER(LEN=*) :: MESSAGE

! Pointer into array (points to start of a k-plane)
      INTEGER :: IJK, K

      character(LEN=4),  allocatable :: array1c(:)

!-----------------------------------------------


      allocate (array1c(ijkmax2))
      call convert_to_io_c(array,array1c,ijkmax2)

      DO K = 1, KMAX2
         IJK = FUNIJK(1,1,K)

         WRITE (UNIT_OUT, 1100) MESSAGE, K
         CALL OUT_ARRAY_KC (ARRAY1C(IJK))
      END DO
 1100 FORMAT(/,1X,A,' at K = ',I4,/)

      deallocate (array1c)

      RETURN
      END SUBROUTINE OUT_ARRAY_C
