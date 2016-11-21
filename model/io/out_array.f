!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_ARRAY (ARRAY, MESSAGE)                             C
!  Author: P.Nicoletti                                Date: 02-DEC-91  C
!                                                                      C
!  Purpose: print out a 3D array to standard output                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_ARRAY(ARRAY, MESSAGE)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE fldvar
      USE physprop
      USE funits
      USE compar
      USE functions
      USE in_binary_512

      IMPLICIT NONE

!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! Array to print out
      DOUBLE PRECISION :: ARRAY(*)

! Message to print out
      CHARACTER(LEN=*) :: MESSAGE

! local variables
!
! pointer into array (points to start of a k-plane)
      INTEGER :: IJK
! loop counter
      INTEGER :: K

      double precision,  allocatable :: array1(:)
!
!-----------------------------------------------

      allocate (array1(ijkmax2))

      call convert_to_io_dp(array,array1,ijkmax2)


      IF(CYCLIC_Z) then
        DO K = 2, KMAX1
           IJK = FUNIJK_IO(1,1,K)
           WRITE (UNIT_OUT, 1100) MESSAGE, K
           CALL OUT_ARRAY_K (ARRAY1(IJK))
        END DO
      ELSE
        DO K = 1, KMAX2
           IJK = FUNIJK_IO(1,1,K)
           WRITE (UNIT_OUT, 1100) MESSAGE, K
           CALL OUT_ARRAY_K (ARRAY1(IJK))
        END DO
      ENDIF
 1100 FORMAT(/,1X,A,' at K = ',I4,/)


      deallocate (array1)

      RETURN
      END SUBROUTINE OUT_ARRAY
