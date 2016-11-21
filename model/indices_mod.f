!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: indices                                                     !
!  Author: M. Syamlal                                 Date: dd-mmm-yy  !
!                                                                      !
!  Purpose: Global arrays for index computations.                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE indices

! +/- increments, shifted for cyclic BCs and DMP partitions.
      INTEGER, DIMENSION(:), ALLOCATABLE :: Im1, Ip1 ! I-1, I+1
      INTEGER, DIMENSION(:), ALLOCATABLE :: Jm1, Jp1 ! J-1, J+1
      INTEGER, DIMENSION(:), ALLOCATABLE :: Km1, Kp1 ! K-1, K+1

      END MODULE indices
