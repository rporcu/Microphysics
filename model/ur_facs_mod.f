      MODULE ur_facs

      use param, only: DIM_EQS

! Under relaxation factors
      DOUBLE PRECISION :: UR_FAC(DIM_EQS)

! Under relaxation factors for coefficient update:
!  [0]  every time step (explicit)
!  [1]  every iteration (implicit)
! (0,1) underrelaxed

! Note that these values need to be temporarily set to 1 before the calc_coeff
! call in time_march.  And after the call reset to their original value.


      END MODULE ur_facs
