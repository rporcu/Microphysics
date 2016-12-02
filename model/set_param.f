!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutinee: SET_PARAMETERS                                         !
!  Purpose: Set parameters used in array allocations.                  !
!                                                                      !
!  Author: J.Musser                                  Date: 17-APR-14   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_PARAMETERS

! Global Varialbes:
!---------------------------------------------------------------------//
! Number of solids
      use physprop, only: MMAX

! Domain indices.
      use geometry, only: IMAX3, JMAX3, KMAX3, IJKMAX3

! Rank specific decompositions.
      USE compar, only: iStart3, iEnd3, iStart4, iEnd4
      USE compar, only: jStart3, jEnd3, jStart4, jEnd4
      USE compar, only: kStart3, kEnd3, kStart4, kEnd4

! Global Parameters:
!---------------------------------------------------------------------//
! Total number of solids phases
      use param, only: DIMENSION_M
! Maximum number of species.
      use param, only: DIMENSION_N_g ! Gas
      use param, only: DIMENSION_N_s ! Solids
! Dimension for the upper triangle of an MxM matrix
      use param1, only: DIMENSION_LM
! Total number of species
      use param1, only: DIMENSION_N_all
! Parameter constants.
      use param1, only: UNDEFINED_I
! Axis decomposition
      USE param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
      USE param, only: DIMENSION_3

! MPI-Domain decompoint and rank flags.
      use compar, only: myPE

      IMPLICIT NONE
!......................................................................!

! The total number of solids.
      DIMENSION_M = MMAX

! Number of gas phase species.
      DIMENSION_N_g = 1

! Number of solids phase species. (Max over all phases)
      DIMENSION_N_s = 1

! Max number of species over all phases
      DIMENSION_N_all = max(DIMENSION_N_g, DIMENSION_N_s)

! Size of MxM upper triangular matrix.
      DIMENSION_LM = (DIMENSION_M * (DIMENSION_M-1)/2)+1

! Set DIMENSION_x variables.
      DIMENSION_I = IMAX3
      DIMENSION_J = JMAX3
      DIMENSION_K = KMAX3

      DIMENSION_3 = (kEnd3-kStart3+1)*(jEnd3-jStart3+1)*(iEnd3-iStart3+1)

      RETURN
      END SUBROUTINE SET_PARAMETERS
