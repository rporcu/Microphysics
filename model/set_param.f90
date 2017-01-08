MODULE SET_PARAMETERS_MODULE
CONTAINS
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
      use constant, only: MMAX

! Domain indices.
      use geometry, only: IMAX2, JMAX2, KMAX2

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of species.
      use param, only: DIMENSION_N_g ! Gas
      use param, only: DIMENSION_N_s ! Solids
! Total number of species
      use param1, only: DIMENSION_N_all
! Axis decomposition
      USE param, only: DIMENSION_3

      IMPLICIT NONE
!......................................................................!

! Number of gas phase species.
      DIMENSION_N_g = 1

! Number of solids phase species. (Max over all phases)
      DIMENSION_N_s = 1

! Max number of species over all phases
      DIMENSION_N_all = max(DIMENSION_N_g, DIMENSION_N_s)

      DIMENSION_3 = IMAX2*JMAX2*KMAX2

      RETURN
      END SUBROUTINE SET_PARAMETERS
END MODULE SET_PARAMETERS_MODULE
