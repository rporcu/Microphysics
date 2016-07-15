!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: physprop                                                    C
!  Purpose: Common block containing physical property data             C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE physprop

! Modules
!---------------------------------------------------------------------//
      Use param, only: dim_m
!---------------------------------------------------------------------//


! Number of solids phases
      INTEGER :: MMAX

! Real number of solids phases for GHD theory
      INTEGER :: SMAX

! Particle diameters
      DOUBLE PRECISION :: D_p0(DIM_M)

! Constant solids phase densities.
      DOUBLE PRECISION :: RO_s0(DIM_M)


      END MODULE physprop
