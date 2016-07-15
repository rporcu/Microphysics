!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: mflux                                                  C
!  Purpose: Module for mass fluxes and densities at faces              C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE mflux

! x-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gE
! y-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gN
! z-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gT

! macroscopic gas density at east face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gE
! macroscopic gas density at north face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gN
! macroscopic gas density at top face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gT


      END MODULE mflux
