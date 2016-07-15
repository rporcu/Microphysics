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

! y-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sN
! x-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sE
! z-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sT

! macroscopic gas density at east face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gE
! macroscopic gas density at north face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gN
! macroscopic gas density at top face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gT
! macroscopic solids density at north face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sN
! macroscopic solids density at east face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sE
! macroscopic solids density at top face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sT


! for GHD Theory
! x-component of solids total number density flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_nE
! y-component of solids total number density flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_nN
! z-component of solids total number density flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_nT
! end GHD Theory modification


      END MODULE mflux
