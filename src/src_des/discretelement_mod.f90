!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!   Module name: DISCRETELEMENT                                        C
!   Purpose: DES mod file                                              C
!            Common Block containing DEM conditions                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE DISCRETELEMENT

!-----------------------------------------------
! Modules
!-----------------------------------------------

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use param, only: dim_m
      use param, only: undefined, one, zero, half
      implicit none
!-----------------------------------------------

! Total number of particles in simulation: read from input or generated
!      integer :: PARTICLES

! End Output/debug controls
!-----------------------------------------------------------------<<<

! Value of solids time step based on particle properties
      real(rt) :: dtsolid

! End neighbor search related quantities
!-----------------------------------------------------------------<<<


! Particle-particle and Particle-wall collision model parameters
!----------------------------------------------------------------->>>
! Spring contants
      real(rt) :: KN, KN_W  !Normal
      real(rt) :: KT, KT_W, KT_FAC, KT_W_FAC

! Damping coefficients
      real(rt) :: ETA_DES_N, ETA_N_W  !Normal
      real(rt) :: ETA_DES_T, ETA_T_W  !Tangential

! Tangential damping factors, eta_t = eta_t_factor * eta_N
      real(rt) :: DES_ETAT_FAC, DES_ETAT_W_FAC

! Damping coeffients in array form
      real(rt) :: DES_ETAN(DIM_M, DIM_M), DES_ETAN_WALL(DIM_M)
      real(rt) :: DES_ETAT(DIM_M, DIM_M), DES_ETAT_WALL(DIM_M)

! Friction coefficients
      real(rt) :: MEW, MEW_W

! coeff of restituion input in one D array, solid solid
      real(rt) :: DES_EN_INPUT(DIM_M+DIM_M*(DIM_M-1)/2)

! coeff of restitution input in one D array, solid wall
      real(rt) :: DES_EN_WALL_INPUT(DIM_M)

      real(rt) :: dp_max(dim_m) = -1.0d0
      real(rt) :: dp_min(dim_m) = -1.0d0
      real(rt) :: dp_avg(dim_m) = -1.0d0

      real(rt) :: ro_max(dim_m) = -1.0d0
      real(rt) :: ro_min(dim_m) = -1.0d0
      real(rt) :: ro_avg(dim_m) = -1.0d0

! End particle-particle and particle-wall collision model parameters
!-----------------------------------------------------------------<<<

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                     !
!  Subroutine: DES_CROSSPRDCT                                         !
!  Purpose: Calculate the cross product of two 3D vectors.            !
!                                                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      FUNCTION DES_CROSSPRDCT (XX,YY)

      implicit none

! Dummy arguments
!---------------------------------------------------------------------//
! Input vectors
      real(rt), DIMENSION(3), INTENT(IN) :: XX, YY
! Result: cross product of vectors
      real(rt), DIMENSION(3) :: DES_CROSSPRDCT
!......................................................................!

      DES_CROSSPRDCT(1) = XX(2)*YY(3) - XX(3)*YY(2)
      DES_CROSSPRDCT(2) = XX(3)*YY(1) - XX(1)*YY(3)
      DES_CROSSPRDCT(3) = XX(1)*YY(2) - XX(2)*YY(1)

      END FUNCTION DES_CROSSPRDCT

      END MODULE DISCRETELEMENT
