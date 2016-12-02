      MODULE residual

      Use param, only: DIM_n

      integer, parameter :: MAX_RESID_INDEX = 8
!
      integer, parameter :: NRESID = 8 + DIM_N

      integer, parameter :: RESID_p  =  1 ! Pressure
      integer, parameter :: RESID_ro =  2 ! Density, volume fraction
      integer, parameter :: RESID_u  =  3 ! U-velocity
      integer, parameter :: RESID_v  =  4 ! V-velocity
      integer, parameter :: RESID_w  =  5 ! W-velocity
      integer, parameter :: RESID_t  =  6 ! Temperature
      integer, parameter :: RESID_th =  7 ! Granular temperature
      integer, parameter :: RESID_sc =  8 ! User-defined scalar
      integer, parameter :: RESID_ke =  9 ! K-epsilon equations
      integer, parameter :: RESID_x  = 10 ! Mass fraction

! Group residuals by equation
      integer, parameter :: HYDRO_GRP   = 1     !hydrodynamics
      integer, parameter :: THETA_GRP   = 2     !Granular Energy
      integer, parameter :: ENERGY_GRP  = 3     !Energy
      integer, parameter :: SCALAR_GRP  = 5     !Scalars
      integer, parameter :: KE_GRP      = 6     !K-Epsilon

! Prefix of Residuals string
      integer, parameter :: NPREFIX  = 10
      CHARACTER, parameter, DIMENSION(NPREFIX) :: RESID_PREFIX = &
        (/ 'P', 'R', 'U', 'V', 'W', 'T', 'G', 'S', 'K', 'X' /)

! Average residual
      DOUBLE PRECISION :: RESID(NRESID)
! Maximum residual
      DOUBLE PRECISION :: MAX_RESID(NRESID)
! Residual Numerator
      DOUBLE PRECISION :: NUM_RESID(NRESID)
! Residual Denominator
      DOUBLE PRECISION :: DEN_RESID(NRESID)
! Residual Packing for Global Operations
      DOUBLE PRECISION :: RESID_PACK(2*NRESID)

! (i,j,k) location of maximum residual
      integer :: i_resid(nresid)
      integer :: j_resid(nresid)
      integer :: k_resid(nresid)

! sum of residuals every 5 iterations
      DOUBLE PRECISION :: SUM5_RESID

! Residual sum within a group of equations
      LOGICAL          :: GROUP_RESID
      DOUBLE PRECISION :: RESID_GRP(6)

! Residuals to be printed out
      CHARACTER(LEN=4) :: RESID_STRING(MAX_RESID_INDEX)
      CHARACTER(LEN=8) :: RESID_GRP_STRING(6)

! Indices of residuals to be printed out
      integer :: RESID_INDEX(MAX_RESID_INDEX, 2)

! For checking the over-all fluid mass balance
      DOUBLE PRECISION :: accum_resid_g

      END MODULE residual
