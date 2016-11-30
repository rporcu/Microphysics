      MODULE residual

      Use param, only: DIM_n

      INTEGER, PARAMETER :: MAX_RESID_INDEX = 8
!
      INTEGER, PARAMETER :: NRESID = 8 + DIM_N

      INTEGER, PARAMETER :: RESID_p  =  1 ! Pressure
      INTEGER, PARAMETER :: RESID_ro =  2 ! Density, volume fraction
      INTEGER, PARAMETER :: RESID_u  =  3 ! U-velocity
      INTEGER, PARAMETER :: RESID_v  =  4 ! V-velocity
      INTEGER, PARAMETER :: RESID_w  =  5 ! W-velocity
      INTEGER, PARAMETER :: RESID_t  =  6 ! Temperature
      INTEGER, PARAMETER :: RESID_th =  7 ! Granular temperature
      INTEGER, PARAMETER :: RESID_sc =  8 ! User-defined scalar
      INTEGER, PARAMETER :: RESID_ke =  9 ! K-epsilon equations
      INTEGER, PARAMETER :: RESID_x  = 10 ! Mass fraction

! Group residuals by equation
      INTEGER, PARAMETER :: HYDRO_GRP   = 1     !hydrodynamics
      INTEGER, PARAMETER :: THETA_GRP   = 2     !Granular Energy
      INTEGER, PARAMETER :: ENERGY_GRP  = 3     !Energy
      INTEGER, PARAMETER :: SCALAR_GRP  = 5     !Scalars
      INTEGER, PARAMETER :: KE_GRP      = 6     !K-Epsilon

! Prefix of Residuals string
      INTEGER, PARAMETER :: NPREFIX  = 10
      CHARACTER, PARAMETER, DIMENSION(NPREFIX) :: RESID_PREFIX = &
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

! IJK location of maximum residual
      INTEGER :: IJK_RESID(NRESID)

! sum of residuals every 5 iterations
      DOUBLE PRECISION :: SUM5_RESID

! Residual sum within a group of equations
      LOGICAL          :: GROUP_RESID
      DOUBLE PRECISION :: RESID_GRP(6)

! Residuals to be printed out
      CHARACTER(LEN=4) :: RESID_STRING(MAX_RESID_INDEX)
      CHARACTER(LEN=8) :: RESID_GRP_STRING(6)

! Indices of residuals to be printed out
      INTEGER :: RESID_INDEX(MAX_RESID_INDEX, 2)

! For checking the over-all fluid mass balance
      DOUBLE PRECISION :: accum_resid_g

      END MODULE residual
