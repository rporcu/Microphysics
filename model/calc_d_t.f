!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_D_T                                                !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: Calculate coefficients linking velocity correction to      !
!           pressure correction -- Top                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_D_T(A_M, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Number of solids phases.
      use physprop, only: MMAX
! Flag: Solve the Z-Momentum Equations
      use run, only: MOMENTUM_Z_EQ
! Flag: Coupled DEM, PIC, or Hybrid simulation
      use discretelement, only: DES_CONTINUUM_COUPLED
! Pressure correction coefficient
      use fldvar, only: D_T
! Volume x average at momentum cell center drag for DEM/PIC
      use discretelement, only: VXF_GDS
! Volume fractions of gas and solids phases.
     use fldvar, only: EP_G
! Pressure scale factor
     use scales, only: P_SCALE
! Flag: Cartesian grid simulation
      use cutcell, only: CARTESIAN_GRID
! Area of cell's XY face.
      use geometry, only: AXY
! Volume of W-momentum cell.
      use geometry, only: VOL_W
! Function to average across Z face.
      use fun_avg, only: AVG_Z
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Flags: Impermeable surface and mass flow at top face, IJK of cell to top
      use functions, only: IP_AT_T, MFLOW_AT_T, TOP_OF
! Indices: K index of cell
      use indices, only: K_OF

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and size of solids phase arrays.
      use param, only: DIMENSION_3, DIMENSION_M
! Double precision parameters.
      use param1, only: ZERO, SMALL_NUMBER, ONE


      IMPLICIT NONE


! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN):: A_m(DIMENSION_3,-3:3,0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local variables:
!---------------------------------------------------------------------//
! Usual Indices
      INTEGER :: IJK, IJKT, K
! Average volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPGA
! local alias for face area
      DOUBLE PRECISION :: AREA_FACE
! Temp variable for double precision values.
      DOUBLE PRECISION :: TMPdp
!......................................................................!

! Initialize the error flag.
      IER = 0

! Determine which calculations are needed
      IF (MOMENTUM_Z_EQ(0)) THEN

         DO IJK = ijkstart3, ijkend3
            IF (IP_AT_T(IJK) .OR. MFLOW_AT_T(IJK)) THEN
               D_T(IJK) = ZERO
               CYCLE
            ENDIF
 
            AREA_FACE = merge(ONE, AXY(IJK), CARTESIAN_GRID)
            K = K_OF(IJK)
            IJKT = TOP_OF(IJK)
            EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)
 
            TMPdp = -A_M(IJK,0,0)
            IF(DES_CONTINUUM_COUPLED) TMPdp = TMPdp + VXF_GDS(IJK)

            IF(abs(TMPdp) > SMALL_NUMBER) THEN
               D_T(IJK) = P_SCALE*AREA_FACE*EPGA/TMPdp
            ELSE
               D_T(IJK) = ZERO
            ENDIF
         ENDDO

      ENDIF

      RETURN
      END SUBROUTINE CALC_D_T
