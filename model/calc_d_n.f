!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_d_n                                                !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!           pressure correction -- North                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_D_N(A_M, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve the Y-Momentum Equations
      use run, only: MOMENTUM_Y_EQ
! Flag: Coupled DEM, PIC, or Hybrid simulation
      use discretelement, only: DES_CONTINUUM_COUPLED
! Volume x average at momentum cell center drag for DEM/PIC
      use discretelement, only: VXF_GDS
! Pressure correction coefficient
      use fldvar, only: D_N
! Volume fractions of gas and solids phases.
      use fldvar, only: EP_G
! Pressure scale factor
      use scales, only: P_SCALE
! Flag: Cartesian grid simulation
      use cutcell, only: CARTESIAN_GRID
! Area of cell's XZ face.
      use geometry, only: AXZ
! Volume of V-momentum cell.
      use geometry, only: VOL_V
! Function to average across Y face.
      use fun_avg, only: AVG_Y
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Flags: Impermeable surface and mass flow at north face, IJK of cell to north
      use functions, only: IP_AT_N, MFLOW_AT_N, NORTH_OF
! Indices: J index of cell
      use indices, only: J_OF

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
! Integer error flag.
      INTEGER, INTENT(OUT) :: IER

! Local variables:
!---------------------------------------------------------------------//
! Usual Indices
      INTEGER :: J, IJK, IJKN
! Average volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPGA
! local alias for face area
      DOUBLE PRECISION :: AREA_FACE
! Temp variable for double precision values.
      DOUBLE PRECISION :: TMPdp
!......................................................................!

! Initialize the error flag.
      IER = 0

      IF(MOMENTUM_Y_EQ(0))THEN

         DO IJK = IJKSTART3, IJKEND3
 
            IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN
               D_N(IJK) = ZERO
               CYCLE
            ENDIF
 
            AREA_FACE = merge(ONE, AXZ(IJK), CARTESIAN_GRID)
            J = J_OF(IJK)
            IJKN = NORTH_OF(IJK)
            EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)

            TMPdp = -A_M(IJK,0,0)
            IF(DES_CONTINUUM_COUPLED) TMPdp = TMPdp + VxF_gds(IJK)

            IF(abs(TMPdp) > SMALL_NUMBER) THEN
               D_N(IJK) = P_SCALE*AREA_FACE*EPGA/TMPdp
            ELSE
               D_N(IJK) = ZERO
            ENDIF
 
         ENDDO

      ENDIF

      RETURN
      END SUBROUTINE CALC_D_N

