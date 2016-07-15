!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_d_e                                                !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!  pressure correction -- East                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_D_E(A_M, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve the X-Momentum Equations
      use run, only: MOMENTUM_X_EQ
! Flag: Coupled DEM, PIC, or Hybrid simulation
      use discretelement, only: DES_CONTINUUM_COUPLED
! Volume x average at momentum cell center drag for DEM/PIC
      use discretelement, only: VXF_GDS
! Pressure correction coefficient
      use fldvar, only: D_E
! Volume fraction of gas phase.
      use fldvar, only: EP_G
! Number of solids phases.
      use physprop, only: MMAX
! Pressure scale factor
      use scales, only: P_SCALE
! Flag: Cartesian grid simulation
      use cutcell, only: CARTESIAN_GRID
! Area of cell's YZ face.
      use geometry, only: AYZ
! Volume of V-momentum cell.
      use geometry, only: VOL_U
! Function to average across Y face.
      use fun_avg, only: AVG_X
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Flags: Impermeable surface and mass flow at east face, IJK of cell to east
      use functions, only: IP_AT_E, MFLOW_AT_E, EAST_OF
! Indices: J index of cell
      use indices, only: I_OF

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and size of solids phase arrays.
      use param, only: DIMENSION_3, DIMENSION_M
! Double precision parameters.
      use param1, only: ZERO, SMALL_NUMBER, ONE


      IMPLICIT NONE


! Dummy arguments
!---------------------------------------------------------------------//
! Error index
      INTEGER, INTENT(INOUT) :: IER
! Septadiagonal matrix A_m.  The center coefficient is negative.
      DOUBLE PRECISION, INTENT(IN):: A_m(DIMENSION_3,-3:3,0:DIMENSION_M)

! Local variables:
!---------------------------------------------------------------------//
! Usual Indices
      INTEGER :: I, IJK, IJKE
! Average volume fraction at momentum cell centers
      DOUBLE PRECISION :: EPGA
! local alias for face area
      DOUBLE PRECISION :: AREA_FACE
! Temp variable for double precision values.
      DOUBLE PRECISION :: TMPdp

!......................................................................!

! Initialize the error flag.
      IER = 0


      IF(MOMENTUM_X_EQ(0)) THEN

         DO IJK = IJKSTART3, IJKEND3
 
            IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN   !impermeable
               D_E(IJK) = ZERO
               CYCLE
            ENDIF
 
            AREA_FACE = merge(ONE, AYZ(IJK), CARTESIAN_GRID)
 
            I = I_OF(IJK)
            IJKE = EAST_OF(IJK)
            EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

! Add DEM temp A_M so they are included in the pressure correciton eq.
            TMPdp = -A_M(IJK,0,0)
            IF(DES_CONTINUUM_COUPLED) TMPdp = TMPdp + VXF_GDS(IJK)

            IF(abs(TMPdp) > SMALL_NUMBER) THEN
               D_E(IJK) = P_SCALE*AREA_FACE*EPGA/TMPdp

            ELSE
               D_E(IJK) = ZERO
            ENDIF
         ENDDO

      ENDIF

      RETURN
      END SUBROUTINE CALC_D_E
