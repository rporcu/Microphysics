MODULE CALC_D_MOD
   use cutcell, only: CARTESIAN_GRID
   use fldvar, only: EP_G
   use fun_avg, only: AVG_X, AVG_Y, AVG_Z
   use geometry, only: AYZ, AXZ, AXY
   use param1, only: ZERO, SMALL_NUMBER, ONE
   use functions, only: IP_AT_E, IP_AT_N, IP_AT_T
   use functions, only: MFLOW_AT_E, MFLOW_AT_N, MFLOW_AT_T
   CONTAINS

      double precision function epga_x(ijk)
         use functions, only: EAST_OF
         use indices, only: I_OF
         implicit none
         integer, intent(in) :: ijk
         integer :: I, IJKE
         DOUBLE PRECISION :: AREA_FACE
         IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN   !impermeable
            EPGA_X = ZERO
         ELSE
            I = I_OF(IJK)
            IJKE = EAST_OF(IJK)
            AREA_FACE = merge(ONE, AYZ(IJK), CARTESIAN_GRID)
            EPGA_X = AREA_FACE*AVG_X(EP_G(IJK),EP_G(IJKE),I)
         ENDIF
      end function epga_x

      double precision function epga_y(ijk)
         use functions, only: NORTH_OF
         use indices, only: J_OF
         implicit none
         integer, intent(in) :: ijk
         integer :: J, IJKN
         DOUBLE PRECISION :: AREA_FACE
         IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN
            EPGA_Y = ZERO
         ELSE
            J = J_OF(IJK)
            IJKN = NORTH_OF(IJK)
            AREA_FACE = merge(ONE, AXZ(IJK), CARTESIAN_GRID)
            EPGA_Y = AREA_FACE*AVG_Y(EP_G(IJK),EP_G(IJKN),J)
         ENDIF
      end function epga_y

      double precision function epga_z(ijk)
         use functions, only: TOP_OF
         use indices, only: K_OF
         implicit none
         integer, intent(in) :: ijk
         integer :: K, IJKT
         DOUBLE PRECISION :: AREA_FACE
         IF (IP_AT_T(IJK) .OR. MFLOW_AT_T(IJK)) THEN
            EPGA_Z = ZERO
         ELSE
            K = K_OF(IJK)
            IJKT = TOP_OF(IJK)
            AREA_FACE = merge(ONE, AXY(IJK), CARTESIAN_GRID)
            EPGA_Z = AREA_FACE*AVG_Z(EP_G(IJK),EP_G(IJKT),K)
         ENDIF

      end function epga_z

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_d_n                                                !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!           pressure correction -- North                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_D(D, AXIS, A_M)

! Global Variables:
!---------------------------------------------------------------------//
! Flag: Coupled DEM simulation
      use discretelement, only: DES_CONTINUUM_COUPLED
! Volume x average at momentum cell center drag for DEM/PIC
      use discretelement, only: VXF_GDS
! Pressure scale factor
      use scales, only: P_SCALE
! Flags: Impermeable surface and mass flow at north face, IJK of cell to north
      use functions, only: IP_AT_N, MFLOW_AT_N, NORTH_OF
      use functions, only: funijk
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3, imap
      USE compar, only: istart2, iend2, jstart2, jend2, kstart2, kend2
      USE compar, only: istart1, iend1, jstart1, jend1, kstart1, kend1

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and size of solids phase arrays.
      use param, only: DIMENSION_3

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Pressure correction
      DOUBLE PRECISION, INTENT(OUT) :: d(:)
! "X", "Y", or "Z"
      CHARACTER, INTENT(IN) :: axis
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN):: A_m(DIMENSION_3,-3:3)

! Local variables:
!---------------------------------------------------------------------//
! Usual Indices
      INTEGER :: I,J,K,IJK
! Temp variable for double precision values.
      DOUBLE PRECISION :: TMPdp
!......................................................................!

      abstract interface
         function epga_t (ijk)
            DOUBLE PRECISION :: epga_t
            integer, intent (in) :: ijk
         end function epga_t
      end interface

! Average volume fraction at momentum cell centers
      procedure (epga_t), pointer :: epga => null ()

      if (axis.eq.'X') then
         EPGA => EPGA_X
      else if (axis.eq.'Y') then
         EPGA => EPGA_Y
      else if (axis.eq.'Z') then
         EPGA => EPGA_Z
      endif

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

         IJK = FUNIJK(i,j,k)

         TMPdp = -A_M(IJK,0)
         IF(DES_CONTINUUM_COUPLED) TMPdp = TMPdp + VxF_gds(IJK)

         IF(abs(TMPdp) > SMALL_NUMBER) THEN
            D(IJK) = P_SCALE*EPGA(IJK)/TMPdp
         ELSE
            D(IJK) = ZERO
         ENDIF

      ENDDO
      ENDDO
      ENDDO

      RETURN
   END SUBROUTINE CALC_D
   END MODULE CALC_D_MOD
