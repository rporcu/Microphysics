MODULE CALC_D_MOD
   use fldvar, only: EP_G
   use fun_avg, only: AVG_X, AVG_Y, AVG_Z
   use geometry, only: AYZ, AXZ, AXY
   use param1, only: ZERO, SMALL_NUMBER, ONE
   use functions, only: ip_at_e, ip_at_n, ip_at_t
   use functions, only: MFLOW_AT_E, MFLOW_AT_N, MFLOW_AT_T
   CONTAINS

      double precision function epga_x(i,j,k)
         use functions, only: FUNIJK, ieast
         implicit none
         integer, intent(in) :: i,j,k
         integer :: IJK,IJKE
         DOUBLE PRECISION :: AREA_FACE

         IJK = FUNIJK(i,j,k)

         IF (ip_at_e(i,j,k) .OR. MFLOW_AT_E(IJK)) THEN   !impermeable
            EPGA_X = ZERO
         ELSE
            IJKE = FUNIJK(ieast(i,j,k),j,k)
            AREA_FACE = AYZ(IJK)
            EPGA_X = AREA_FACE*AVG_X(EP_G(IJK),EP_G(IJKE),I)
         ENDIF
      end function epga_x

      double precision function epga_y(i,j,k)
         use functions, only: FUNIJK, jnorth
         implicit none
         integer, intent(in) :: i,j,k
         integer :: IJK, IJKN
         DOUBLE PRECISION :: AREA_FACE

         IJK = FUNIJK(i,j,k)

         IF (ip_at_n(i,j,k) .OR. MFLOW_AT_N(IJK)) THEN
            EPGA_Y = ZERO
         ELSE
            IJKN = FUNIJK(i,jnorth(i,j,k),k)
            AREA_FACE = AXZ(IJK)
            EPGA_Y = AREA_FACE*AVG_Y(EP_G(IJK),EP_G(IJKN),J)
         ENDIF
      end function epga_y

      double precision function epga_z(i,j,k)
         use functions, only: FUNIJK, ktop
         implicit none
         integer, intent(in) :: i,j,k
         integer :: IJK, IJKT
         DOUBLE PRECISION :: AREA_FACE

         IJK = FUNIJK(i,j,k)
         IF (ip_at_t(i,j,k) .OR. MFLOW_AT_T(IJK)) THEN
            EPGA_Z = ZERO
         ELSE
            IJKT = FUNIJK(i,j,ktop(i,j,k))
            AREA_FACE = AXY(IJK)
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
         function epga_t (i,j,k)
            DOUBLE PRECISION :: epga_t
            integer, intent (in) :: i,j,k
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
            D(IJK) = P_SCALE*EPGA(i,j,k)/TMPdp
         ELSE
            D(IJK) = ZERO
         ENDIF

      ENDDO
      ENDDO
      ENDDO

      RETURN
   END SUBROUTINE CALC_D
   END MODULE CALC_D_MOD
