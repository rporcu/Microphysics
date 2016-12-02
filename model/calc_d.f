!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_d_n                                                !
!  Author: M. Syamlal                                 Date: 21-JUN-96  !
!                                                                      !
!  Purpose: calculate coefficients linking velocity correction to      !
!           pressure correction -- North                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
MODULE CALC_D_MOD

   CONTAINS

      SUBROUTINE CALC_D(D, AXIS, A_M, ep_g)

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

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and size of solids phase arrays.

      use compar, only: istart2,iend2,jstart2,jend2,kstart2,kend2
      use compar, only: istart3,iend3,jstart3,jend3,kstart3,kend3

      use fun_avg, only: AVG
      use geometry, only: AYZ, AXZ, AXY
      use param1, only: ZERO, SMALL_NUMBER
      use functions, only: ip_at_e, ip_at_n, ip_at_t
      use functions, only: ieast, jnorth, ktop
      use functions, only: MFLOW_AT_E, MFLOW_AT_N, MFLOW_AT_T

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Pressure correction
      DOUBLE PRECISION, INTENT(OUT) :: d(:,:,:)
! "X", "Y", or "Z"
      CHARACTER, INTENT(IN) :: axis
! Septadiagonal matrix A_M
      DOUBLE PRECISION, INTENT(IN):: A_M&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)

      DOUBLE PRECISION, INTENT(IN):: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

! Local variables:
!---------------------------------------------------------------------//
! Usual Indices
      INTEGER :: I,J,K

! Temp variable for double precision values.
      DOUBLE PRECISION :: TMPdp
      DOUBLE PRECISION :: EPGA
!......................................................................!

      DO K = kstart2, kend2
        DO J = jstart2, jend2
          DO I = istart2, iend2

         TMPdp = -A_M(I,J,K,0)

         IF(DES_CONTINUUM_COUPLED) TMPdp = TMPdp + VxF_gds(i,j,k)

         IF(abs(TMPdp) > SMALL_NUMBER) THEN

            if (axis.eq.'X') then

               IF (ip_at_e(i,j,k) .OR. MFLOW_AT_E(i,j,k)) THEN   !impermeable
                  EPGA = ZERO
               ELSE
                  EPGA = AYZ*AVG(EP_G(I,J,K),EP_G(ieast(i,j,k),j,k))
               ENDIF

            else if (axis.eq.'Y') then

               IF (ip_at_n(i,j,k) .OR. MFLOW_AT_N(i,j,k)) THEN
                  EPGA = ZERO
               ELSE
                  EPGA = AXZ*AVG(EP_G(I,J,K),EP_G(i,jnorth(i,j,k),k))
               ENDIF

            else if (axis.eq.'Z') then
               IF (ip_at_t(i,j,k) .OR. MFLOW_AT_T(i,j,k)) THEN
                  EPGA = ZERO
               ELSE
                  EPGA = AXY*AVG(EP_G(I,J,K),EP_G(i,j,ktop(i,j,k)))
               ENDIF
            endif

            D(I,J,K) = P_SCALE*EPGA/TMPdp

         ELSE
            D(I,J,K) = ZERO
         ENDIF

      ENDDO
      ENDDO
      ENDDO

   END SUBROUTINE CALC_D

   END MODULE CALC_D_MOD
