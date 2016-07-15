!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: VF_gs_X                                                 !
!  Author: M. Syamlal                                 Date: 20-MAY-96  !
!                                                                      !
!  Purpose: Calculate the average drag coefficient at i+1/2, j, k and  !
!           multiply with u-momentum cell volume.                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE VF_GS_X

      USE compar
      USE discretelement
      USE drag
      USE fun_avg
      USE functions
      USE geometry
      USE indices

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, IJK, IJKE

! Calculate the combined effect for all discrete solids.
      IF(DISCRETE_ELEMENT .AND. .NOT.DES_ONEWAY_COUPLED) THEN
         DO IJK = IJKSTART3, IJKEND3
            IF(IP_AT_E(IJK)) THEN
               VXF_GDS(IJK) = ZERO
            ELSE
               I = I_OF(IJK)
               IJKE = EAST_OF(IJK)
               VXF_GDS(IJK) = VOL_U(IJK) *                             &
                  AVG_X(F_GDS(IJK),F_GDS(IJKE),I)
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE VF_GS_X

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: VF_gs_Y                                                 C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
!                                                                      C
!  Purpose: Calculate the average drag coefficient at i, j+1/2, k and  C
!           multiply with v-momentum cell volume.                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE VF_GS_Y

      USE geometry
      USE indices
      USE compar
      USE drag
      USE discretelement
      USE fun_avg
      USE functions

      IMPLICIT NONE

! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: J, IJK, IJKN
!-----------------------------------------------

      IF(DISCRETE_ELEMENT .AND. .NOT.DES_ONEWAY_COUPLED) THEN
         DO IJK = IJKSTART3, IJKEND3
            IF(IP_AT_N(IJK)) THEN
               VXF_GDS(IJK) = ZERO
            ELSE
               J = J_OF(IJK)
               IJKN = NORTH_OF(IJK)
               VXF_GDS(IJK) = VOL_V(IJK) *                             &
                  AVG_Y(F_GDS(IJK),F_GDS(IJKN),J)
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE VF_GS_Y

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: VF_gs_Z                                                 C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
!                                                                      C
!  Purpose: Calculate the average drag coefficient at i, j, k+1/2 and  C
!           multiply with W-momentum cell volume.                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE VF_GS_Z

      USE geometry
      USE indices
      USE compar
      USE drag
      USE discretelement
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: K, IJK, IJKT
!-----------------------------------------------

      IF(DISCRETE_ELEMENT .AND. .NOT.DES_ONEWAY_COUPLED) THEN
         DO IJK = ijkstart3, ijkend3
            IF (IP_AT_T(IJK)) THEN
               VXF_GDS(IJK) = ZERO
            ELSE
               K = K_OF(IJK)
               IJKT = TOP_OF(IJK)
               VXF_GDS(IJK) = VOL_W(IJK) *                             &
                  AVG_Z(F_GDS(IJK),F_GDS(IJKT),K)
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE VF_GS_Z

