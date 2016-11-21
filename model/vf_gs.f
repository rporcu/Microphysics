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
      INTEGER :: i,j,k, IJK, IJKE

! Calculate the combined effect for all discrete solids.
      IF(DISCRETE_ELEMENT .AND. .NOT.DES_ONEWAY_COUPLED) THEN
         do k = kstart3, kend3
           do j = jstart3, jend3
             do i = istart3, iend3

             ijk = funijk(i,j,k)
            IF(ip_at_e(i,j,k)) THEN
               VXF_GDS(IJK) = ZERO
            ELSE
               IJKE = FUNIJK(ieast(i,j,k),j,k)
               VXF_GDS(IJK) = VOL_U(IJK) *                             &
                  AVG_X(F_GDS(IJK),F_GDS(IJKE),I)
            ENDIF
             end do
           end do
         end do
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
      INTEGER :: i,j,k, IJK, IJKN
!-----------------------------------------------

      IF(DISCRETE_ELEMENT .AND. .NOT.DES_ONEWAY_COUPLED) THEN
         do k = kstart3, kend3
           do j = jstart3, jend3
             do i = istart3, iend3

             ijk = funijk(i,j,k)
            IF(ip_at_n(i,j,k)) THEN
               VXF_GDS(IJK) = ZERO
            ELSE
               IJKN = FUNIJK(i,jnorth(i,j,k),k)
               VXF_GDS(IJK) = VOL_V(IJK) *                             &
                  AVG_Y(F_GDS(IJK),F_GDS(IJKN),J)
            ENDIF
             end do
           end do
         end do
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
      INTEGER :: i, j, k, IJK, IJKT
!-----------------------------------------------

      IF(DISCRETE_ELEMENT .AND. .NOT.DES_ONEWAY_COUPLED) THEN
         do k = kstart3, kend3
           do j = jstart3, jend3
             do i = istart3, iend3

             ijk = funijk(i,j,k)

            IF (ip_at_t(i,j,k)) THEN
               VXF_GDS(IJK) = ZERO
            ELSE
               IJKT = FUNIJK(i,j,ktop(i,j,k))
               VXF_GDS(IJK) = VOL_W(IJK) *                             &
                  AVG_Z(F_GDS(IJK),F_GDS(IJKT),K)
            ENDIF
             end do
           end do
         end do
      ENDIF

      RETURN
      END SUBROUTINE VF_GS_Z

