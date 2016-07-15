!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_RO_g                                                C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!                                                                      C
!  Purpose: Initialize the gas densities                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_RO_G

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE eos, only: EOSG
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: IJK
!-----------------------------------------------
      IF (RO_G0 == UNDEFINED) THEN   ! compressible case

!!$omp parallel do private(IJK)
         DO IJK = ijkstart3, ijkend3
! calculate ro_g and rop_g in all fluid and flow boundary cells
            IF (.NOT.WALL_AT(IJK)) THEN
! set_bc0 will have already defined ro_g and rop_g in MI and PI
! boundary cells (redundant-remove in set_bc0?)
               RO_G(IJK) = EOSG(MW_AVG,P_G(IJK),295.15d0)
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK)
            ENDIF
         ENDDO

      ELSE   ! incompressible case

!!$omp   parallel do private(ijk)
         DO IJK = ijkstart3, ijkend3
            IF (.NOT.WALL_AT(IJK)) THEN
! assign ro_g and calculate rop_g in all fluid and flow boundary cells
! set_constprop will have already defined ro_g in fluid and flow
! boundary cells (redundant- remove here?)
               RO_G(IJK) = RO_G0
! set_bc0 will have already defined rop_g in MI and PI boundary cells
! (redundant-remove in set_bc0?)
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK)
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE SET_RO_G


