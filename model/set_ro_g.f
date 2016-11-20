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
      USE compar, only:  istart3, jstart3, kstart3, iend3, jend3, kend3
      USE constant
      USE eos, only: EOSG
      USE fldvar
      USE param1   , only: UNDEFINED
      USE functions, only: funijk, wall_cell
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer :: ijk,i,j,k
!-----------------------------------------------
      IF (RO_G0 == UNDEFINED) THEN   ! compressible case

         ! Calculate ro_g and rop_g in all fluid and flow boundary cells
         do k = kstart3, kend3
           do j = jstart3, jend3
             do i = istart3, iend3

             ijk = funijk(i,j,k)

! set_bc0 will have already defined ro_g and rop_g in MI and PI
! boundary cells (redundant-remove in set_bc0?)
             IF (.NOT.wall_cell(i,j,k)) THEN
               RO_G(IJK) = EOSG(MW_AVG,P_G(IJK),295.15d0)
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK)
            ENDIF

         end do
         end do
         end do

      ELSE   ! incompressible case

         do k = kstart3, kend3
           do j = jstart3, jend3
             do i = istart3, iend3

             ijk = funijk(i,j,k)

            IF (.NOT.wall_cell(i,j,k)) THEN
! assign ro_g and calculate rop_g in all fluid and flow boundary cells
! set_constprop will have already defined ro_g in fluid and flow
! boundary cells (redundant- remove here?)
               RO_G(IJK) = RO_G0
! set_bc0 will have already defined rop_g in MI and PI boundary cells
! (redundant-remove in set_bc0?)
               ROP_G(IJK) = EP_G(IJK)*RO_G(IJK)
            ENDIF
         end do
         end do
         end do
      ENDIF

      RETURN
      END SUBROUTINE SET_RO_G


