!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_RO_g                                                C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!                                                                      C
!  Purpose: Initialize the gas densities                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_RO_G(ro_g,rop_g,p_g,ep_g)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar   , only: istart3, jstart3, kstart3, iend3, jend3, kend3
      USE eos      , only: EOSG
      USE fld_const, only: mw_avg, ro_g0
      USE param1   , only: UNDEFINED
      USE functions, only: wall_at

      IMPLICIT NONE

      double precision, intent(inout) ::  ro_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) ::   p_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) ::  ep_g(istart3:iend3,jstart3:jend3,kstart3:kend3)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer :: i,j,k
!-----------------------------------------------
      IF (RO_G0 == UNDEFINED) THEN   ! compressible case

         ! Calculate ro_g and rop_g in all fluid and flow boundary cells
         do k = kstart3, kend3
           do j = jstart3, jend3
             do i = istart3, iend3

! set_bc0 will have already defined ro_g and rop_g in MI and PI
! boundary cells (redundant-remove in set_bc0?)
             IF (.NOT.wall_at(i,j,k)) THEN
               RO_G(i,j,k) = EOSG(mw_avg,P_G(i,j,k),295.15d0)
               ROP_G(i,j,k) = EP_G(I,J,K)*RO_G(i,j,k)
            ENDIF

         end do
         end do
         end do

      ELSE   ! incompressible case

         do k = kstart3, kend3
           do j = jstart3, jend3
             do i = istart3, iend3

            IF (.NOT.wall_at(i,j,k)) THEN
! assign ro_g and calculate rop_g in all fluid and flow boundary cells
! set_constprop will have already defined ro_g in fluid and flow
! boundary cells (redundant- remove here?)
               RO_G(i,j,k) = RO_G0
! set_bc0 will have already defined rop_g in MI and PI boundary cells
! (redundant-remove in set_bc0?)
               ROP_G(i,j,k) = EP_G(I,J,K)*RO_G(i,j,k)
            ENDIF
         end do
         end do
         end do
      ENDIF

      RETURN
      END SUBROUTINE SET_RO_G
