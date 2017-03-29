module conv_pp_g_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CONV_Pp_g                                               !
!  Purpose: Determine convection terms for Pressure correction         !
!           equation.                                                  !
!                                                                      !
!  Notes: The off-diagonal coefficients calculated here must be        !
!         positive. The center coefficient and the source vector are   !
!         negative.                                                    !
!                                                                      !
!         Multiplication with factors d_e, d_n, and d_t are carried    !
!         out in source_pp_g.  Constant pressure boundaries are        !
!         handled by holding the Pp_g at the boundaries zero.  For     !
!         specified mass flow boundaries (part of) a's are calculated  !
!         here since b is calculated from a's in source_pp_g.  After   !
!         calculating b, a's are multiplied by d and at the flow       !
!         boundaries are set to zero.                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine conv_pp_g(ulo, uhi, vlo, vhi, wlo, whi, alo, ahi, &
         A_m, rop_ge, rop_gn, rop_gt, dx, dy, dz)

! Modules
!-----------------------------------------------
      use matrix   , only: e, w, s, n, t, b

      implicit none

      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)
      integer(c_int), intent(in   ) :: alo(3),ahi(3)
      real(c_real), intent(in   ) :: dx,dy,dz

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      real(c_real), intent(in   ) :: rop_ge&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: rop_gn&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: rop_gt&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

! Local variables
!-----------------------------------------------
! Indices
      integer ::  i,j,k
      real(c_real) :: axy, axz, ayz
!-----------------------------------------------

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

      ! Calculate convection fluxes through each of the faces
      do k = alo(3),ahi(3)
         do j = alo(2),ahi(2)
            do i = alo(1),ahi(1)

               A_m(i,j,k,e) = rop_ge(i  ,j,k)*ayz
               A_m(i,j,k,w) = rop_ge(i-1,j,k)*ayz

               A_m(i,j,k,n) = rop_gn(i,j  ,k)*axz
               A_m(i,j,k,s) = rop_gn(i,j-1,k)*axz

               A_m(i,j,k,t) = rop_gt(i,j,k  )*axy
               A_m(i,j,k,b) = rop_gt(i,j,k-1)*axy

            enddo
         enddo
      enddo

   end subroutine conv_pp_g

end module conv_pp_g_module
