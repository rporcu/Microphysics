MODULE CHECK_CONVERGENCE_MODULE

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_CONVERGENCE                                       C
!  Author: M. Syamlal                                 Date: 8-JUL-96   C
!                                                                      C
!  Purpose: Monitor convergence                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   integer(c_int) function check_convergence(nit, resid) &
      bind(C, name="check_convergence")

      use residual, only: nresid, sum5_resid
      use residual, only: resid_p, resid_u, resid_v, resid_w
      use residual, only: tol_resid, tol_diverge
      use param, only: zero, is_undefined, large_number

      implicit none

      ! Iteration number
      integer(c_int), intent(in) :: nit

      real(c_real), intent(inout) :: resid(8,2)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! sum of residuals
      real(c_real) :: SUM_RESID
! index
      integer :: lc
!-----------------------------------------------

! Normalize residuals
      do lc=1, nresid
         if (resid(lc,2) > zero) then
            resid(lc,1) = resid(lc,1)/resid(lc,2)
         elseif (abs(resid(lc,1)) < epsilon(resid(lc,1))) then
            resid(lc,1) = zero
         else
            resid(lc,1) = large_number
         endif
      enddo

! add pressure correction residual to momentum residuals
      sum_resid = resid(resid_p,1) + resid(resid_u,1) + &
         resid(resid_v,1) + resid(resid_w,1)

! Require at least two iterations.
      if(nit == 1) then
         check_convergence = 0
         return
      endif

! total residual
      if(sum_resid<=tol_resid) then
         check_convergence = 1          ! converged
      elseif (sum_resid>=tol_diverge ) then
         if (nit /= 1) then
            check_convergence = 2       ! diverged
         else
            check_convergence = 0       ! not converged
         endif
      else
         check_convergence = 0          ! not converged
      endif

   end function check_convergence
end module check_convergence_module
