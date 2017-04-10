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

      use residual, only: max_resid_index, nresid
      use residual, only: resid_p, resid_u, resid_v, resid_w
      use residual, only: resid_index, resid_string, resid_x
      use residual, only: sum5_resid, group_resid, resid_prefix, resid_grp, hydro_grp
      use run, only: detect_stall
      use residual, only: tol_resid, tol_diverge

      use param1, only: zero, undefined_i, is_undefined, large_number

      implicit none

      ! Iteration number
      integer(c_int), intent(in) :: nit

      real(c_real), intent(inout) :: resid(8,2)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! sum of residuals
      real(c_real) :: SUM_RESID
! max of residuals
      real(c_real) :: maxres
! index
      integer :: lc, maxl
!-----------------------------------------------

! Normalize residuals
      do lc=1, 8
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

! find the variable with maximum residual
      if (is_undefined(resid_index(max_resid_index,1))) then
         maxres = -1.0d0
         maxl = 0
         do lc = 1, 8
            if (resid(lc,1) >= maxres) then
               maxres = resid(lc,1)
               maxl = lc
            endif
         enddo
         write (resid_string(max_resid_index), '(a1,i1)') &
            resid_prefix(maxl), 0
      endif
      if (group_resid) resid_grp(hydro_grp) = sum_resid

! Every 5 iterations detect whether the run is stalled by checking
! that the total residual has decreased.
      if(detect_stall .and. mod(nit,5) == 0) then
         if(nit > 10) then
            if(sum5_resid <= sum_resid) then
! The run is stalled. Reduce the time step.
               check_convergence = 2
               return
            endif
         endif
         sum5_resid = sum_resid
      endif

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
