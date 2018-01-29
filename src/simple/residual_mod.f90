      module residual

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      ! Tolerance in residuals allowed for convergence
      real(c_real) :: tol_resid
      ! Minimum residual for declaring divergence
      real(c_real) :: tol_diverge
      ! Factor for normalizing the residual of gas cont. eq.
      real(c_real) :: norm_g

      integer, parameter :: nresid = 8

      integer, parameter :: resid_p  =  1 ! Pressure
      integer, parameter :: resid_u  =  2 ! U-velocity
      integer, parameter :: resid_v  =  3 ! V-velocity
      integer, parameter :: resid_w  =  4 ! W-velocity

! sum of residuals every 5 iterations
      real(c_real) :: SUM5_resid

   contains


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_resid_vel                                          !
!  Purpose: Calculate residuals for momentum equations                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine calc_resid_vel(dir, lo, hi, alo, ahi, &
         v0lo, v0hi, v1lo, v1hi, v2lo, v2hi, &
         vel, vels1, vels2, A_m, b_m, mask, num, den)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use param  , only: small_number
      use matrix  , only: e, w, s, n, t, b

      implicit none

      integer     , intent(in   ) :: dir
      integer     , intent(in   ) ::  lo(3), hi(3)
      integer     , intent(in   ) :: alo(3),ahi(3)
      integer     , intent(in   ) :: v0lo(3),v0hi(3)
      integer     , intent(in   ) :: v1lo(3),v1hi(3)
      integer     , intent(in   ) :: v2lo(3),v2hi(3)

      ! primary velocity component
      real(c_real), intent(IN) :: vel&
         (v0lo(1):v0hi(1),v0lo(2):v0hi(2),v0lo(3):v0hi(3))

      ! other components used here only for scaling
      real(c_real), intent(IN) :: vels1&
         (v1lo(1):v1hi(1),v1lo(2):v1hi(2),v1lo(3):v1hi(3))
      real(c_real), intent(IN) :: vels2&
         (v2lo(1):v2hi(1),v2lo(2):v2hi(2),v2lo(3):v2hi(3))

      ! Septadiagonal matrix A_m
      real(c_real) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)

      ! Vector b_m
      real(c_real) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      real(c_real) :: mask&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      ! Residual ID, numerator and denominator
      real(c_real), intent(  out) :: num, den

!-----------------------------------------------
!     Local variables
!-----------------------------------------------
      ! Velocity magnitude
      real(c_real) :: magvel

      ! Tangential velocity averaged onto normal face
      real(c_real) :: uavg, vavg, wavg

      ! Numerators and denominators
      real(c_real)   :: num1, den1

      ! Indices
      integer(c_int) :: i, j, k
      integer(c_int) :: llo(3),lhi(3) 

      llo = alo; lhi = ahi

!     Evaluate the residual at cell (i,j,k):
!     RESp = B-sum(Anb*VARnb)-Ap*VARp
!       (where nb = neighbor cells and p = center/0 cell)

      ! Note these are the dimensions of the normal velocity array which is face-based
      do k = lo(3),  hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               num1 = b_m(i,j,k) - (&
                  a_m(i,j,k,0)*vel(i,j,k)+&
                  a_m(i,j,k,e)*vel(i+1,j,k) + &
                  a_m(i,j,k,w)*vel(i-1,j,k) + &
                  a_m(i,j,k,n)*vel(i,j+1,k) + &
                  a_m(i,j,k,s)*vel(i,j-1,k) + &
                  a_m(i,j,k,t)*vel(i,j,k+1) + &
                  a_m(i,j,k,b)*vel(i,j,k-1) )

               ! Ignore momentum residual in stagnant regions.  Need an alternative
               ! criteria for residual scaling for such cases.
               if (dir.eq.1) then
                  vavg = 0.25d0 * (vels1(i  ,j,k) + vels1(i  ,j+1,k) + & 
                                   vels1(i-1,j,k) + vels1(i-1,j+1,k) )
                  wavg = 0.25d0 * (vels2(i  ,j,k) + vels2(i  ,j,k+1) + & 
                                   vels2(i-1,j,k) + vels2(i-1,j,k+1) )
                  magvel = sqrt(vel(i,j,k)**2 + vavg**2 + wavg**2)
               else if (dir.eq.2) then
                  wavg = 0.25d0 * (vels1(i,j  ,k) + vels1(i,j  ,k+1) + & 
                                   vels1(i,j-1,k) + vels1(i,j-1,k+1) )
                  uavg = 0.25d0 * (vels2(i,j  ,k) + vels2(i+1,j  ,k) + & 
                                   vels2(i,j-1,k) + vels2(i+1,j-1,k) )
                  magvel = sqrt(vel(i,j,k)**2 + wavg**2 + uavg**2)
               else if (dir.eq.3) then
                  uavg = 0.25d0 * (vels1(i,j,k  ) + vels1(i+1,j,k  ) + & 
                                   vels1(i,j,k-1) + vels1(i+1,j,k-1) )
                  vavg = 0.25d0 * (vels2(i,j,k  ) + vels2(i,j+1,k  ) + & 
                                   vels2(i,j,k-1) + vels2(i,j+1,k-1) )
                  magvel = sqrt(vel(i,j,k)**2 + uavg**2 + vavg**2)
               end if

               if (magvel > small_number) then
                  num1 = abs(num1)
                  den1 = abs(a_m(i,j,k,0)*magvel)

                  ! Adding to terms that are accumulated
                  num = num + num1/mask(i,j,k)
                  den = den + den1/mask(i,j,k)
               endif
            enddo
         enddo
      enddo

   end subroutine calc_resid_vel


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_RESID_pp                                           !
!  Purpose: Calculate residuals for pressure correction equation       !
!                                                                      !
!  Comments:                                                           !
!  for correction equations the convergence for the corrections must   !
!  go to zero, therefore the vector b must go to zero. this value      !
!  cannot be normalized as the other equations are since the           !
!  denominator here will vanish.  thus the residual is normalized      !
!  based on its value in the first iteration                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine calc_resid_pp(alo, ahi, lo, hi, b_m, b_mmax, num, den)

      implicit none

      integer, intent(in   ) :: alo(3),ahi(3),lo(3),hi(3)
      real(c_real), intent(inout) :: num, den

!   Vector b_m
      real(c_real) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      real(c_real) :: b_mmax&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))

      integer :: i, j, k

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               num = num  + abs(b_m(i,j,k))
            enddo
         enddo
      enddo

      if(norm_g <= epsilon(0.0)) then
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  den = den + 10.d0*abs(b_mmax(i,j,k))
               enddo
            enddo
         enddo
      else
         den = den + (hi(3)-lo(3)+1)*(hi(2)-lo(2)+1)*(hi(1)-lo(1)+1)
      endif

      return
   end subroutine calc_resid_pp
end module residual
