      module residual

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      Use param, only: DIM_n

      integer, parameter :: MAX_resid_INDEX = 8
!
      integer, parameter :: Nresid = 8 + DIM_N

      integer, parameter :: resid_p  =  1 ! Pressure
      integer, parameter :: resid_ro =  2 ! Density, volume fraction
      integer, parameter :: resid_u  =  3 ! U-velocity
      integer, parameter :: resid_v  =  4 ! V-velocity
      integer, parameter :: resid_w  =  5 ! W-velocity
      integer, parameter :: resid_t  =  6 ! Temperature
      integer, parameter :: resid_th =  7 ! Granular temperature
      integer, parameter :: resid_sc =  8 ! User-defined scalar
      integer, parameter :: resid_ke =  9 ! K-epsilon equations
      integer, parameter :: resid_x  = 10 ! Mass fraction

! Group residuals by equation
      integer, parameter :: HYDRO_GRP   = 1     !hydrodynamics
      integer, parameter :: THETA_GRP   = 2     !Granular Energy
      integer, parameter :: ENERGY_GRP  = 3     !Energy
      integer, parameter :: SCALAR_GRP  = 5     !Scalars
      integer, parameter :: KE_GRP      = 6     !K-Epsilon

! Prefix of Residuals string
      integer, parameter :: NPREFIX  = 10
      CHARACTER, parameter, DIMENSION(NPREFIX) :: resid_PREFIX = &
        (/ 'P', 'R', 'U', 'V', 'W', 'T', 'G', 'S', 'K', 'X' /)

      ! Average residual
      real(c_real) :: resid(Nresid)

      ! Residual Numerator
      real(c_real) :: num_resid(Nresid)=0.0d0

      ! Residual Denominator
      real(c_real) :: den_resid(Nresid)=0.0d0

! sum of residuals every 5 iterations
      real(c_real) :: SUM5_resid

! Residual sum within a group of equations
      LOGICAL          :: GROUP_resid
      real(c_real) :: resid_GRP(6)

! Residuals to be printed out
      CHARACTER(LEN=4) :: resid_STRING(MAX_resid_INDEX)
      CHARACTER(LEN=8) :: resid_GRP_STRING(6)

! Indices of residuals to be printed out
      integer :: resid_INDEX(MAX_resid_INDEX, 2)

! For checking the over-all fluid mass balance
      real(c_real) :: accum_resid_g

   contains

!---------------------------------------------------------------------

      subroutine set_resid_p(val) &
          bind(C, name="set_resid_p")

      real(c_real)  , intent(in) :: val

      resid(resid_p) = val

      end subroutine set_resid_p


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: calc_resid_vel                                          C
!  Purpose: Calculate residuals for momentum equations                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine calc_resid_vel(alo, ahi, &
         v0lo, v0hi, v1lo, v1hi, v2lo, v2hi, &
         vel, vels1, vels2, A_m, b_m, &
         eq_id, num, den)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use param1  , only: small_number, zero
      use matrix  , only: e, w, s, n, t, b
      use geometry, only: domlo

      implicit none

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

      ! Residual ID, numerator and denominator
      integer,      intent(in   ) :: eq_id
      real(c_real), intent(  out) :: num, den

!-----------------------------------------------
!     Local variables
!-----------------------------------------------
      ! Velocity magnitude
      real(c_real) :: magvel

      ! Indices
      INTEGER :: i, j, k

      ! Numerators and denominators
      real(c_real) :: num1, den1

      ! Number of fluid cells

!-----------------------------------------------
      integer :: llo(3),lhi(3)

      llo = alo
      lhi = ahi

      if(eq_id == resid_u .and. alo(1) /= domlo(1)-1) llo(1) = alo(1)+1
      if(eq_id == resid_v .and. alo(2) /= domlo(2)-1) llo(2) = alo(2)+1
      if(eq_id == resid_w .and. alo(3) /= domlo(3)-1) llo(3) = alo(3)+1

! initializing

      do k = llo(3),lhi(3)
         do j = llo(2),lhi(2)
            do i = llo(1),lhi(1)

! evaluating the residual at cell (i,j,k):
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)

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
               magvel = sqrt(vel(i,j,k)**2 + vels1(i,j,k)**2+ vels2(i,j,k)**2)

               if (magvel > small_number) then
                  num1 = abs(num1)
                  den1 = abs(a_m(i,j,k,0)*magvel)

                  ! Adding to terms that are accumulated
                  num = num + num1
                  den = den + den1
               endif
            enddo
         enddo
      enddo

   end subroutine calc_resid_vel

end module residual
