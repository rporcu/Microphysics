      MODULE residual

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      Use param, only: DIM_n

      integer, parameter :: MAX_RESID_INDEX = 8
!
      integer, parameter :: NRESID = 8 + DIM_N

      integer, parameter :: RESID_p  =  1 ! Pressure
      integer, parameter :: RESID_ro =  2 ! Density, volume fraction
      integer, parameter :: RESID_u  =  3 ! U-velocity
      integer, parameter :: RESID_v  =  4 ! V-velocity
      integer, parameter :: RESID_w  =  5 ! W-velocity
      integer, parameter :: RESID_t  =  6 ! Temperature
      integer, parameter :: RESID_th =  7 ! Granular temperature
      integer, parameter :: RESID_sc =  8 ! User-defined scalar
      integer, parameter :: RESID_ke =  9 ! K-epsilon equations
      integer, parameter :: RESID_x  = 10 ! Mass fraction

! Group residuals by equation
      integer, parameter :: HYDRO_GRP   = 1     !hydrodynamics
      integer, parameter :: THETA_GRP   = 2     !Granular Energy
      integer, parameter :: ENERGY_GRP  = 3     !Energy
      integer, parameter :: SCALAR_GRP  = 5     !Scalars
      integer, parameter :: KE_GRP      = 6     !K-Epsilon

! Prefix of Residuals string
      integer, parameter :: NPREFIX  = 10
      CHARACTER, parameter, DIMENSION(NPREFIX) :: RESID_PREFIX = &
        (/ 'P', 'R', 'U', 'V', 'W', 'T', 'G', 'S', 'K', 'X' /)

! Average residual
      real(c_real) :: RESID(NRESID)
! Maximum residual
      real(c_real) :: MAX_RESID(NRESID)
! Residual Numerator
      real(c_real) :: NUM_RESID(NRESID)
! Residual Denominator
      real(c_real) :: DEN_RESID(NRESID)

! (i,j,k) location of maximum residual
      integer :: i_resid(nresid)
      integer :: j_resid(nresid)
      integer :: k_resid(nresid)

! sum of residuals every 5 iterations
      real(c_real) :: SUM5_RESID

! Residual sum within a group of equations
      LOGICAL          :: GROUP_RESID
      real(c_real) :: RESID_GRP(6)

! Residuals to be printed out
      CHARACTER(LEN=4) :: RESID_STRING(MAX_RESID_INDEX)
      CHARACTER(LEN=8) :: RESID_GRP_STRING(6)

! Indices of residuals to be printed out
      integer :: RESID_INDEX(MAX_RESID_INDEX, 2)

! For checking the over-all fluid mass balance
      real(c_real) :: accum_resid_g

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_pp                                           C
!  Purpose: Calculate residuals for pressure correction equation       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments:                                                           C
!  for correction equations the convergence for the corrections must   C
!  go to zero, therefore the vector b must go to zero. this value      C
!  cannot be normalized as the other equations are since the           C
!  denominator here will vanish.  thus the residual is normalized      C
!  based on its value in the first iteration                           C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_RESID_PP(slo,shi,lo,hi,&
         B_M, NORM, NUM, DEN, RESID, MAX_RESID, &
         i_resid, j_resid, k_resid)

      use param1  , only: large_number, zero, one

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      !   Vector b_m
      real(c_real) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Normalization factor
      real(c_real), INTENT(IN) :: NORM

      ! Numerator and denominator
      real(c_real), INTENT(OUT) :: NUM, DEN

      ! Average value of Residual
      real(c_real), INTENT(OUT) :: RESID

      ! Maximum value of Residual
      real(c_real), INTENT(OUT) :: MAX_RESID

      ! (i,j,k) of Maximum value of Residual
      INTEGER, INTENT(OUT) :: i_resid, j_resid, k_resid
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: i, j, k
! Number of fluid cells
      INTEGER :: NCELLS
! Numerators and denominators
      real(c_real) :: NUM1, DEN1
!-----------------------------------------------

! initializing values
      num = zero
      den = zero
      max_resid = -one
      ncells = 0
      den1 = one

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

! evaluating the residual at cell (i,j,k):
               num1 = abs(b_m(i,j,k))
               if (num1 > max_resid) then
                  max_resid = num1
                  i_resid = i
                  j_resid = j
                  k_resid = k
               endif

! adding to terms that are accumulated
               ncells = ncells + 1
               num = num + num1
               den = den + den1

            enddo
         enddo
      enddo

! Normalizing the residual
      if (abs(den*norm) > epsilon(den*norm)) then
         resid = num/(den*norm)
         max_resid = ncells*max_resid/(den*norm)
      elseif (abs(num) < epsilon(num)) then
         resid = zero
         max_resid = zero
         i_resid = 0
         j_resid = 0
         k_resid = 0
      else
         resid = large_number
         max_resid = large_number
      endif

      return

      end subroutine calc_resid_pp


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_VEL                                          C
!  Purpose: Calculate residuals for momentum equations                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_RESID_VEL(slo, shi, lo, hi, &
         vel, vels1, vels2, A_M, B_M, NUM, DEN, &
         RESID, MAX_RESID, i_resid, j_resid, k_resid, axis)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use param1  , only: large_number, small_number, zero, one
      use matrix  , only: e, w, s, n, t, b

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      ! primary velocity component
      real(c_real), INTENT(IN) :: vel&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! other components used here only for scaling
      real(c_real), INTENT(IN) :: vels1&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN) :: vels2&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Septadiagonal matrix A_m
      real(c_real) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      ! Numerator and denominator
      real(c_real), INTENT(OUT) :: NUM, DEN

      ! Average value of Residual
      real(c_real), INTENT(OUT) :: RESID

      ! Maximum value of Residual
      real(c_real), INTENT(OUT) :: MAX_RESID

      ! (i,j,k) of Maximum value of Residual
      INTEGER, INTENT(OUT) :: i_resid, j_resid, k_resid

      character, intent(in) :: axis

!-----------------------------------------------
!     Local variables
!-----------------------------------------------
      ! Velocity magnitude
      real(c_real) :: magvel

      ! Indices
      INTEGER :: i, j, k

      ! Numerators and denominators
      real(c_real) :: NUM1, DEN1

      ! Number of fluid cells
      INTEGER :: NCELLS

      ! New local variables for DMP version
      real(c_real) :: resid_ijk&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))


      integer :: llo(3),lhi(3)
!-----------------------------------------------

      llo = lo
      lhi = hi

      if(axis == 'U') then
         lhi(1) = hi(1) -1
      elseif(axis == 'V') then
         lhi(2) = hi(2) -1
      elseif(axis == 'W') then
         lhi(3) = hi(3) -1
      else
         stop
      endif

! initializing
      num = zero
      den = zero
      max_resid = -one
      ncells = 0

      do k = llo(3),lhi(3)
         do j = llo(2),lhi(2)
            do i = llo(1),lhi(1)
               resid_ijk(i,j,k) = zero

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

! Storing value of residual at each (i,j,k) location
                  resid_ijk(i,j,k) = num1

! Adding to terms that are accumulated
                  ncells = ncells + 1
                  num = num + num1
                  den = den + den1
               endif
            enddo
         enddo
      enddo


      i_resid = 1
      j_resid = 1
      k_resid = 1

      max_resid = resid_ijk( i_resid, j_resid, k_resid)

      do k = llo(3),lhi(3)
         do j = llo(2),lhi(2)
            do i = llo(1),lhi(1)
               if (resid_ijk( i,j,k ) > max_resid) then
                  i_resid = i
                  j_resid = j
                  k_resid = k
                  max_resid = resid_ijk(i,j,k)
               endif
            enddo
         enddo
      enddo

      ! Normalizing the residual
      if (den > zero) then
         resid = num/den
         max_resid = ncells * max_resid / den
      elseif (abs(num) < epsilon(num)) then
         resid = zero
         max_resid = zero
         i_resid = 0
         j_resid = 0
         k_resid = 0
      else
         resid = large_number
         max_resid = large_number
      endif

      return

   end subroutine calc_resid_vel
end module residual
