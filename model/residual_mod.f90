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
         i_resid, j_resid, k_resid, flag)

      use param1  , only: large_number, zero, one
      use run     , only: debug_resid

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
      INTEGER, INTENT(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

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
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0
      DEN1 = ONE

      DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
          DO I = slo(1),shi(1)
             IF(flag(i,j,k,1)==1) THEN

               ! evaluating the residual at cell (i,j,k):
               NUM1 = ABS(B_M(I,J,K))
               IF (NUM1 > MAX_RESID) THEN
                  MAX_RESID = NUM1
                  i_resid = i
                  j_resid = j
                  k_resid = k
            ENDIF

! adding to terms that are accumulated
            NCELLS = NCELLS + 1
            NUM = NUM + NUM1
            DEN = DEN + DEN1
         ENDIF

      ENDDO
      ENDDO
      ENDDO


      IF (.not.debug_resid) THEN
         ! Normalizing the residual
         IF (DEN*NORM > ZERO) THEN
            ! if norm=1 then this simply becomes an unscaled 'average' residual
            resid = NUM/(DEN*NORM)
         ELSEIF (abs(NUM) < epsilon(NUM)) THEN
            resid = zero
         ELSE
            resid = LARGE_NUMBER
         ENDIF
      ELSE   ! if(debug_resid) branch
         ! Normalizing the residual
         IF (abs(DEN*NORM) > epsilon(DEN*NORM)) THEN
            resid = NUM/(DEN*NORM)
            max_resid = NCELLS*MAX_RESID/(DEN*NORM)
         ELSEIF (abs(NUM) < epsilon(NUM)) THEN
            resid = zero
            max_resid = zero
            i_resid = 0
            j_resid = 0
            k_resid = 0
         ELSE
            resid = LARGE_NUMBER
            max_resid = LARGE_NUMBER
         ENDIF
      ENDIF   ! end if/else debug_resid branch

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_PP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_U                                            C
!  Purpose: Calculate residuals for u-momentum equation.               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_RESID_VEL(slo, shi, lo, hi, &
         vel,   vnlo,  vnhi, &
         vels1, vt1lo, vt1hi, &
         vels2, vt2lo, vt2hi, A_M, B_M, NUM, DEN, &
         RESID, MAX_RESID, i_resid, j_resid, k_resid, flag)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use param1  , only: large_number, small_number, zero, one
      use matrix  , only: e, w, s, n, t, b
      use run     , only: debug_resid

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: slo(3),shi(3),lo(3),hi(3)
      integer(c_int), intent(in   ) :: vnlo(3),  vnhi(3)
      integer(c_int), intent(in   ) :: vt1lo(3), vt1hi(3)
      integer(c_int), intent(in   ) :: vt2lo(3), vt2hi(3)

      ! primary velocity component
      real(c_real), INTENT(IN) :: vel&
         (vnlo(1):vnhi(1),vnlo(2):vnhi(2),vnlo(3):vnhi(3))

      ! other components used here only for scaling
      real(c_real), INTENT(IN) :: vels1&
         (vt1lo(1):vt1hi(1),vt1lo(2):vt1hi(2),vt1lo(3):vt1hi(3))
      real(c_real), INTENT(IN) :: vels2&
         (vt2lo(1):vt2hi(1),vt2lo(2):vt2hi(2),vt2lo(3):vt2hi(3))

      ! Septadiagonal matrix A_m
      real(c_real) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      INTEGER, INTENT(in) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      ! Numerator and denominator
      real(c_real), INTENT(OUT) :: NUM, DEN

      ! Average value of Residual
      real(c_real), INTENT(OUT) :: RESID

      ! Maximum value of Residual
      real(c_real), INTENT(OUT) :: MAX_RESID

      ! (i,j,k) of Maximum value of Residual
      INTEGER, INTENT(OUT) :: i_resid, j_resid, k_resid

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

!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

      DO K = vnlo(3)+1,vnhi(3)-1
         DO J = vnlo(2)+1,vnhi(2)-1
            DO I = vnlo(1)+1,vnhi(1)-1
               RESID_IJK(i,j,k) = ZERO


! Skip walls where some values are undefined.
               IF(flag(i,j,k,1)>=100) cycle

! This isn't correct. The flag "2" should be axis dependent.
               IF (flag(i,j,k,2) > 1000) THEN

! evaluating the residual at cell (i,j,k):
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)

                  NUM1 = B_M(I,J,K) - (&
                     A_M(I,J,K,0)*vel(i,j,k)+&
                     A_M(I,J,K,E)*vel( iplus(i,j,k),j,k)+&
                     A_M(I,J,K,W)*vel(iminus(i,j,k),j,k)+&
                     A_M(I,J,K,N)*vel(i, jplus(i,j,k),k)+&
                     A_M(I,J,K,S)*vel(i,jminus(i,j,k),k)+&
                     A_M(I,J,K,T)*vel(i,j, kplus(i,j,k))+&
                     A_M(I,J,K,B)*vel(i,j,kminus(i,j,k)))

! Ignore momentum residual in stagnant regions.  Need an alternative
! criteria for residual scaling for such cases.
                  magvel = SQRT(vel(i,j,k)**2 + vels1(i,j,k)**2+ vels2(i,j,k)**2)

                  IF (magvel > SMALL_NUMBER) THEN
                     NUM1 = ABS(NUM1)
                     DEN1 = ABS(A_M(I,J,K,0)*magvel)

! Storing value of residual at each (i,j,k) location
                     RESID_IJK(i,j,k) = NUM1

! Adding to terms that are accumulated
                     NCELLS = NCELLS + 1
                     NUM = NUM + NUM1
                     DEN = DEN + DEN1
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF (.not.debug_resid) RETURN

      i_RESID = 1
      j_RESID = 1
      k_RESID = 1

      max_resid = resid_ijk( i_resid, j_resid, k_resid)

      DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
          DO I = slo(1),shi(1)
            IF (resid_ijk( i,j,k ) > max_resid) THEN
              i_resid = i
              j_resid = j
              k_resid = k
              max_resid = resid_ijk(i,j,k)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ! Normalizing the residual
      IF (DEN > ZERO) THEN
         resid = NUM/DEN
         max_resid = NCELLS * max_resid / DEN
      ELSEIF (abs(NUM) < epsilon(NUM)) THEN
         resid = ZERO
         max_resid = ZERO
         i_resid = 0
         j_resid = 0
         k_resid = 0
      ELSE
         resid = LARGE_NUMBER
         max_resid = LARGE_NUMBER
      ENDIF

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

   END SUBROUTINE CALC_RESID_VEL

END MODULE residual
