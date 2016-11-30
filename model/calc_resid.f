!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_C                                            C
!  Purpose: Calculate residuals for continuity equations               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_C(VAR, A_M, B_M, M, NUM, DEN, &
         RESID, MAX_RESID, IJK_RESID)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use param, ONLY: DIMENSION_3
      use param1, ONLY: ZERO, ONE, UNDEFINED
      use matrix, only: e, w, s, n, t, b
      use geometry , only: do_k, ijkmax2
      use compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3, myPE, numPEs
      use compar   , only: ijksize3_all
      use run      , only: debug_resid

      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable being evaluated
      DOUBLE PRECISION, INTENT(IN) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
!   Vector b_m
      DOUBLE PRECISION :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Phase index
      INTEGER, INTENT(IN) :: M
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IJKW, IJKS, IJKB, IJKE, IJKN, IJKT
      INTEGER :: I, J, K
      INTEGER :: i_resid, j_resid, k_resid
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc
!-----------------------------------------------

! initializing values
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
          RESID_IJK(IJK) = ZERO
      ENDDO
      ENDDO
      ENDDO

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
         IF (fluid_at(i,j,k)) THEN
            IJKW = FUNIJK(iwest(i,j,k),j,k)
            IJKE = FUNIJK(ieast(i,j,k),j,k)
            IJKS = FUNIJK(i,jsouth(i,j,k),k)
            IJKN = FUNIJK(i,jnorth(i,j,k),k)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(I,J,K) - (A_M(I,J,K,0)*VAR(IJK)+A_M(I,J,K,E)*VAR(IJKE)+&
               A_M(I,J,K,W)*VAR(IJKW)+A_M(I,J,K,N)*VAR(IJKN)+A_M(I,J,K,S)*VAR(&
               IJKS))

            IF (DO_K) THEN
               IJKB = FUNIJK(i,j,kbot(i,j,k))
               IJKT = FUNIJK(i,j,ktop(i,j,k))
               NUM1 = NUM1 - (A_M(I,J,K,T)*VAR(IJKT)+A_M(I,J,K,B)*VAR(IJKB))
            ENDIF

            NUM1 = ABS(NUM1)
            DEN1 = ABS(A_M(I,J,K,0)*VAR(IJK))
! storing value of residual at each ijk location
            RESID_IJK(IJK) = NUM1

! adding to terms that are accumulated
            NCELLS = NCELLS + 1
            NUM = NUM + NUM1
            DEN = DEN + DEN1
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      IF(.not.debug_resid) RETURN

! Collecting all the information among all the procesors -
! determining the global sum
      ! call global_all_sum(NUM)
      ! call global_all_sum(DEN)
      ! call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
         IF (RESID_IJK(IJK) > MAX_RESID) THEN
            IJK_RESID = IJK
            i_RESID = i
            j_RESID = j
            k_RESID = k
            MAX_RESID = RESID_IJK( IJK_RESID )
         ENDIF
      ENDDO
      ENDDO
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK(i_resid,j_resid,k_resid)
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      ! call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      max_resid_gl = max_resid_l
      ijk_resid_gl = ijk_resid_l
      ! call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      ! call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Calling to determine the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2
      do nproc=0,NumPEs-1
         if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
            IJK_RESID = IJK_RESID_GL(nproc)
         endif
      enddo

! Normalizing the residual
      IF (DEN > ZERO) THEN
         RESID = NUM/DEN
         MAX_RESID = NCELLS*MAX_RESID/DEN
      ELSEIF (NUM == ZERO) THEN
         RESID = ZERO
         MAX_RESID = ZERO
         IJK_RESID = 0
      ELSE
         RESID = UNDEFINED
         MAX_RESID = UNDEFINED
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!         CALL WRITE_ERROR ('CALC_RESID_C', LINE, 1)
      ENDIF

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_C

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_S                                            C
!  Purpose: Calculate residuals for scalar equations                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_S(VAR, A_M, B_M, M, NUM, DEN, &
                              RESID, MAX_RESID, IJK_RESID, TOL)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3, myPE, numPEs
      use compar  , only: ijksize3_all
      use param   , only: dimension_3
      use param1  , only: undefined, zero, one
      use matrix  , only: e, w, s, n, t, b
      use geometry, only: do_k, ijkmax2
      use run     , only: debug_resid

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable being evaluated
      DOUBLE PRECISION, INTENT(IN) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Phase index
      INTEGER, INTENT(IN) :: M
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
! Ignore residual calculation for scalar values below this
      DOUBLE PRECISION, INTENT(IN) :: TOL
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
      INTEGER :: i, j, k, i_resid, j_resid, k_resid
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc
!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
         RESID_IJK(IJK) = ZERO
      ENDDO
      ENDDO
      ENDDO

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)

         IF (fluid_at(i,j,k) .AND. ABS(VAR(IJK)) > TOL) THEN

            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IPJK = FUNIJK(iplus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IJPK = FUNIJK(i,jplus(i,j,k),k)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(I,J,K) - (A_M(I,J,K,0)*VAR(IJK)+A_M(I,J,K,E)*VAR(IPJK)+&
               A_M(I,J,K,W)*VAR(IMJK)+A_M(I,J,K,N)*VAR(IJPK)+A_M(I,J,K,S)*VAR(&
               IJMK))
            IF (DO_K) THEN
               IJKM = FUNIJK(i,j,kminus(i,j,k))
               IJKP = FUNIJK(i,j,kplus(i,j,k))
               NUM1 = NUM1 - (A_M(I,J,K,T)*VAR(IJKP)+A_M(I,J,K,B)*VAR(IJKM))
            ENDIF

            NUM1 = ABS(NUM1)
            DEN1 = ABS(A_M(I,J,K,0)*VAR(IJK))
! storing value of residual at each ijk location
            RESID_IJK(IJK) = NUM1

! adding to terms that are accumulated
            NCELLS = NCELLS + 1
            NUM = NUM + NUM1
            DEN = DEN + DEN1
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      IF(.not.debug_resid) RETURN

! Collecting all the information among all the procesors -
! determining the global sum
      ! call global_all_sum(NUM)
      ! call global_all_sum(DEN)
      ! call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
         IF (RESID_IJK(IJK) > MAX_RESID) THEN
               IJK_RESID = IJK
               i_resid = i
               j_resid = j
               k_resid = k
               MAX_RESID = RESID_IJK( IJK_RESID )
         ENDIF
      ENDDO
      ENDDO
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK(i_resid,j_resid,k_resid)
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      ! call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      max_resid_gl = max_resid_l
      ! call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      ijk_resid_gl = ijk_resid_l
      ! call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Determining the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2
      do nproc=0,NumPEs-1
        if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
           IJK_RESID = IJK_RESID_GL(nproc)
        endif
      enddo

! Normalizing the residual
      IF (DEN > ZERO) THEN
         RESID = NUM/DEN
         MAX_RESID = NCELLS*MAX_RESID/DEN
      ELSEIF (NUM == ZERO) THEN
         RESID = ZERO
         MAX_RESID = ZERO
         IJK_RESID = 0
      ELSE
         RESID = UNDEFINED
         MAX_RESID = UNDEFINED
!         WRITE(LINE,*)'Message: All center coefficients are zero.'
!         CALL WRITE_ERROR('CALC_RESID_S', LINE, 1)
      ENDIF

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_S

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

      SUBROUTINE CALC_RESID_PP(B_M, NORM, NUM, DEN, RESID, MAX_RESID, &
         IJK_RESID)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3, myPE, numPEs
      use param1  , only: large_number, zero, one
      use matrix  , only: e, w, s, n, t, b
      use geometry, only: ijkmax2
      use run     , only: debug_resid

      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
!   Vector b_m
      DOUBLE PRECISION :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Normalization factor
      DOUBLE PRECISION, INTENT(IN) :: NORM
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK
      INTEGER :: i, j, k, i_resid, j_resid, k_resid
! Number of fluid cells
      INTEGER :: NCELLS
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! New local variables for DMP version
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc
!-----------------------------------------------

! initializing values
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0
      DEN1 = ONE
      IJK_RESID = 1

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
         IF (fluid_at(i,j,k)) THEN

! evaluating the residual at cell ijk:
            NUM1 = ABS(B_M(I,J,K))
            IF (NUM1 > MAX_RESID) THEN
               MAX_RESID = NUM1
               IJK_RESID = IJK
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


      IF(.not.debug_resid) THEN
! Collecting all the information among all the procesors
         ! call global_all_sum(NUM)
         ! call global_all_sum(DEN)

! Normalizing the residual
         IF (DEN*NORM > ZERO) THEN
! if norm=1 then this simply becomes an unscaled 'average' residual
            RESID = NUM/(DEN*NORM)
         ELSEIF (NUM == ZERO) THEN
            RESID = ZERO
         ELSE
            RESID = LARGE_NUMBER
         ENDIF
      ELSE   ! if(debug_resid) branch

! Collecting all the information among all the procesors -
! determining the global sum
         ! call global_all_sum(NUM)
         ! call global_all_sum(DEN)
         ! call global_all_sum(NCELLS)

! Determining the max residual
         do nproc=0,NumPEs-1
            if(nproc.eq.myPE) then
               MAX_RESID_L(nproc) = MAX_RESID
               IJK_RESID_L(nproc) = FUNIJK(i_resid,j_resid,k_resid)
            else
               MAX_RESID_L(nproc) = 0.0
               IJK_RESID_L(nproc) = 0
            endif
         enddo

! Determining the maximum among all the procesors
         ! call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
         max_resid_gl = max_resid_l
         ! call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
         ijk_resid_gl = ijk_resid_l
         ! call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Determining the global IJK location w.r.t. serial version
         IJK_RESID = IJKMAX2
         do nproc=0,NumPEs-1
            if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
               IJK_RESID = IJK_RESID_GL(nproc)
            endif
         enddo

! Normalizing the residual
         IF (DEN*NORM > ZERO) THEN
            RESID = NUM/(DEN*NORM)
            MAX_RESID = NCELLS*MAX_RESID/(DEN*NORM)
         ELSEIF (NUM == ZERO) THEN
            RESID = ZERO
            MAX_RESID = ZERO
            IJK_RESID = 0
         ELSE
            RESID = LARGE_NUMBER
            MAX_RESID = LARGE_NUMBER
!            WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!             CALL WRITE_ERROR ('CALC_RESID_pp', LINE, 1)
         ENDIF

      ENDIF   ! end if/else debug_resid branch

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_PP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_U                                            C
!  Purpose: Calculate residuals for u-momentum equations               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_U(U_M, V_M, W_M, A_M, B_M, M, NUM, DEN, &
         RESID, MAX_RESID, IJK_RESID)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use param   , only: dimension_3
      use param1  , only: large_number, small_number, zero, one
      use matrix  , only: e, w, s, n, t, b
      use geometry, only: do_k, ijkmax2
      use compar  , only: ijksize3_all
      use compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3, myPE, numPEs
      use run     , only: debug_resid

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! U velocity (x-dir)
      DOUBLE PRECISION, INTENT(IN) :: U_m(DIMENSION_3)
! V velocity (y-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: V_m(DIMENSION_3)
! W velocity (z-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: W_m(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Phase index
      INTEGER, INTENT(IN) :: M
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
!-----------------------------------------------
!     Local variables
!-----------------------------------------------
! Velocity magnitude
      DOUBLE PRECISION :: VEL
! Indices
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
      INTEGER :: i, j, k, i_resid, j_resid, k_resid
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc

!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
        RESID_IJK(IJK) = ZERO

! Skip walls where some values are undefined.
        IF(wall_at(i,j,k)) cycle

         IF (.NOT.ip_at_e(i,j,k)) THEN

            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IPJK = FUNIJK(iplus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IJPK = FUNIJK(i,jplus(i,j,k),k)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(I,J,K) - (A_M(I,J,K,0)*U_M(IJK)+&
               A_M(I,J,K,E)*U_M(IPJK)+A_M(I,J,K,W)*U_M(IMJK)+&
               A_M(I,J,K,N)*U_M(IJPK)+A_M(I,J,K,S)*U_M(IJMK))
            IF (DO_K) THEN
               IJKM = funijk(i,j,kminus(i,j,k))
               IJKP = funijk(i,j,kplus(i,j,k))
               NUM1 = NUM1 - (A_M(I,J,K,T)*U_M(IJKP)+A_M(I,J,K,B)*U_M(IJKM))
            ENDIF

! Ignore momentum residual in stagnant regions.  Need an alternative
! criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2)
            IF (VEL > SMALL_NUMBER) THEN
               NUM1 = ABS(NUM1)
               DEN1 = ABS(A_M(I,J,K,0)*VEL)
! storing value of residual at each ijk location
               RESID_IJK(IJK) = NUM1
! adding to terms that are accumulated
               NCELLS = NCELLS + 1
               NUM = NUM + NUM1
               DEN = DEN + DEN1
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      IF(.not.debug_resid) RETURN

! Collecting all the information among all the procesors -
! determining the global sum
      ! call global_all_sum(NUM)
      ! call global_all_sum(DEN)
      ! call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
         IF (RESID_IJK( IJK ) > MAX_RESID) THEN
            IJK_RESID = IJK
            i_resid = i
            j_resid = j
            k_resid = k
            MAX_RESID = RESID_IJK( IJK_RESID )
         ENDIF
      ENDDO
      ENDDO
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK(i_resid,j_resid,k_resid)
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      ! call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      max_resid_gl = max_resid_l
      ! call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      ijk_resid_gl = ijk_resid_l
      ! call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Determining the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2
      do nproc=0,NumPEs-1
        IF(MAX_RESID_GL(nproc).eq.MAX_RESID.and.&
           IJK_RESID_GL(nproc).lt.IJK_RESID) THEN
          IJK_RESID = IJK_RESID_GL(nproc)
        ENDIF
      ENDDO

! Normalizing the residual
      IF (DEN > ZERO) THEN
         RESID = NUM/DEN
         MAX_RESID = NCELLS*MAX_RESID/DEN
      ELSEIF (NUM == ZERO) THEN
         RESID = ZERO
         MAX_RESID = ZERO
         IJK_RESID = 0
      ELSE
         RESID = LARGE_NUMBER
         MAX_RESID = LARGE_NUMBER
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!         CALL WRITE_ERROR ('CALC_RESID_U', LINE, 1)
      ENDIF

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_U

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_V                                            C
!  Purpose: Calculate residuals for v-momentum equations               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_V(U_M, V_M, W_M, A_M, B_M, M, NUM, DEN, &
         RESID, MAX_RESID, IJK_RESID)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use param   , only: dimension_3
      use param1  , only: large_number, small_number, zero, one
      use matrix  , only: e, w, s, n, t, b
      use geometry, only: do_k, ijkmax2
      use compar  , only: ijksize3_all
      use compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3, myPE, numPEs
      use run     , only: debug_resid

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! U velocity (x-dir)
      DOUBLE PRECISION, INTENT(IN) :: U_m(DIMENSION_3)
! V velocity (y-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: V_m(DIMENSION_3)
! W velocity (z-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: W_m(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Phase index
      INTEGER, INTENT(IN) :: M
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
!-----------------------------------------------
!     Local variables
!-----------------------------------------------
! Velocity magnitude
      DOUBLE PRECISION :: VEL
! Indices
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
      INTEGER :: i, j, k, i_resid, j_resid, k_resid
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc

!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
        RESID_IJK(IJK) = ZERO

! Skip walls where some values are undefined.
        IF(wall_at(i,j,k)) cycle

         IF (.NOT.ip_at_n(i,j,k)) THEN

            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IPJK = FUNIJK(iplus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IJPK = FUNIJK(i,jplus(i,j,k),k)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(I,J,K) - (A_M(I,J,K,0)*V_M(IJK)+&
               A_M(I,J,K,E)*V_M(IPJK)+A_M(I,J,K,W)*V_M(IMJK)+&
               A_M(I,J,K,N)*V_M(IJPK)+A_M(I,J,K,S)*V_M(IJMK))
            IF (DO_K) THEN
               IJKM = funijk(i,j,kminus(i,j,k))
               IJKP = funijk(i,j,kplus(i,j,k))
               NUM1 = NUM1 - (A_M(I,J,K,T)*V_M(IJKP)+A_M(I,J,K,B)*V_M(IJKM))
            ENDIF

! Ignore momentum residual in stagnant regions.  Need an alternative
! criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2)
            IF (VEL > SMALL_NUMBER) THEN
               NUM1 = ABS(NUM1)
               DEN1 = ABS(A_M(I,J,K,0)*VEL)
! storing value of residual at each ijk location
               RESID_IJK(IJK) = NUM1
! adding to terms that are accumulated
               NCELLS = NCELLS + 1
               NUM = NUM + NUM1
               DEN = DEN + DEN1
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      if(.not.debug_resid) return

! Collecting all the information among all the procesors -
! determining the global sum
      ! call global_all_sum(NUM)
      ! call global_all_sum(DEN)
      ! call global_all_sum(NCELLS)

      MAX_RESID = RESID_IJK( 1 )
      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
            IJK = FUNIJK(i,j,k)
              IF (RESID_IJK( IJK ) > MAX_RESID) THEN
                IJK_RESID = IJK
                i_resid = i
                j_resid = j
                k_resid = k
                MAX_RESID = RESID_IJK( IJK_RESID )
             ENDIF
          ENDDO
        ENDDO
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK(i_resid,j_resid,k_resid)
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determine the maximum among all the procesors
      ! call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      max_resid_gl = max_resid_l
      ! call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      ijk_resid_gl = ijk_resid_l
      ! call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Determining the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2
      do nproc=0,NumPEs-1
        if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
          IJK_RESID = IJK_RESID_GL(nproc)
        endif
      enddo

! Normalizing the residual
      IF (DEN > ZERO) THEN
         RESID = NUM/DEN
         MAX_RESID = NCELLS*MAX_RESID/DEN
      ELSEIF (NUM == ZERO) THEN
         RESID = ZERO
         MAX_RESID = ZERO
         IJK_RESID = 0
      ELSE
         RESID = LARGE_NUMBER
         MAX_RESID = LARGE_NUMBER
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!         CALL WRITE_ERROR ('CALC_RESID_V', LINE, 1)
      ENDIF

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_V

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RESID_W                                            C
!  Purpose: Calculate residuals for w-momentum equations               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_RESID_W(U_M, V_M, W_M, A_M, B_M, M, NUM, DEN, &
         RESID, MAX_RESID, IJK_RESID)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use param   , only: dimension_3
      use param1  , only: large_number, small_number, zero, one
      use matrix  , only: e, w, s, n, t, b
      use geometry, only: do_k, ijkmax2
      use compar  , only: ijksize3_all
      use compar  , only: istart3, iend3, jstart3, jend3, kstart3, kend3, myPE, numPEs
      use run     , only: debug_resid

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! U velocity (x-dir)
      DOUBLE PRECISION, INTENT(IN) :: U_m(DIMENSION_3)
! V velocity (y-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: V_m(DIMENSION_3)
! W velocity (z-dir), used here for scaling
      DOUBLE PRECISION, INTENT(IN) :: W_m(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
!   Vector b_m
      DOUBLE PRECISION :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Phase index
      INTEGER, INTENT(IN) :: M
! Numerator and denominator
      DOUBLE PRECISION, INTENT(OUT) :: NUM, DEN
! Average value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: RESID
! Maximum value of Residual
      DOUBLE PRECISION, INTENT(OUT) :: MAX_RESID
! IJK of Maximum value of Residual
      INTEGER, INTENT(OUT) :: IJK_RESID
!-----------------------------------------------
!     Local variables
!-----------------------------------------------
! Velocity magnitude
      DOUBLE PRECISION :: VEL
! Indices
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
      INTEGER :: i, j, k, i_resid, j_resid, k_resid
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc

!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
        RESID_IJK(IJK) = ZERO

! Skip walls where some values are undefined.
        IF(wall_at(i,j,k)) cycle

         IF (.NOT.ip_at_t(i,j,k)) THEN
            IMJK = FUNIJK(iminus(i,j,k),j,k)
            IPJK = FUNIJK(iplus(i,j,k),j,k)
            IJMK = FUNIJK(i,jminus(i,j,k),k)
            IJPK = FUNIJK(i,jplus(i,j,k),k)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(I,J,K) - (A_M(I,J,K,0)*W_M(IJK)+&
               A_M(I,J,K,E)*W_M(IPJK)+A_M(I,J,K,W)*W_M(IMJK)+&
               A_M(I,J,K,N)*W_M(IJPK)+A_M(I,J,K,S)*W_M(IJMK))
            IF (DO_K) THEN
               IJKM = funijk(i,j,kminus(i,j,k))
               IJKP = funijk(i,j,kplus(i,j,k))
               NUM1 = NUM1 - (A_M(I,J,K,T)*W_M(IJKP)+A_M(I,J,K,B)*W_M(IJKM))
            ENDIF

! Ignore momentum residual in stagnant regions.  Need an alternative
! criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2)
            IF (VEL > SMALL_NUMBER) THEN
               NUM1 = ABS(NUM1)
               DEN1 = ABS(A_M(I,J,K,0)*VEL)
! storing value of residual at each ijk location
               RESID_IJK(IJK) = NUM1
! adding to terms that are accumulated
               NCELLS = NCELLS + 1
               NUM = NUM + NUM1
               DEN = DEN + DEN1
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      if(.not.debug_resid) return

! Collecting all the information among all the procesors -
! determining the global sum
      ! call global_all_sum(NUM)
      ! call global_all_sum(DEN)
      ! call global_all_sum(NCELLS)

      MAX_RESID = RESID_IJK( 1 )
      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
            IJK = FUNIJK(i,j,k)
            IF (RESID_IJK( IJK ) > MAX_RESID) THEN
              IJK_RESID = IJK
              MAX_RESID = RESID_IJK( IJK_RESID )
            ENDIF
          ENDDO
        ENDDO
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK(i_resid,j_resid,k_resid)
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      ! call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      max_resid_gl = max_resid_l
      ! call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      ijk_resid_gl = ijk_resid_l
      ! call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

! Determining the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2

      do nproc=0,NumPEs-1
        if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
          IJK_RESID = IJK_RESID_GL(nproc)
        endif
      enddo

! Normalizing the residual
      IF (DEN > ZERO) THEN
         RESID = NUM/DEN
         MAX_RESID = NCELLS*MAX_RESID/DEN
      ELSE IF (NUM == ZERO) THEN
         RESID = ZERO
         MAX_RESID = ZERO
         IJK_RESID = 0
      ELSE
         RESID = LARGE_NUMBER
         MAX_RESID = LARGE_NUMBER
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.'
!         CALL WRITE_ERROR ('CALC_RESID_W', LINE, 1)
      ENDIF

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_RESID_W
