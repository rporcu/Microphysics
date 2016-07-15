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
      USE param, ONLY: DIMENSION_3, DIMENSION_M
      USE param1, ONLY: ZERO, ONE, UNDEFINED
      use ambm, only: e, w, s, n, t, b
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable being evaluated
      DOUBLE PRECISION, INTENT(IN) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
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

!!$omp parallel do private( IJK )
      DO IJK = ijkstart3, ijkend3
          RESID_IJK(IJK) = ZERO
      ENDDO

!!$omp  parallel do private( IJK, IJKW, IJKS, IJKB, IJKE, IJKN, IJKT,  &
!!$omp&  NUM1, DEN1) &
!!$omp&  REDUCTION(+:NUM,DEN,NCELLS)
      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (FLUID_AT(IJK)) THEN
            IJKW = WEST_OF(IJK)
            IJKS = SOUTH_OF(IJK)
            IJKE = EAST_OF(IJK)
            IJKN = NORTH_OF(IJK)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*VAR(IJK)+A_M(IJK,E,M)*VAR(IJKE)+&
               A_M(IJK,W,M)*VAR(IJKW)+A_M(IJK,N,M)*VAR(IJKN)+A_M(IJK,S,M)*VAR(&
               IJKS))

            IF (DO_K) THEN
               IJKB = BOTTOM_OF(IJK)
               IJKT = TOP_OF(IJK)
               NUM1 = NUM1 - (A_M(IJK,T,M)*VAR(IJKT)+A_M(IJK,B,M)*VAR(IJKB))
            ENDIF

            NUM1 = ABS(NUM1)
            DEN1 = ABS(A_M(IJK,0,M)*VAR(IJK))
! storing value of residual at each ijk location
            RESID_IJK(IJK) = NUM1

! adding to terms that are accumulated
            NCELLS = NCELLS + 1
            NUM = NUM + NUM1
            DEN = DEN + DEN1
         ENDIF
      ENDDO

      IF(.not.debug_resid) RETURN

! Collecting all the information among all the procesors -
! determining the global sum
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (RESID_IJK(IJK) > MAX_RESID) THEN
            IJK_RESID = IJK
            MAX_RESID = RESID_IJK( IJK_RESID )
         ENDIF
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

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
      USE param
      USE param1
      use ambm, only: e, w, s, n, t, b
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run

      USE fldvar
      USE physprop
      USE toleranc

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable being evaluated
      DOUBLE PRECISION, INTENT(IN) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
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

!!$omp parallel do private( IJK )
      DO IJK = ijkstart3, ijkend3
         RESID_IJK(IJK) = ZERO
      ENDDO

!!$omp    parallel do &
!!$omp&   private(   IJK, IMJK, IJMK, IPJK, IJPK, IJKM, IJKP, &
!!$omp&   NUM1, DEN1) &
!!$omp&   REDUCTION(+:NUM, DEN,NCELLS)

      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

         IF (FLUID_AT(IJK) .AND. ABS(VAR(IJK)) > TOL) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*VAR(IJK)+A_M(IJK,E,M)*VAR(IPJK)+&
               A_M(IJK,W,M)*VAR(IMJK)+A_M(IJK,N,M)*VAR(IJPK)+A_M(IJK,S,M)*VAR(&
               IJMK))
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IJKP = KP_OF(IJK)
               NUM1 = NUM1 - (A_M(IJK,T,M)*VAR(IJKP)+A_M(IJK,B,M)*VAR(IJKM))
            ENDIF

            NUM1 = ABS(NUM1)
            DEN1 = ABS(A_M(IJK,0,M)*VAR(IJK))
! storing value of residual at each ijk location
            RESID_IJK(IJK) = NUM1

! adding to terms that are accumulated
            NCELLS = NCELLS + 1
            NUM = NUM + NUM1
            DEN = DEN + DEN1
         ENDIF
      ENDDO

      IF(.not.debug_resid) RETURN

! Collecting all the information among all the procesors -
! determining the global sum
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (RESID_IJK(IJK) > MAX_RESID) THEN
               IJK_RESID = IJK
               MAX_RESID = RESID_IJK( IJK_RESID )
         ENDIF
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

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
      USE param
      USE param1
      use ambm, only: e, w, s, n, t, b
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
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

      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (FLUID_AT(IJK)) THEN

! evaluating the residual at cell ijk:
            NUM1 = ABS(B_M(IJK,0))
            IF (NUM1 > MAX_RESID) THEN
               MAX_RESID = NUM1
               IJK_RESID = IJK
            ENDIF

! adding to terms that are accumulated
            NCELLS = NCELLS + 1
            NUM = NUM + NUM1
            DEN = DEN + DEN1
         ENDIF
      ENDDO


      IF(.not.debug_resid) THEN
! Collecting all the information among all the procesors
         call global_all_sum(NUM)
         call global_all_sum(DEN)

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
         call global_all_sum(NUM)
         call global_all_sum(DEN)
         call global_all_sum(NCELLS)

! Determining the max residual
         do nproc=0,NumPEs-1
            if(nproc.eq.myPE) then
               MAX_RESID_L(nproc) = MAX_RESID
               IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
            else
               MAX_RESID_L(nproc) = 0.0
               IJK_RESID_L(nproc) = 0
            endif
         enddo

! Determining the maximum among all the procesors
         call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
         call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
         call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

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
      USE param
      USE param1
      use ambm, only: e, w, s, n, t, b
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run
      USE fldvar
      USE physprop
      USE toleranc
      USE fun_avg

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
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
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
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc

! Solids volume fraction at face
      DOUBLE PRECISION :: EPSA

!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

!$omp  parallel default(none) &
!$omp  private( IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
!$omp           VEL,  NUM1, DEN1,EPSA) &
!$omp  shared (ijkstart3, ijkend3, resid_ijk, i_of, j_of, k_of,m,a_m,b_m,w_m,do_k,u_m,v_m,num,den,ncells)
!$omp do reduction(+:num, DEN, NCELLS )
      DO IJK = ijkstart3, ijkend3
        RESID_IJK(IJK) = ZERO

        IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

! Skip walls where some values are undefined.
        IF(WALL_AT(IJK)) cycle

         IF (.NOT.IP_AT_E(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*U_M(IJK)+&
               A_M(IJK,E,M)*U_M(IPJK)+A_M(IJK,W,M)*U_M(IMJK)+&
               A_M(IJK,N,M)*U_M(IJPK)+A_M(IJK,S,M)*U_M(IJMK))
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IJKP = KP_OF(IJK)
               NUM1 = NUM1 - (A_M(IJK,T,M)*U_M(IJKP)+A_M(IJK,B,M)*U_M(IJKM))
            ENDIF

! Ignore momentum residual in stagnant regions.  Need an alternative
! criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2)
            IF (VEL > SMALL_NUMBER) THEN
               NUM1 = ABS(NUM1)
               DEN1 = ABS(A_M(IJK,0,M)*VEL)
! storing value of residual at each ijk location
               RESID_IJK(IJK) = NUM1
! adding to terms that are accumulated
               NCELLS = NCELLS + 1
               NUM = NUM + NUM1
               DEN = DEN + DEN1
            ENDIF
         ENDIF
      ENDDO
!$omp end parallel

      IF(.not.debug_resid) RETURN

! Collecting all the information among all the procesors -
! determining the global sum
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (RESID_IJK( IJK ) > MAX_RESID) THEN
            IJK_RESID = IJK
            MAX_RESID = RESID_IJK( IJK_RESID )
         ENDIF
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

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
      USE param
      USE param1
      use ambm, only: e, w, s, n, t, b
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run
      USE fldvar
      USE physprop
      USE toleranc
      USE fun_avg

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
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
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
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc
! Solids volume fraction at face
      DOUBLE PRECISION :: EPSA

!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

!$omp  parallel default(none) &
!$omp  private( IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
!$omp           VEL,  NUM1, DEN1,EPSA) &
!$omp  shared (ijkstart3, ijkend3, resid_ijk, i_of, j_of, k_of,m,a_m,b_m,w_m,do_k,u_m,v_m,num,den,ncells)
!$omp do reduction(+:num, DEN, NCELLS )
      DO IJK = ijkstart3, ijkend3
        RESID_IJK(IJK) = ZERO

        IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

! Skip walls where some values are undefined.
        IF(WALL_AT(IJK)) cycle


         IF (.NOT.IP_AT_N(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*V_M(IJK)+&
               A_M(IJK,E,M)*V_M(IPJK)+A_M(IJK,W,M)*V_M(IMJK)+&
               A_M(IJK,N,M)*V_M(IJPK)+A_M(IJK,S,M)*V_M(IJMK))
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IJKP = KP_OF(IJK)
               NUM1 = NUM1 - (A_M(IJK,T,M)*V_M(IJKP)+A_M(IJK,B,M)*V_M(IJKM))
            ENDIF

! Ignore momentum residual in stagnant regions.  Need an alternative
! criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2)
            IF (VEL > SMALL_NUMBER) THEN
               NUM1 = ABS(NUM1)
               DEN1 = ABS(A_M(IJK,0,M)*VEL)
! storing value of residual at each ijk location
               RESID_IJK(IJK) = NUM1
! adding to terms that are accumulated
               NCELLS = NCELLS + 1
               NUM = NUM + NUM1
               DEN = DEN + DEN1
            ENDIF
         ENDIF
      ENDDO
!$omp end parallel

      if(.not.debug_resid) return

! Collecting all the information among all the procesors -
! determining the global sum
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
          IF (RESID_IJK( IJK ) > MAX_RESID) THEN
              IJK_RESID = IJK
              MAX_RESID = RESID_IJK( IJK_RESID )
          ENDIF
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determine the maximum among all the procesors
      call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

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
      USE param
      USE param1
      use ambm, only: e, w, s, n, t, b
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE run
      USE fldvar
      USE physprop
      USE toleranc
      USE fun_avg

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
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m(DIMENSION_3, 0:DIMENSION_M)
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
! Numerators and denominators
      DOUBLE PRECISION :: NUM1, DEN1
! Number of fluid cells
      INTEGER :: NCELLS
! New local variables for DMP version
      DOUBLE PRECISION, DIMENSION(ijksize3_all(myPE)) :: RESID_IJK
      DOUBLE PRECISION :: MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER :: IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER :: nproc
! Solids volume fraction at face
      DOUBLE PRECISION :: EPSA

!-----------------------------------------------

! initializing
      NUM = ZERO
      DEN = ZERO
      MAX_RESID = -ONE
      NCELLS = 0

!$omp  parallel default(none) &
!$omp  private( IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
!$omp           VEL,  NUM1, DEN1,EPSA) &
!$omp  shared (ijkstart3, ijkend3, resid_ijk, i_of, j_of, k_of,m,a_m,b_m,w_m,do_k,u_m,v_m,num,den,ncells)
!$omp do reduction(+:num, DEN, NCELLS )
      DO IJK = ijkstart3, ijkend3
        RESID_IJK(IJK) = ZERO

        IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

! Skip walls where some values are undefined.
        IF(WALL_AT(IJK)) cycle


         IF (.NOT.IP_AT_T(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)

! evaluating the residual at cell ijk:
!   RESp = B-sum(Anb*VARnb)-Ap*VARp
!   (where nb = neighbor cells and p = center/0 cell)
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*W_M(IJK)+&
               A_M(IJK,E,M)*W_M(IPJK)+A_M(IJK,W,M)*W_M(IMJK)+&
               A_M(IJK,N,M)*W_M(IJPK)+A_M(IJK,S,M)*W_M(IJMK))
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IJKP = KP_OF(IJK)
               NUM1 = NUM1 - (A_M(IJK,T,M)*W_M(IJKP)+A_M(IJK,B,M)*W_M(IJKM))
            ENDIF

! Ignore momentum residual in stagnant regions.  Need an alternative
! criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2)
            IF (VEL > SMALL_NUMBER) THEN
               NUM1 = ABS(NUM1)
               DEN1 = ABS(A_M(IJK,0,M)*VEL)
! storing value of residual at each ijk location
               RESID_IJK(IJK) = NUM1
! adding to terms that are accumulated
               NCELLS = NCELLS + 1
               NUM = NUM + NUM1
               DEN = DEN + DEN1
            ENDIF
         ENDIF
      ENDDO
!$omp end parallel

      if(.not.debug_resid) return

! Collecting all the information among all the procesors -
! determining the global sum
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

          IF (RESID_IJK( IJK ) > MAX_RESID) THEN
              IJK_RESID = IJK
              MAX_RESID = RESID_IJK( IJK_RESID )
          ENDIF
      ENDDO

! Determining the max residual
      do nproc=0,NumPEs-1
         if(nproc.eq.myPE) then
            MAX_RESID_L(nproc) = MAX_RESID
            IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
         else
            MAX_RESID_L(nproc) = 0.0
            IJK_RESID_L(nproc) = 0
         endif
      enddo

! Determining the maximum among all the procesors
      call global_all_max(MAX_RESID)

! Collecting all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

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
