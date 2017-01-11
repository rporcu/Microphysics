module leqsol

   use compar, ONLY: istart, iend, jstart, jend, kstart, kend
   use compar, only: iend, jend, kend
   use compar, only: iend1, jend1, kend1
   use compar, only: iend2, jend2, kend2
   use compar, only: iend3, jend3, kend3
   use compar, only: istart, jstart, kstart
   use compar, only: istart1, jstart1, kstart1
   use compar, only: istart2, jstart2, kstart2
   use compar, only: mype
   use error_manager, only: ival, flush_err_msg, err_msg
   use exit_mod, only: mfix_exit
   use functions, only: iplus, jplus, kplus, iminus, jminus, kminus
   use funits, only: dmp_log, unit_log
   use param, only: DIM_EQS
   use param1, only: zero

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

! Maximum number of outer iterations
  INTEGER :: MAX_NIT

! Automatic adjustment of leq parameters possible (set in iterate after
! the completion of first iteration).
  LOGICAL :: LEQ_ADJUST

! Maximum number of linear equation solver iterations
  INTEGER :: LEQ_IT(DIM_EQS)

! Linear equation solver method
  INTEGER :: LEQ_METHOD(DIM_EQS)

! Total Iterations
  INTEGER :: ITER_TOT(DIM_EQS) = 0

! Linear equation solver sweep direction
  CHARACTER(LEN=4) :: LEQ_SWEEP(DIM_EQS)

! Linear equation solver tolerance
  real(c_real) :: LEQ_TOL(DIM_EQS)

! Preconditioner option
  CHARACTER(LEN=4) :: LEQ_PC(DIM_EQS)

! Option to minimize dot products
  LOGICAL :: MINIMIZE_DOTPRODUCTS

! Option to transpose A_m
  LOGICAL :: DO_TRANSPOSE

! Frequency of convergence check in BiCGStab
  INTEGER :: ICHECK_BICGS

! Optimize for massively parallel machine
  LOGICAL :: OPT_PARALLEL

! Linear and non-linear solver statistics
  LOGICAL :: SOLVER_STATISTICS

  LOGICAL :: USE_DOLOOP
  LOGICAL :: IS_SERIAL

CONTAINS

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
  SUBROUTINE REPORT_SOLVER_STATS(TNIT, STEPS)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: TNIT, STEPS

    INTEGER :: LC

    WRITE(ERR_MSG,1100) iVal(TNIT), iVal(TNIT/STEPS)
    CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

1100 FORMAT(/2x,'Total number of non-linear iterations: ', A,/2x,&
         'Average number per time-step: ',A)

    WRITE(ERR_MSG,1200)
    CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

1200 FORMAT(2x,'|',10('-'),'|',13('-'),'|',14('-'),'|',/&
         2x,'| Equation |  Number of  |  Avg Solves  |',/&
         2x,'|  Number  |   Solves    |   for NIT    |',/&
         2x,'|',10('-'),'|',13('-'),'|',14('-'),'|')

    DO LC = 1, DIM_EQS
       WRITE(ERR_MSG,1201) LC, ITER_TOT(LC), ITER_TOT(LC)/TNIT
       CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
    ENDDO

1201 FORMAT(2x,'|',3x,I3,4x,'|',2x,I9,2x,'|',2x,I10,2x,'|',/ &
         2x,'|',10('-'),'|',13('-'),'|',14('-'),'|')


    RETURN
  END SUBROUTINE REPORT_SOLVER_STATS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MSOLVE                                              C
!  Purpose:                                                            C
!  Notes: if leq_method is biggs or cg then this subroutine is         C
!         invoked when leq_pc='line'. if leq_method is gmres then      C
!         this subroutine is invoked (leq_pc setting does not matter)  C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE LEQ_MSOLVE(slo, shi, B_m, A_M, Var, sweep_type)

    use solver_params, only: sweep_rsrs, sweep_isis, sweep_asas
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: slo(3),shi(3)

! Vector b_m
    real(c_real), INTENT(IN) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
! Septadiagonal matrix A_m
    real(c_real), INTENT(IN) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
! Variable
    real(c_real), INTENT(INOUT) :: Var&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
    integer         , intent(in) :: sweep_type
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
    LOGICAL, PARAMETER :: USE_IKLOOP = .FALSE.
    LOGICAL, PARAMETER :: SETGUESS = .TRUE.
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!
    INTEGER :: ITER, NITER
    INTEGER :: I , J, K
    INTEGER :: I1, J1, K1, I2, J2, K2, IK, JK, IJ
    INTEGER :: ISIZE, JSIZE, KSIZE
    INTEGER :: ICASE

    CHARACTER :: CH
    LOGICAL :: DO_ISWEEP
    LOGICAL :: DO_REDBLACK, DO_ALL

!-----------------------------------------------

    IF (SETGUESS) THEN
       do k = slo(3),shi(3)
          do j = slo(2),shi(2)
             do i = slo(1),shi(1)
                VAR(i,j,k) = B_M(I,J,K)
             enddo
          enddo
       enddo
    ENDIF

!   NITER = LEN( CMETHOD )
    NITER = 4

    DO ITER=1,NITER

       DO_ISWEEP   = .FALSE.
       DO_ALL      = .FALSE.
       DO_REDBLACK = .FALSE.

       if ((iter .eq. 1 .or. iter .eq. 3) .and. (sweep_type .eq. sweep_rsrs)) &
         do_redblack = .true.
       if ((iter .eq. 1 .or. iter .eq. 3) .and. (sweep_type .eq. sweep_isis)) &
         do_isweep   = .true.
       if ((iter .eq. 1 .or. iter .eq. 3) .and. (sweep_type .eq. sweep_asas)) &
         do_all      = .true.

       ! Perform sweeps
!      CH = CMETHOD( ITER:ITER )

!      DO_ISWEEP = (CH .EQ. 'I') .OR. (CH .EQ. 'i')
!      DO_JSWEEP = (CH .EQ. 'J') .OR. (CH .EQ. 'j')
!      DO_ALL = (CH .EQ. 'A') .OR. (CH .EQ. 'a')
!      DO_REDBLACK = (CH .EQ. 'R') .OR. (CH .EQ. 'r')
!      DO_SENDRECV = (CH .EQ. 'S') .OR. (CH .EQ. 's')

! do_all true only for leq_pc='asas'
! ---------------------------------------------------------------->>>
       IF(DO_ALL) THEN        ! redblack for all sweeps, not used by default
! JK Loop
! --------------------------------
          j1 = jstart
          k1 = kstart
          j2 = jend
          k2 = kend
          jsize = j2-j1+1
          ksize = k2-k1+1
          DO icase = 1, 2
             DO JK=icase, ksize*jsize, 2
                if (mod(jk,jsize).ne.0) then
                   k = int( jk/jsize ) + k1
                else
                   k = int( jk/jsize ) + k1 -1
                endif
                j = (jk-1-(k-k1)*jsize) + j1 + mod(k,2)
                if(j.gt.j2) j=j-j2 + j1 -1
                CALL LEQ_JKSWEEP(J, K, Var, A_m, B_m, slo, shi)
             ENDDO
          ENDDO

! IJ Loop
! --------------------------------
          i1 = istart
          j1 = jstart
          i2 = iend
          j2 = jend
          isize = i2-i1+1
          jsize = j2-j1+1
          DO icase = 1, 2
             DO IJ=icase, jsize*isize, 2
                if (mod(ij,isize).ne.0) then
                   j = int( ij/isize ) + j1
                else
                   j = int( ij/isize ) + j1 -1
                endif
                i = (ij-1-(j-j1)*isize) + i1 + mod(j,2)
                if(i.gt.i2) i=i-i2 + i1 -1
                CALL LEQ_IJSWEEP(I, J, Var, A_m, B_m, slo, shi)
             ENDDO

          ENDDO

! IK Loop
! --------------------------------
          i1 = istart
          k1 = kstart
          i2 = iend
          k2 = kend
          isize = i2-i1+1
          ksize = k2-k1+1

          DO icase = 1, 2
             DO IK=icase, ksize*isize, 2
                if (mod(ik,isize).ne.0) then
                   k = int( ik/isize ) + k1
                else
                   k = int( ik/isize ) + k1 -1
                endif
                i = (ik-1-(k-k1)*isize) + i1 + mod(k,2)
                if(i.gt.i2) i=i-i2 + i1 -1
                CALL LEQ_IKSWEEP(I, K, Var, A_m, B_m, slo, shi)
             ENDDO

          ENDDO
       ENDIF ! end DO_ALL
! ----------------------------------------------------------------<<<

! do_redblack only true leq_pc='rsrs'
! ---------------------------------------------------------------->>>
       IF(DO_REDBLACK) THEN
          DO k=kstart,kend
             IF(mod(k,2).ne.0)THEN
                DO I=istart+1,iend,2
                   CALL LEQ_IKSWEEP(I, K, Var, A_m, B_m, slo, shi)
                ENDDO
             ELSE
                DO I=istart,iend,2
                   CALL LEQ_IKSWEEP(I, K, Var, A_m, B_m, slo, shi)
                ENDDO
             ENDIF
          ENDDO
          DO k=kstart,kend
             IF(mod(k,2).ne.0)THEN
                DO I=istart,iend,2
                   CALL LEQ_IKSWEEP(I, K, Var, A_m, B_m, slo, shi)
                ENDDO
             ELSE
                DO I=istart+1,iend,2
                   CALL LEQ_IKSWEEP(I, K, Var, A_m, B_m, slo, shi)
                ENDDO
             ENDIF
          ENDDO

       ENDIF       ! end if(do_redblack)
! ----------------------------------------------------------------<<<

!  Not sure the purpose of us_ikloop
!  The SMP directives below need review                        !Tingwen Jan 2012
! use_ikloop is currently hard-wired to false (so goto else branch)
! ---------------------------------------------------------------->>>
       IF(USE_IKLOOP) THEN
          i1 = istart
          k1 = kstart
          i2 = iend
          k2 = kend
          isize = i2-i1+1
          ksize = k2-k1+1
          IF (DO_ISWEEP) THEN
             DO IK=1, ksize*isize
                if (mod(ik,isize).ne.0) then
                   k = int( ik/isize ) + k1
                else
                   k = int( ik/isize ) + k1 -1
                endif
                i = (ik-1-(k-k1)*isize) + i1
                CALL LEQ_IKSWEEP(I, K, Var, A_m, B_m, slo, shi)
             ENDDO
          ENDIF

! ----------------------------------------------------------------<<<
       ELSE   ! else branch of if(use_ikloop)
!  Not sure the purpose of us_ikloop
!  The SMP directives below need review                        !Tingwen Jan 2012
! ---------------------------------------------------------------->>>
          IF (DO_ISWEEP) THEN
             DO K=kstart,kend
                DO I=istart,iend
                   CALL LEQ_IKSWEEP(I, K, Var, A_m, B_m, slo, shi)
                ENDDO
             ENDDO
          ENDIF
       ENDIF   ! end if/else (use(ikloop)
    ENDDO   ! end do iter=1,niter

  END SUBROUTINE LEQ_MSOLVE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MSOLVE0                                             C
!  Notes: do nothing or no preconditioning (leq_pc='none')             C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE LEQ_MSOLVE0(slo, shi, B_m, A_M, Var, sweep_type)

    IMPLICIT NONE

  integer, intent(in) :: slo(3),shi(3)

! Vector b_m
    real(c_real), INTENT(IN) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
! Septadiagonal matrix A_m
    real(c_real), INTENT(IN) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
! Variable
    real(c_real), INTENT(OUT) :: Var&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    ! sweep direction
    integer         , intent(in) :: sweep_type
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    integer :: i,j,k
!-----------------------------------------------

! do nothing or no preconditioning
      DO K = slo(3),shi(3)
        DO J = slo(2),shi(2)
          DO I = slo(1),shi(1)
          var(i,j,k) = b_m(i,j,k)
       enddo
       enddo
       enddo

  end subroutine leq_msolve0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_MSOLVE1                                             C
!  Notes: diagonal scaling (leq_pc='diag')                             C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE LEQ_MSOLVE1(slo, shi, B_m, A_M, Var, sweep_type)

    IMPLICIT NONE

    integer, intent(in) :: slo(3),shi(3)

    ! Vector b_m
    real(c_real), INTENT(IN) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    ! Septadiagonal matrix A_m
    real(c_real), INTENT(IN) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

    real(c_real), INTENT(OUT) :: Var&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    ! sweep direction
    integer         , intent(in) :: sweep_type

    integer :: i,j,k

    do k = slo(3),shi(3)
       do i = slo(1),shi(1)
          do j = slo(2),shi(2)
             var(i,j,k) = b_m(i,j,k)/A_m(i,j,k,0)
          enddo
       enddo
    enddo

  end subroutine leq_msolve1

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_IKSWEEP                                             C
!  Purpose: Perform line sweep at coordinate I, K                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_IKSWEEP(I, K, VAR, A_M, B_M, slo, shi)

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Line position
      INTEGER, INTENT(IN) :: I, K
      INTEGER, INTENT(IN) :: slo(3), shi(3)
! Variable
      real(c_real), INTENT(INOUT) :: Var&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
! Septadiagonal matrix A_m
      real(c_real), INTENT(IN) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
! Vector b_m
      real(c_real), INTENT(IN) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      real(c_real), DIMENSION(slo(2):shi(2)) :: CC, DD, EE, BB
      INTEGER :: NSTART, NEND, INFO
      INTEGER :: J
!-----------------------------------------------


      NSTART = slo(2)+2
      NEND = shi(2)-2

      DO J = slo(2), shi(2)
         DD(J) = A_M(I,J,K,  0)
         CC(J) = A_M(I,J,K, -2)
         EE(J) = A_M(I,J,K,  2)
         BB(J) = B_M(I,J,K) -  A_M(I,J,K,-1) * Var(iminus(i,j,k),j,k) &
                            -  A_M(I,J,K, 1) * Var( iplus(i,j,k),j,k) &
                            -  A_M(I,J,K,-3) * Var(i,j,kminus(i,j,k)) &
                            -  A_M(I,J,K, 3) * Var(i,j, kplus(i,j,k))
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)

      IF (INFO.NE.0) THEN
         write(*,*) 'leq_iksweep',INFO, myPE
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'ROUTINE = ', ' IKSWEEP'
         RETURN
      ENDIF

      DO J = slo(2)+2, shi(2)-2
         Var(i,j,k) = BB(J)
      ENDDO

      END SUBROUTINE LEQ_IKSWEEP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_JKSWEEP                                             C
!  Purpose: Perform line sweep at coordinate I, K                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_JKSWEEP(J, K, VAR, A_M, B_M, slo, shi)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: slo(3), shi(3)
! Line position
      INTEGER, INTENT(IN) :: J, K
! Variable
      real(c_real), INTENT(INOUT) :: Var&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
! Septadiagonal matrix A_m
      real(c_real), INTENT(IN) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
! Vector b_m
      real(c_real), INTENT(IN) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      real(c_real), DIMENSION (ISTART:IEND) :: CC, DD, EE, BB
      INTEGER :: NSTART, NEND, INFO, I
!-----------------------------------------------

      NEND = IEND
      NSTART = ISTART

      DO I=NSTART,NEND
         DD(I) = A_M(I,J,K,  0)
         CC(I) = A_M(I,J,K, -1)
         EE(I) = A_M(I,J,K,  1)
         BB(I) = B_M(I,J,K)  -  A_M(I,J,K,-2) * Var(i,jminus(i,j,k),k) &
                             -  A_M(I,J,K, 2) * Var(i,jplus(i,j,k) ,k) &
                             -  A_M(I,J,K,-3) * Var(i,j,kminus(i,j,k)) &
                             -  A_M(I,J,K, 3) * Var(i,j ,kplus(i,j,k))
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)

      IF (INFO.NE.0) THEN
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'ROUTINE = ', ' JKSWEEP'
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'DGTSV RETURNS INFO = ', INFO
         call mfix_exit(myPE)
      ENDIF

      DO I=NSTART, NEND
         Var(i,j,k) = BB(I)
      ENDDO

      RETURN
      END SUBROUTINE LEQ_JKSWEEP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_IJSWEEP                                             C
!  Purpose: Perform line sweep at coordinate I, K                      C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_IJSWEEP(I, J, VAR, A_M, B_M, slo, shi)

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Line position
      INTEGER, INTENT(IN) :: I, J
      INTEGER, INTENT(IN) :: slo(3),shi(3)
! Variable
      real(c_real), INTENT(INOUT) :: Var&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
! Septadiagonal matrix A_m
      real(c_real), INTENT(IN) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
! Vector b_m
      real(c_real), INTENT(IN) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      real(c_real), DIMENSION (KSTART:KEND) :: CC, DD, EE, BB
      INTEGER :: NEND, NSTART, INFO, K
!-----------------------------------------------

      NEND = KEND
      NSTART = KSTART

      DO K=NSTART, NEND
         DD(K) = A_M(I,J,K,  0)
         CC(K) = A_M(I,J,K, -3)
         EE(K) = A_M(I,J,K,  3)
         BB(K) = B_M(I,J,K)  -  A_M(I,J,K,-2) * Var(i,jminus(i,j,k),k) &
                             -  A_M(I,J,K, 2) * Var(i, jplus(i,j,k),k) &
                             -  A_M(I,J,K,-1) * Var(iminus(i,j,k),j,k) &
                             -  A_M(I,J,K, 1) * Var( iplus(i,j,k),j,k)
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)

      DO K=NSTART, NEND
         Var(i,j,k) = BB(k)
      ENDDO

      END SUBROUTINE LEQ_IJSWEEP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  real(c_real) function dot_product_par(r1,r2,slo,shi)

    use bl_fort_module, only : c_real

    implicit none
    integer, intent(in) :: slo(3),shi(3)

    real(c_real), intent(in) :: r1&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(c_real), intent(in) :: r2&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    real(c_real) :: prod
    integer :: i, j, k

    prod = 0.0d0

    do k = slo(3),shi(3)
       do i = slo(1),shi(1)
          do j = slo(2),shi(2)
               prod = prod + r1(i,j,k)*r2(i,j,k)
          enddo
       enddo
    enddo

    dot_product_par = prod

  end function dot_product_par

end module leqsol
