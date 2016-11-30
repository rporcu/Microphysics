MODULE leqsol

  use param, only: DIM_EQS

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
  DOUBLE PRECISION :: LEQ_TOL(DIM_EQS)

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

    use error_manager, only: ival, flush_err_msg, err_msg

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
!  Subroutine: LEQ_MATVEC                                              C
!  Purpose: Compute matrix vector multiplication                       C
!           (for linear equation Ax=b compute Ax                       C
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

  SUBROUTINE LEQ_MATVEC(VNAME, VAR, A_M, Avar)

!-----------------------------------------------
! Modules
!-----------------------------------------------
    USE compar, ONLY: istart, iend, jstart, jend, kstart, kend
    USE geometry, ONLY: do_k
    USE param, ONLY: DIMENSION_3
      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3
    IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable name
    CHARACTER(LEN=*), INTENT(IN) :: Vname
! Variable
    DOUBLE PRECISION, INTENT(IN) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
    DOUBLE PRECISION, INTENT(IN) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector AVar
    DOUBLE PRECISION, INTENT(OUT) :: AVar(DIMENSION_3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Variable
    INTEGER :: I, J, K, IJK
!-----------------------------------------------

      if (do_k) then
         do k = kstart,kend
            do i = istart,iend
               do j = jstart,jend
                  ijk = funijk(i,j,k)
                  AVar(ijk) =  A_m(i,j,k,-3) * Var(funijk(i,j,kminus(i,j,k)))   &
                             + A_m(i,j,k,-2) * Var(funijk(i,jminus(i,j,k),k))   &
                             + A_m(i,j,k,-1) * Var(funijk(iminus(i,j,k),j,k))   &
                             + A_m(i,j,k, 0) * Var(funijk(i,j,k)            )   &
                             + A_m(i,j,k, 1) * Var(funijk(iplus(i,j,k),j,k) )   &
                             + A_m(i,j,k, 2) * Var(funijk(i,jplus(i,j,k),k) )   &
                             + A_m(i,j,k, 3) * Var(funijk(i,j,kplus(i,j,k)) )
               enddo
            enddo
         enddo

      else
         k = 1
         do i = istart,iend
            do j = jstart,jend
               ijk = funijk(i,j,k)
               AVar(ijk) =  A_m(i,j,k,-2) * Var(funijk(i,jminus(i,j,k),k))   &
                          + A_m(i,j,k,-1) * Var(funijk(iminus(i,j,k),j,k))   &
                          + A_m(i,j,k, 0) * Var(funijk(i,j,k)            )   &
                          + A_m(i,j,k, 1) * Var(funijk(iplus(i,j,k),j,k) )   &
                          + A_m(i,j,k, 2) * Var(funijk(i,jplus(i,j,k),k) )
            enddo
         enddo

      endif


       ! call send_recv(Avar,nlayers_bicgs)
    RETURN

  CONTAINS

    INCLUDE 'functions.inc'

  END SUBROUTINE LEQ_MATVEC


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

  SUBROUTINE LEQ_MSOLVE(VNAME, B_m, A_M, Var, CMETHOD)

!-----------------------------------------------
! Modules
!-----------------------------------------------
    USE param
    USE param1
    USE geometry
    USE compar
    USE functions
    IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable name
    CHARACTER(LEN=*), INTENT(IN) :: Vname
! Vector b_m
    DOUBLE PRECISION, INTENT(IN) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Septadiagonal matrix A_m
    DOUBLE PRECISION, INTENT(IN) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Variable
    DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
    CHARACTER(LEN=4), INTENT(IN) :: CMETHOD
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
    INTEGER :: IJK, I , J, K
    INTEGER :: I1, J1, K1, I2, J2, K2, IK, JK, IJ
    INTEGER :: ISIZE, JSIZE, KSIZE
    INTEGER :: ICASE

!     CHARACTER(LEN=4), PARAMETER :: CMETHOD = 'II'
    CHARACTER :: CH
    LOGICAL :: DO_ISWEEP, DO_JSWEEP, DO_KSWEEP
    LOGICAL :: DO_SENDRECV, DO_REDBLACK, DO_ALL

!-----------------------------------------------

    IF (SETGUESS) THEN
       do k = kstart3,kend3
          do i = istart3,iend3
             do j = jstart3,jend3
                IJK = funijk(i,j,k)
                VAR(IJK) = B_M(I,J,K)
             enddo
          enddo
       enddo

       ! call send_recv(var,nlayers_bicgs)
    ENDIF

    NITER = LEN( CMETHOD )

    DO ITER=1,NITER

! Perform sweeps
       CH = CMETHOD( ITER:ITER )
       DO_ISWEEP = (CH .EQ. 'I') .OR. (CH .EQ. 'i')
       DO_JSWEEP = (CH .EQ. 'J') .OR. (CH .EQ. 'j')
       DO_KSWEEP = (CH .EQ. 'K') .OR. (CH .EQ. 'k')
       DO_ALL = (CH .EQ. 'A') .OR. (CH .EQ. 'a')
       DO_REDBLACK = (CH .EQ. 'R') .OR. (CH .EQ. 'r')
       DO_SENDRECV = (CH .EQ. 'S') .OR. (CH .EQ. 's')

       IF (NO_K) THEN   ! two dimensional
! 2D run no need to enable openmp parallel
          IF ( DO_ISWEEP ) THEN
             DO I=istart,iend,1
                CALL LEQ_ISWEEP( I, Vname, Var, A_m, B_m )
             ENDDO
          ENDIF
! ----------------------------------------------------------------<<<
! Handan Liu added 2D RSRS sweep and parallelized this loop on Jan 22 2013:
          IF (DO_REDBLACK) THEN
             DO I=istart,iend,2
                CALL LEQ_ISWEEP( I, Vname, Var, A_m, B_m )
             ENDDO
             DO I=istart+1,iend,2
                CALL LEQ_ISWEEP( I, Vname, Var, A_m, B_m )
             ENDDO
          ENDIF
! ---------------------------------------------------------------->>>
       ELSE   ! three dimensional


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
                   CALL LEQ_JKSWEEP(J, K, Vname, Var, A_m, B_m)
                ENDDO

             ENDDO
             ! call send_recv(var,nlayers_bicgs)

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
                   CALL LEQ_IJSWEEP(I, J, Vname, Var, A_m, B_m)
                ENDDO

             ENDDO
             ! call send_recv(var,nlayers_bicgs)

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
                   CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                ENDDO

             ENDDO
          ENDIF ! end DO_ALL
! ----------------------------------------------------------------<<<

! do_redblack only true leq_pc='rsrs'
! ---------------------------------------------------------------->>>
          IF(DO_REDBLACK) THEN
!i1 = istart
!k1 = kstart
!i2 = iend
!k2 = kend
!isize = i2-i1+1
!ksize = k2-k1+1
!               DO icase = 1, 2
!                  DO IK=icase, ksize*isize, 2
!                     if (mod(ik,isize).ne.0) then
!                        k = int( ik/isize ) + k1
!                     else
!                        k = int( ik/isize ) + k1 -1
!                     endif
!                     i = (ik-1-(k-k1)*isize) + i1 + mod(k,2)
!                     if(i.gt.i2) i=i-i2 + i1 -1
!                     CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
!                  ENDDO
!               ENDDO
!             ELSE
! Handan Liu split above loop for OpenMP at May 22 2013, modified at July 17
             DO k=kstart,kend
                IF(mod(k,2).ne.0)THEN
                   DO I=istart+1,iend,2
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ELSE
                   DO I=istart,iend,2
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ENDIF
             ENDDO
             DO k=kstart,kend
                IF(mod(k,2).ne.0)THEN
                   DO I=istart,iend,2
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ELSE
                   DO I=istart+1,iend,2
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ENDIF
             ENDDO

          ENDIF       ! end if(do_redblack)
! ----------------------------------------------------------------<<<

!  Not sure the purpose of us_ikloop
!  The SMP directives below need review                        !Tingwen Jan 2012
          IF(USE_IKLOOP) THEN
! use_ikloop is currently hard-wired to false (so goto else branch)
! ---------------------------------------------------------------->>>
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
                   CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                ENDDO
             ENDIF
             IF (DO_KSWEEP) THEN
                DO IK=1, ksize*isize
                   if (mod(ik,ksize).ne.0) then
                      i = int( ik/ksize ) + i1
                   else
                      i = int( ik/ksize ) + i1 -1
                   endif
                   k = (ik-1-(i-i1)*ksize) + k1
                   CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
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
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ENDDO
             ENDIF
             IF (DO_KSWEEP) THEN
                DO I=istart,iend
                   DO K=kstart,kend
                      CALL LEQ_IKSWEEP(I, K, Vname, Var, A_m, B_m)
                   ENDDO
                ENDDO
             ENDIF
          ENDIF   ! end if/else (use(ikloop)
! ----------------------------------------------------------------<<<

       ENDIF   ! end if/else (do_k)


! this is called for all settings of leq_pc
       ! IF (DO_SENDRECV) call send_recv(var,nlayers_bicgs)


    ENDDO   ! end do iter=1,niter

    RETURN
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

  SUBROUTINE LEQ_MSOLVE0(VNAME, B_m, A_M, Var, CMETHOD)

!-----------------------------------------------
! Modules
!-----------------------------------------------
    USE param
    USE param1
    USE geometry
    USE compar
    USE functions
    IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable name
    CHARACTER(LEN=*), INTENT(IN) :: Vname
! Vector b_m
    DOUBLE PRECISION, INTENT(IN) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Septadiagonal matrix A_m
    DOUBLE PRECISION, INTENT(IN) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Variable
    DOUBLE PRECISION, INTENT(OUT) :: Var(DIMENSION_3)
! sweep direction
    CHARACTER(LEN=4), INTENT(IN) :: CMETHOD
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    integer :: ijk,i,j,k
!-----------------------------------------------

! do nothing or no preconditioning
      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
          var(ijk) = b_m(i,j,k)
       enddo
       enddo
       enddo
    ! call send_recv(var,nlayers_bicgs)

    return
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

  SUBROUTINE LEQ_MSOLVE1(VNAME, B_m, A_M, Var, CMETHOD)

!-----------------------------------------------
! Modules
!-----------------------------------------------
    USE param
    USE param1
    USE geometry
    USE compar
    USE functions
    IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Variable name
    CHARACTER(LEN=*), INTENT(IN) :: Vname
! Vector b_m
    DOUBLE PRECISION, INTENT(IN) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Septadiagonal matrix A_m
    DOUBLE PRECISION, INTENT(IN) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Variable
    DOUBLE PRECISION, INTENT(OUT) :: Var(DIMENSION_3)
! sweep direction
    CHARACTER(LEN=4), INTENT(IN) :: CMETHOD
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    integer :: i,j,k, ijk
!-----------------------------------------------

       do k=kstart2,kend2
          do i=istart2,iend2
             do j=jstart2,jend2
                ijk = funijk( i,j,k )
                var(ijk) = b_m(i,j,k)/A_m(i,j,k,0)
             enddo
          enddo
       enddo

    ! call send_recv(var,nlayers_bicgs)

    return
  end subroutine leq_msolve1


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_ISWEEP                                              C
!  Purpose: Perform line sweep at coordinate I                         C
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

      SUBROUTINE LEQ_ISWEEP(I, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE funits
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
!  Line position
      INTEGER, INTENT(IN) :: I
! Variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION (JSTART:JEND) :: CC, DD, EE, BB
      INTEGER :: NSTART, NEND, INFO
      INTEGER :: IJK, J, K, IM1JK, IP1JK
!-----------------------------------------------

      NEND = JEND
      NSTART = JSTART
      K = 1

      DO J=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         IM1JK = funijk(iminus(i,j,k),j,k)
         IP1JK = funijk(iplus(i,j,k),j,k)
         DD(J) = A_M(I,J,K,  0)
         CC(J) = A_M(I,J,K, -2)
         EE(J) = A_M(I,J,K,  2)
         BB(J) = B_M(I,J,K) -  A_M(I,J,K,-1) * Var( IM1JK )  &
                          -  A_M(I,J,K, 1) * Var( IP1JK )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
!     CALL DGTSL( JEND-JSTART+1, CC, DD, EE, BB, INFO )
      CALL DGTSV( JEND-JSTART+1, 1, CC(JSTART+1), DD, EE, BB, JEND-JSTART+1, INFO )

      IF (INFO.NE.0) THEN
         RETURN
      ENDIF

      DO J=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         Var(IJK) =  BB(J)
      ENDDO

      RETURN
      END SUBROUTINE LEQ_ISWEEP

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

      SUBROUTINE LEQ_IKSWEEP(I, K, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE funits
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Line position
      INTEGER, INTENT(IN) :: I, K
! Variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION(JSTART:JEND) :: CC, DD, EE, BB
      INTEGER :: NSTART, NEND, INFO
      INTEGER :: IJK, J
!-----------------------------------------------

      NEND = JEND
      NSTART = JSTART

      DO J=NSTART, NEND
         IJK = funijk(i,j,k)
         DD(J) = A_M(I,J,K,  0)
         CC(J) = A_M(I,J,K, -2)
         EE(J) = A_M(I,J,K,  2)
         BB(J) = B_M(I,J,K) -  A_M(I,J,K,-1) * Var( funijk(iminus(i,j,k),j,k) ) &
                          -  A_M(I,J,K, 1) * Var( funijk(iplus(i,j,k),j,k) ) &
                          -  A_M(I,J,K,-3) * Var( funijk(i,j,kminus(i,j,k)) ) &
                          -  A_M(I,J,K, 3) * Var( funijk(i,j,kplus(i,j,k)) )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
!     CALL DGTSL( JEND-JSTART+1, CC, DD, EE, BB, INFO )
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)

      IF (INFO.NE.0) THEN
         write(*,*) 'leq_iksweep',INFO, myPE
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'ROUTINE = ', ' IKSWEEP'
         RETURN
      ENDIF

      DO J=NSTART, NEND
         Var(funijk(i,j,k)) = BB(J)
      ENDDO

      RETURN
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

      SUBROUTINE LEQ_JKSWEEP(J, K, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE funits
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Line position
      INTEGER, INTENT(IN) :: J, K
! Variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION (ISTART:IEND) :: CC, DD, EE, BB
      INTEGER :: NSTART, NEND, INFO, IJK, I
!-----------------------------------------------

      NEND = IEND
      NSTART = ISTART

      DO I=NSTART,NEND
         IJK = FUNIJK(I,J,K)
         DD(I) = A_M(I,J,K,  0)
         CC(I) = A_M(I,J,K, -1)
         EE(I) = A_M(I,J,K,  1)
         BB(I) = B_M(I,J,K)    -  A_M(I,J,K,-2) * Var( funijk(i,jminus(i,j,k),k) ) &
                             -  A_M(I,J,K, 2) * Var( funijk(i,jplus(i,j,k) ,k) ) &
                             -  A_M(I,J,K,-3) * Var( funijk(i,j,kminus(i,j,k)) ) &
                             -  A_M(I,J,K, 3) * Var( funijk(i,j ,kplus(i,j,k))  )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)

      IF (INFO.NE.0) THEN
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'VNAME = ', VNAME
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'ROUTINE = ', ' JKSWEEP'
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'DGTSV RETURNS INFO = ', INFO
         call mfix_exit(myPE)
      ENDIF

      DO I=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         Var(IJK) = BB(I)
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

      SUBROUTINE LEQ_IJSWEEP(I, J, Vname, VAR, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE funits
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Line position
      INTEGER, INTENT(IN) :: I, J
! Variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! Variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(IN) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION (KSTART:KEND) :: CC, DD, EE, BB
      INTEGER :: NEND, NSTART, INFO, IJK, K
!-----------------------------------------------

      NEND = KEND
      NSTART = KSTART

      DO K=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         DD(K) = A_M(I,J,K,  0)
         CC(K) = A_M(I,J,K, -3)
         EE(K) = A_M(I,J,K,  3)
         BB(K) = B_M(I,J,K)    -  A_M(I,J,K,-2) * Var( funijk(i,jminus(i,j,k),k) ) &
                             -  A_M(I,J,K, 2) * Var( funijk(i,jplus(i,j,k),k) ) &
                             -  A_M(I,J,K,-1) * Var( funijk(iminus(i,j,k),j,k) ) &
                             -  A_M(I,J,K, 1) * Var( funijk(iplus(i,j,k),j,k) )
      ENDDO

      CC(NSTART) = ZERO
      EE(NEND) = ZERO
      INFO = 0
      CALL DGTSV(NEND-NSTART+1, 1, CC(NSTART+1), DD, EE, BB, NEND-NSTART+1, INFO)

      IF (INFO.NE.0) THEN
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'VNAME = ', VNAME
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'ROUTINE = ', ' IJSWEEP'
         IF(DMP_LOG)WRITE (UNIT_LOG,*) 'DGTSV RETURNS INFO = ', INFO
         call mfix_exit(myPE)
      ENDIF

      DO K=NSTART, NEND
         IJK = FUNIJK(I,J,K)
         Var(IJK) = BB(K)
      ENDDO

      RETURN
      END SUBROUTINE LEQ_IJSWEEP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  double precision function dot_product_par(r1,r2)

!-----------------------------------------------
! Modules
!-----------------------------------------------
    use geometry
    use compar
    use functions
    implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
!      double precision, intent(in), dimension(DIMENSION_3) :: r1,r2
    double precision, intent(in), dimension(DIMENSION_3) :: r1,r2
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    double precision :: prod
    integer :: i, j, k, ijk
!-----------------------------------------------

       prod = 0.0d0

          do k = kstart1, kend1
             do i = istart1, iend1
                do j = jstart1, jend1
                   ijk = funijk_map_c (i,j,k)
                   prod = prod + r1(ijk)*r2(ijk)
                enddo
             enddo
          enddo

          dot_product_par = prod

  end function dot_product_par


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  function dot_product_par2(r1,r2,r3,r4)

!-----------------------------------------------
! Modules
!-----------------------------------------------
    use geometry
    use compar
    use functions
    implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
    double precision, intent(in), dimension(DIMENSION_3) :: r1,r2,r3,r4
!-----------------------------------------------
! Local variables
!-----------------------------------------------
    double precision, Dimension(2) :: prod, dot_product_par2
    integer :: i, j, k, ijk
!-----------------------------------------------

       prod(:) = 0.0d0

       do k = kstart1, kend1
          do i = istart1, iend1
             do j = jstart1, jend1

                ijk = funijk_map_c (i,j,k)
                prod(1) = prod(1) + r1(ijk)*r2(ijk)
                prod(2) = prod(2) + r3(ijk)*r4(ijk)
             enddo
          enddo
       enddo

       dot_product_par2 = prod

  end function dot_product_par2

END MODULE leqsol
