      SUBROUTINE LEQ_BICGS(VNO, VAR, A_M, B_m, sweep_type, &
                           TOL, pc_type, ITMAX, IER , slo, shi, lo, hi) &
         bind(C, name="mfix_solve_lin_eq")

         use compar, only: mype, pe_io
         use exit_mod, only: mfix_exit
         use matvec_module, only: leq_matvec, leq_residual, leq_scale
         use leqsol, only: icheck_bicgs, minimize_dotproducts
         use leqsol, only: leq_msolve0, leq_msolve1, dot_product_par
         use param1, only: zero, one, small_number

         use solver_params, only: pc_line, pc_none, pc_diag

         use bl_fort_module, only : c_real
         use iso_c_binding , only: c_int

         IMPLICIT NONE

        INTEGER(c_int), INTENT(IN) :: slo(3),shi(3),lo(3),hi(3)
        INTEGER(c_int), INTENT(IN) :: vno

        ! variable
        real(c_real), intent(inout) :: Var&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

        ! Septadiagonal matrix A_m
        real(c_real), INTENT(INOUT) :: A_m&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

        ! Vector b_m
        real(c_real), INTENT(INOUT) :: B_m&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

        ! Sweep direction of leq solver
        !     e.g., options = 'isis', 'rsrs' (default), 'asas'
        integer(c_int), intent(in) :: sweep_type

        ! Type of preconditioner
        integer(c_int), intent(in) :: pc_type

! convergence tolerance (generally leq_tol)
      real(c_real), INTENT(IN) :: TOL

! maximum number of iterations (generally leq_it)
      INTEGER(c_int), INTENT(IN) :: ITMAX

! error indicator
      INTEGER(c_int), INTENT(INOUT) :: IER

! dummy arguments/procedures set as indicated
!     matvec->leq_matvec
! for preconditioner (leq_pc)
!    'line' msolve->leq_msolve  (default)
!    'diag' msolve->leq_msolve1
!    'none' msolve->leq_msolve0
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      INTEGER, PARAMETER :: idebugl = 0
      real(c_real), PARAMETER :: ratiotol = 0.2
      LOGICAL, PARAMETER :: do_unit_scaling = .true.
!-----------------------------------------------
! Local variables
!-----------------------------------------------

      real(c_real), DIMENSION(:,:,:), allocatable :: R,Rtilde, P,Phat, Svec, Shat, Tvec,V

      real(c_real), DIMENSION(0:ITMAX+1) :: alpha, beta, omega, rho
      real(c_real) :: TxS, TxT, RtildexV
      real(c_real) :: Rnorm, Rnorm0, Snorm, TOLMIN, pnorm
      LOGICAL :: isconverged
      INTEGER :: i, ii, jj, kk
      INTEGER :: iter
!-----------------------------------------------

      allocate(R(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)))
      allocate(Rtilde(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)))
      allocate(P(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)))
      allocate(Phat(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)))
      allocate(Svec(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)))
      allocate(Shat(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)))
      allocate(Tvec(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)))
      allocate(V(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)))

! these scalars should not be necessary to initialize but done as failsafe
      rnorm = ZERO
      rnorm0 = ZERO
      snorm = ZERO
      pnorm = ZERO

! initializing
      alpha(:)  = zero
      beta(:)   = zero
      omega(:)  = zero
      rho(:)    = zero

! Adding all by Handan Liu
! zero out R, Rtilde, P, Phat, Svec, Shat, Tvec, V
! --------------------------------
      R(:,:,:)      = zero
      Rtilde(:,:,:) = zero
      P(:,:,:)      = zero
      Phat(:,:,:)   = zero
      Svec(:,:,:)   = zero
      Shat(:,:,:)   = zero
      Tvec(:,:,:)   = zero
      V(:,:,:)      = zero

      TOLMIN = EPSILON( one )

! Scale matrix to have unit diagonal
! ---------------------------------------------------------------->>>
      if (do_unit_scaling) call leq_scale(b_m, A_M, slo, shi, lo, hi)


! Compute initial residual (R = b-A*x) for Ax=b
!    assume initial guess in Var
!    rtilde = r
! ---------------------------------------------------------------->>>
      call leq_residual(b_m, Var, A_M, R, slo, shi)   ! returns R=A*Var

      Rnorm0 = sqrt( dot_product_par( R,R,slo,shi) )

! determine an initial guess for the residual = residual + small random
! number (so it could be set to anything). note that since random_number
! is used to supply the guess, this line could potentially be the source
! of small differences between runs.  the random number is shifted below
! between -1 and 1 and then scaled by factor 1.0D-6*Rnorm0
      !call random_number(Rtilde(:,:,:))
      Rtilde(:,:,:) = R(:,:,:)! + (2.0d0*Rtilde(:,:,:)-1.0d0)*1.0d-6*Rnorm0

      if (idebugl >= 1) then
         if(myPE.eq.0) print*,'leq_bicgs, initial resid ', Rnorm0
      endif
! ----------------------------------------------------------------<<<


! Main loop
! ---------------------------------------------------------------->>>
      iter = 1
      do i=1,itmax

          rho(i-1) = dot_product_par(Rtilde,R,slo,shi)

         if (rho(i-1) .eq. zero) then
            if(i /= 1)then
! Method fails
! --------------------------------
!              print*, 'leq_bicgs: rho(i-1) == 0 '
               ier = -2
            else
! Method converged.  residual is already zero
! --------------------------------
               ier = 0
            endif
            deallocate(R)
            deallocate(Rtilde)
            deallocate(P)
            deallocate(Phat)
            deallocate(Svec)
            deallocate(Shat)
            deallocate(Tvec)
            deallocate(V)
            return
         endif ! rho(i-1).eq.0

         if (i .eq. 1) then
            P(:,:,:) = R(:,:,:)
         else
            beta(i-1) = ( rho(i-1)/rho(i-2) )*( alpha(i-1) / omega(i-1) )
               P(:,:,:) = R(:,:,:) + beta(i-1)*( P(:,:,:) - omega(i-1)*V(:,:,:) )
         endif ! i.eq.1

! Solve A*Phat(:,:,:) = P(:,:,:)
! V(:,:,:) = A*Phat(:,:,:)
! --------------------------------
         if (pc_type.eq.pc_line) then   ! default
            write(6,*) 'No LINE PC support.. How did you get here?'
            stop  23456
         else if (pc_type .eq. pc_none) then
            call LEQ_MSOLVE0(slo, shi, P, Phat) ! returns Phat
         else if (pc_type .eq. pc_diag) then
            call LEQ_MSOLVE1(slo, shi, P, A_m, Phat) ! returns Phat
         end if

         call LEQ_MATVEC(Phat, A_m, V, slo, shi)   ! returns V=A*Phat

         RtildexV = dot_product_par(Rtilde,V,slo,shi)

! compute alpha
! --------------------------------
         alpha(i) = rho(i-1) / RtildexV

! compute Svec
! --------------------------------
         Svec(:,:,:) = R(:,:,:) - alpha(i) * V(:,:,:)


! Solve A*Shat(:) = Svec(:,:,:)
! Tvec(:) = A*Shat(:)
! --------------------------------
         if (pc_type.eq.pc_line) then   ! default
            write(6,*) 'No LINE PC support.. How did you get here?'
            stop  23456
         else if (pc_type .eq. pc_none) then
            call LEQ_MSOLVE0(slo, shi, Svec, Shat)  ! returns Shat
         else if (pc_type .eq. pc_diag) then
            call LEQ_MSOLVE1(slo, shi, Svec, A_m, Shat)  ! returns Shat
         end if


         call LEQ_MATVEC(Shat, A_m, Tvec, slo, shi)   ! returns Tvec=A*Shat

         TxS = dot_product_par(Tvec,Svec,slo,shi)
         TxT = dot_product_par(Tvec,Tvec,slo,shi)

         IF(TxT.eq.Zero) TxT = SMALL_NUMBER

! compute omega
! --------------------------------
         omega(i) = TxS / TxT

! compute new guess for Var
! --------------------------------
         Var(:,:,:) = Var(:,:,:) +  &
            alpha(i)*Phat(:,:,:) + omega(i)*Shat(:,:,:)

         R(:,:,:) = Svec(:,:,:) - omega(i)*Tvec(:,:,:)


! --------------------------------
         if(.not.minimize_dotproducts.or.(mod(iter,icheck_bicgs).eq.0)) then

            Rnorm = sqrt( dot_product_par(R,R,slo,shi) )

            if (idebugl.ge.1) then
               if (myPE.eq.PE_IO) then
                  print*,'     '
                  print*,'     '
                  print*,'iter, Rnorm ', iter, Rnorm
                  print*,'alpha(i), omega(i) ', alpha(i), omega(i)
                  print*,'TxS, TxT ', TxS, TxT
                  print*,'RtildexV, rho(i-1) ', RtildexV, rho(i-1)
               endif
            endif

! Check convergence; continue if necessary
! for continuation, it is necessary that omega(i) .ne. 0
            isconverged = (Rnorm <= TOL*Rnorm0)

            if (isconverged) EXIT
         endif                  ! end if(.not.minimize_dotproducts)

! Advance the iteration count
         iter = iter + 1

      enddo   ! end do i=1,itmax
! end of linear solver loop
! ----------------------------------------------------------------<<<


      if (idebugl >= 1) then
         call LEQ_MATVEC(Var, A_m, R, slo, shi)   ! returns R=A*Var
         do kk = slo(3),shi(3)
            do jj = slo(2),shi(2)
               do ii = slo(1),shi(1)
                  R(ii,jj,kk) = R(ii,jj,kk) - B_M(iI,jJ,kK)
               enddo
            enddo
         enddo

         Rnorm = sqrt( dot_product_par(R,R,slo,shi) )

         if(myPE.eq.0) print*,'leq_bicgs: final Rnorm ', Rnorm

         if(myPE.eq.0)  print*,'leq_bicgs ratio ', iter,  &
         ' L-2', Rnorm/Rnorm0
      endif   ! end if(idebugl >=1)

!      isconverged = (real(Rnorm) <= TOL*Rnorm0);
      if(.NOT.isConverged) isconverged = (real(Rnorm) <= TOL*Rnorm0);
!     write(*,*) '***',iter, isconverged, Rnorm, TOL, Rnorm0, myPE
      IER = 0
      if (.not.isconverged) then
         IER = -1
         if (real(Rnorm) >= ratiotol*real(Rnorm0)) then
            IER = -2
         endif
      endif

      deallocate(R)
      deallocate(Rtilde)
      deallocate(P)
      deallocate(Phat)
      deallocate(Svec)
      deallocate(Shat)
      deallocate(Tvec)
      deallocate(V)

      return
      end subroutine LEQ_BICGS
