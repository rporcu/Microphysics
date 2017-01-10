      SUBROUTINE LEQ_BICGS(VNO, VAR, A_M, B_m, sweep_type, &
                           TOL, pc_type, ITMAX, IER , slo, shi, lo, hi) &
         bind(C, name="mfix_solve_lin_eq")

         use compar, only: istart2, iend2
         use compar, only: jstart2, jend2
         use compar, only: kstart2, kend2
         use compar, only: mype, pe_io, numpes
         use exit_mod, only: mfix_exit
         use functions, only: funijk
         use matvec_module, only: leq_matvec
         use leqsol, only: is_serial, icheck_bicgs, minimize_dotproducts
         use leqsol, only: leq_msolve, leq_msolve0, leq_msolve1, dot_product_par, dot_product_par2, iter_tot
         use param, only: dimension_3
         use param1, only: zero, one, small_number

         use solver_params, only: pc_line, pc_none, pc_diag

         use bl_fort_module, only : c_real
         use iso_c_binding , only: c_int

         IMPLICIT NONE

        INTEGER, INTENT(IN) :: slo(3),shi(3),lo(3),hi(3)
        INTEGER, INTENT(IN) :: vno

        ! variable
        real(c_real), DIMENSION(DIMENSION_3), INTENT(INOUT) :: Var

        ! Septadiagonal matrix A_m
        real(c_real), INTENT(INOUT) :: A_m&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

        ! Vector b_m
        real(c_real), INTENT(INOUT) :: B_m&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

        ! Sweep direction of leq solver
        !     e.g., options = 'isis', 'rsrs' (default), 'asas'
        integer         , intent(in) :: sweep_type

        ! Type of preconditioner
        integer         , intent(in) :: pc_type

! convergence tolerance (generally leq_tol)
      real(c_real), INTENT(IN) :: TOL

! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX

! error indicator
      INTEGER, INTENT(INOUT) :: IER

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

      real(c_real), DIMENSION(:), allocatable :: R,Rtilde, P,Phat, Svec, Shat, Tvec,V

      real(c_real), DIMENSION(0:ITMAX+1) :: &
                        alpha, beta, omega, rho
      real(c_real) :: TxS, TxT, RtildexV, &
                          aijmax, oam
      real(c_real) :: Rnorm, Rnorm0, Snorm, TOLMIN, pnorm
      LOGICAL :: isconverged
      INTEGER :: i, j, k, ijk, ii, jj, kk
      INTEGER :: iter
      real(c_real), DIMENSION(2) :: TxS_TxT
!-----------------------------------------------

      allocate(R(DIMENSION_3))
      allocate(Rtilde(DIMENSION_3))
      allocate(P(DIMENSION_3))
      allocate(Phat(DIMENSION_3))
      allocate(Svec(DIMENSION_3))
      allocate(Shat(DIMENSION_3))
      allocate(Tvec(DIMENSION_3))
      allocate(V(DIMENSION_3))

      is_serial = numPEs.eq.1.and.is_serial

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
         R(:) = zero
         Rtilde(:) = zero
         P(:) = zero
         Phat(:) = zero
         Svec(:) = zero
         Shat(:) = zero
         Tvec(:) = zero
         V(:) = zero

      TOLMIN = EPSILON( one )

! Scale matrix to have unit diagonal
! ---------------------------------------------------------------->>>
      if (do_unit_scaling) then

         do k = lo(3)-1, hi(3)+1
            do i = lo(1)-1, hi(1)+1
               do j = lo(2)-1, hi(2)+1
                  IJK = funijk(i,j,k)
                  aijmax = maxval(abs(A_M(i,j,k,:)) )
                  OAM = one/aijmax
                  A_M(I,J,K,:) = A_M(I,J,K,:)*OAM
                  B_M(I,J,K) = B_M(I,J,K)*OAM
               enddo
            enddo
         enddo
      endif
! ----------------------------------------------------------------<<<


! Compute initial residual (R = b-A*x) for Ax=b
!    assume initial guess in Var
!    rtilde = r
! ---------------------------------------------------------------->>>
      call LEQ_MATVEC(B_M, Var, A_M, R, slo, shi, lo, hi)   ! returns R=A*Var

      if(is_serial) then
         Rnorm0 = zero
         Rnorm0 = dot_product(R,R)
         Rnorm0 = sqrt( Rnorm0 )
      else
         Rnorm0 = sqrt( dot_product_par( R, R ) )
      endif

! determine an initial guess for the residual = residual + small random
! number (so it could be set to anything). note that since random_number
! is used to supply the guess, this line could potentially be the source
! of small differences between runs.  the random number is shifted below
! between -1 and 1 and then scaled by factor 1.0D-6*Rnorm0
      call random_number(Rtilde(:))
      Rtilde(:) = R(:) + (2.0d0*Rtilde(:)-1.0d0)*1.0d-6*Rnorm0

      if (idebugl >= 1) then
         if(myPE.eq.0) print*,'leq_bicgs, initial resid ', Rnorm0
      endif
! ----------------------------------------------------------------<<<


! Main loop
! ---------------------------------------------------------------->>>
      iter = 1
      do i=1,itmax

         if(is_serial) then
            rho(i-1) = dot_product( Rtilde, R )
         else
            rho(i-1) = dot_product_par( Rtilde, R )
         endif ! is_serial

!         print*,'leq_bicgs, initial rho(i-1) ', rho(i-1)

         if (rho(i-1) .eq. zero) then
            if(i /= 1)then
! Method fails
! --------------------------------
!               print*, 'leq_bicgs: rho(i-1) == 0 '
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
            P(:) = R(:)
         else
            beta(i-1) = ( rho(i-1)/rho(i-2) )*( alpha(i-1) / omega(i-1) )
               P(:) = R(:) + beta(i-1)*( P(:) - omega(i-1)*V(:) )
         endif ! i.eq.1


! Solve A*Phat(:) = P(:)
! V(:) = A*Phat(:)
! --------------------------------
         if (pc_type.eq.pc_line) then   ! default
            call LEQ_MSOLVE(slo, shi, P, A_m, Phat, sweep_type) ! returns Phat
         else if (pc_type .eq. pc_none) then
            call LEQ_MSOLVE0(slo, shi, P, A_m, Phat, sweep_type) ! returns Phat
         else if (pc_type .eq. pc_diag) then
            call LEQ_MSOLVE1(slo, shi, P, A_m, Phat, sweep_type) ! returns Phat
         end if

         call LEQ_MATVEC(Phat, A_m, V, slo, shi, lo, hi)   ! returns V=A*Phat

         if(is_serial) then
            RtildexV = dot_product( Rtilde, V )
         else
            RtildexV = dot_product_par( Rtilde, V )
         endif ! is_serial

!         print*,'leq_bicgs, initial RtildexV: ', RtildexV

! compute alpha
! --------------------------------
         alpha(i) = rho(i-1) / RtildexV

! compute Svec
! --------------------------------
         Svec(:) = R(:) - alpha(i) * V(:)


! Check norm of Svec(:); if small enough:
! set X(:) = X(:) + alpha(i)*Phat(:) and stop
! --------------------------------
         if(.not.minimize_dotproducts) then
            if(is_serial) then
               Snorm = dot_product( Svec, Svec )
               Snorm = sqrt( Snorm )
            else
               Snorm = sqrt( dot_product_par( Svec, Svec ) )
            endif               ! is_serial
!            print*,'leq_bicgs, initial Snorm: ', real(Snorm)


            if (Snorm <= TOLMIN) then
               Var(:) = Var(:) + alpha(i)*Phat(:)

! Recompute residual norm
! --------------------------------
               if (idebugl >= 1) then
                  call LEQ_MATVEC(Var, A_m, R, slo, shi, lo, hi)   ! returns R=A*Var
!                  Rnorm = sqrt( dot_product_par( Var, Var ) )
!                  print*,'leq_bicgs, initial Vnorm: ', Rnorm

                  do kk = slo(3),shi(3)
                     do jj = slo(2),shi(2)
                        do ii = slo(1),shi(1)
                           IJK = funijk(ii,jj,kk)
                           R(IJK) = B_M(iI,jJ,kK) - R(IJK)
                        enddo
                     enddo
                  enddo

                  if(is_serial) then
                     Rnorm =  dot_product( R, R )
                     Rnorm = sqrt( Rnorm )
                  else
                     Rnorm = sqrt( dot_product_par( R, R ) )
                  endif
!                  print*,'leq_bicgs, initial Rnorm: ', Rnorm
               endif            ! idebugl >= 1
                           isConverged = .TRUE.
                           EXIT
            endif               ! end if (Snorm <= TOLMIN)
         endif                  ! end if (.not.minimize_dotproducts)



! Solve A*Shat(:) = Svec(:)
! Tvec(:) = A*Shat(:)
! --------------------------------
         if (pc_type.eq.pc_line) then   ! default
            call LEQ_MSOLVE(slo, shi, Svec, A_m, Shat, sweep_type)  ! returns Shat
         else if (pc_type .eq. pc_none) then
            call LEQ_MSOLVE0(slo, shi, Svec, A_m, Shat, sweep_type)  ! returns Shat
         else if (pc_type .eq. pc_diag) then
            call LEQ_MSOLVE1(slo, shi, Svec, A_m, Shat, sweep_type)  ! returns Shat
         end if

         call LEQ_MATVEC(Shat, A_m, Tvec, slo, shi, lo, hi )   ! returns Tvec=A*Shat

         if(is_serial) then
            TxS = dot_product( Tvec, Svec )
            TxT = dot_product( Tvec, Tvec )
         else
            if(.not.minimize_dotproducts) then
               TxS = dot_product_par( Tvec, Svec )
               TxT = dot_product_par( Tvec, Tvec )
            else
               TxS_TxT = dot_product_par2(Tvec, Svec, Tvec, Tvec )
               TxS = TxS_TxT(1)
               TxT = TxS_TxT(2)
            endif
         endif

         IF(TxT.eq.Zero) TxT = SMALL_NUMBER

! compute omega
! --------------------------------
         omega(i) = TxS / TxT

! compute new guess for Var
! --------------------------------
            Var(:) = Var(:) +  &
               alpha(i)*Phat(:) + omega(i)*Shat(:)
               R(:) = Svec(:) - omega(i)*Tvec(:)


! --------------------------------
         if(.not.minimize_dotproducts.or.(mod(iter,icheck_bicgs).eq.0)) then
            if(is_serial) then
               Rnorm =  dot_product(R, R )
               Rnorm = sqrt( Rnorm )
            else
               Rnorm = sqrt( dot_product_par(R, R) )
            endif               ! is_serial

            if (idebugl.ge.1) then
               if (myPE.eq.PE_IO) then
                  print*,'iter, Rnorm ', iter, Rnorm, Snorm
                  print*,'alpha(i), omega(i) ', alpha(i), omega(i)
                  print*,'TxS, TxT ', TxS, TxT
                  print*,'RtildexV, rho(i-1) ', RtildexV, rho(i-1)
               endif
            endif

!           call mfix_exit(myPE)

! Check convergence; continue if necessary
! for continuation, it is necessary that omega(i) .ne. 0
            isconverged = (Rnorm <= TOL*Rnorm0)

            if (isconverged) then
               iter_tot(vno) = iter_tot(vno) + iter + 1
               EXIT
            endif
         endif                  ! end if(.not.minimize_dotproducts)

! Advance the iteration count
         iter = iter + 1

      enddo   ! end do i=1,itmax
! end of linear solver loop
! ----------------------------------------------------------------<<<


      if (idebugl >= 1) then
         call LEQ_MATVEC(Var, A_m, R, slo, shi, lo, hi)   ! returns R=A*Var
         do kk = slo(3),shi(3)
            do jj = slo(2),shi(2)
               do ii = slo(1),shi(1)
                  IJK = funijk(ii,jj,kk)
                  R(IJK) = R(IJK) - B_M(iI,jJ,kK)
               enddo
            enddo
         enddo

         if(is_serial) then
            Rnorm = dot_product( R,R)
            Rnorm = sqrt( Rnorm )
         else
            Rnorm = sqrt( dot_product_par( R,R) )
         endif

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
         iter_tot(vno) = iter_tot(vno) + iter
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
