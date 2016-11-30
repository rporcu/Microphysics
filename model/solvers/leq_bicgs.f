!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine LEQ_BICGS                                                C
!  Purpose: Solve system of linear system using BICGS method           C
!           Biconjugate gradients stabilized                           C
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
      SUBROUTINE LEQ_BICGS(VNAME, VNO, VAR, A_M, B_m, cmethod, &
                           TOL, PC, ITMAX, IER)

         use compar, only: istart3, iend3
         use compar, only: jstart3, jend3
         use compar, only: kstart3, kend3
         use compar, only: mype
         use exit_mod, only: mfix_exit
         use funits, only: unit_log, dmp_log
         use leqsol, only: leq_matvec, leq_msolve, leq_msolve0, leq_msolve1
         use param, only: DIM_EQS
         use param, only: dimension_3
         use param1, only: zero

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number (not really used here; see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: Var
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
! Note: this setting only seems to matter when leq_pc='line'
      CHARACTER(LEN=*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
! preconditioner (leq_pc)
!     options = 'line' (default), 'diag', 'none'
      CHARACTER(LEN=4), INTENT(IN) ::  PC
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! error indicator
      INTEGER, INTENT(INOUT) :: IER
!-------------------------------------------------
! Local Variables
!-------------------------------------------------

      if(PC.eq.'LINE') then   ! default
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,  &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE, IER )
      elseif(PC.eq.'DIAG') then
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,   &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE1, IER )
      elseif(PC.eq.'NONE') then
         call LEQ_BICGS0( Vname, Vno, Var, A_m, B_m,   &
            cmethod, TOL, ITMAX, LEQ_MATVEC, LEQ_MSOLVE0, IER )
      else
         IF(DMP_LOG)WRITE (UNIT_LOG,*) &
           'preconditioner option not found - check mfix.dat and readme'
         call mfix_exit(myPE)
      endif

      return
      END SUBROUTINE LEQ_BICGS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_BICGS0                                              C
!  Purpose: Compute residual of linear system                          C
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

      SUBROUTINE LEQ_BICGS0(VNAME, VNO, VAR, A_M, B_m, cmethod, &
                            TOL, ITMAX, MATVEC, MSOLVE, IER )

         use compar, only: istart2, iend2
         use compar, only: istart3, iend3
         use compar, only: jstart2, jend2
         use compar, only: jstart3, jend3
         use compar, only: kstart2, kend2
         use compar, only: kstart3, kend3
         use compar, only: mype, pe_io, numpes
         use exit_mod, only: mfix_exit
         use functions, only: funijk
         use funits, only: unit_log, dmp_log
         use leqsol, only: is_serial, icheck_bicgs, minimize_dotproducts
         use leqsol, only: leq_matvec, leq_msolve, leq_msolve0, leq_msolve1, dot_product_par, dot_product_par2, iter_tot
         use param, only: DIM_EQS
         use param, only: dimension_3
         use param1, only: zero, one, small_number

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments/procedure
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number (not really used here-see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable
      DOUBLE PRECISION, DIMENSION(DIMENSION_3), INTENT(INOUT) :: Var
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! Sweep direction of leq solver (leq_sweep)
!     e.g., options = 'isis', 'rsrs' (default), 'asas'
      CHARACTER(LEN=*), INTENT(IN) :: CMETHOD
! convergence tolerance (generally leq_tol)
      DOUBLE PRECISION, INTENT(IN) :: TOL
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
      DOUBLE PRECISION, PARAMETER :: ratiotol = 0.2
      LOGICAL, PARAMETER :: do_unit_scaling = .true.
!-----------------------------------------------
! Local variables
!-----------------------------------------------

      DOUBLE PRECISION, DIMENSION(:), allocatable :: R,Rtilde, P,Phat, Svec, Shat, Tvec,V

      DOUBLE PRECISION, DIMENSION(0:ITMAX+1) :: &
                        alpha, beta, omega, rho
      DOUBLE PRECISION :: TxS, TxT, RtildexV, &
                          aijmax, oam
      DOUBLE PRECISION :: Rnorm, Rnorm0, Snorm, TOLMIN, pnorm
      LOGICAL :: isconverged
      INTEGER :: i, j, k, ijk, ii, jj, kk
      INTEGER :: iter
      DOUBLE PRECISION, DIMENSION(2) :: TxS_TxT
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

         do k = kstart2,kend2
            do i = istart2,iend2
               do j = jstart2,jend2
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
      call MATVEC(Vname, Var, A_M, R)   ! returns R=A*Var

      do k = kstart3,kend3
         do i = istart3,iend3
            do j = jstart3,jend3
               IJK = funijk(i,j,k)
               R(IJK) = B_M(I,J,K) - R(IJK)
            enddo
         enddo
      enddo

      ! call send_recv(R,nlayers_bicgs)

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
         if(myPE.eq.0) print*,'leq_bicgs, initial: ', Vname,' resid ', Rnorm0
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

!         print*,'leq_bicgs, initial: ', Vname,' rho(i-1) ', rho(i-1)

         if (rho(i-1) .eq. zero) then
            if(i /= 1)then
! Method fails
! --------------------------------
!               print*, 'leq_bicgs,',Vname,': rho(i-1) == 0 '
               ier = -2
            else
! Method converged.  residual is already zero
! --------------------------------
               ier = 0
            endif
            ! call send_recv(var,2)
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
         call MSOLVE(Vname, P, A_m, Phat, CMETHOD) ! returns Phat
         call MATVEC(Vname, Phat, A_m, V)   ! returns V=A*Phat

         if(is_serial) then
            RtildexV = dot_product( Rtilde, V )
         else
            RtildexV = dot_product_par( Rtilde, V )
         endif ! is_serial

!         print*,'leq_bicgs, initial: ', Vname,' RtildexV ', RtildexV

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
!            print*,'leq_bicgs, initial: ', Vname,' Snorm ', real(Snorm)


            if (Snorm <= TOLMIN) then
               Var(:) = Var(:) + alpha(i)*Phat(:)

! Recompute residual norm
! --------------------------------
               if (idebugl >= 1) then
                  call MATVEC(Vname, Var, A_m, R)   ! returns R=A*Var
!                  Rnorm = sqrt( dot_product_par( Var, Var ) )
!                  print*,'leq_bicgs, initial: ', Vname,' Vnorm ', Rnorm

                  do kk = kstart3,kend3
                     do jj = jstart3,jend3
                        do ii = istart3,iend3
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
!                  print*,'leq_bicgs, initial: ', Vname,' Rnorm ', Rnorm
               endif            ! idebugl >= 1
                           isConverged = .TRUE.
                           EXIT
            endif               ! end if (Snorm <= TOLMIN)
         endif                  ! end if (.not.minimize_dotproducts)



! Solve A*Shat(:) = Svec(:)
! Tvec(:) = A*Shat(:)
! --------------------------------
         call MSOLVE( Vname, Svec, A_m, Shat, CMETHOD)  ! returns Shat
         call MATVEC( Vname, Shat, A_m, Tvec )   ! returns Tvec=A*Shat

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
         call MATVEC(Vname, Var, A_m, R)   ! returns R=A*Var
         do kk = kstart3,kend3
            do jj = jstart3,jend3
               do ii = istart3,iend3
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

         if(myPE.eq.0)  print*,'leq_bicgs ratio : ', Vname,' ',iter,  &
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

      ! call send_recv(var,2)

      deallocate(R)
      deallocate(Rtilde)
      deallocate(P)
      deallocate(Phat)
      deallocate(Svec)
      deallocate(Shat)
      deallocate(Tvec)
      deallocate(V)

      return
      end subroutine LEQ_BICGS0
