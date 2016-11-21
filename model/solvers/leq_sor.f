!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LEQ_SOR                                                 C
!  Purpose: Solve system of linear system using SOR method             C
!           Successive over-relaxation                                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-AUG-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LEQ_SOR(VNAME, VNO, VAR, A_M, B_M, ITMAX, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE sendrecv
      USE leqsol
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! variable name
      CHARACTER(LEN=*), INTENT(IN) :: Vname
! variable number (not really used here; see calling subroutine)
      INTEGER, INTENT(IN) :: VNO
! variable
      DOUBLE PRECISION, INTENT(INOUT) :: Var(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: &
                        A_m(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: &
                        B_m(DIMENSION_3)
! maximum number of iterations (generally leq_it)
      INTEGER, INTENT(IN) :: ITMAX
! error indicator
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! OVERRELAXATION FACTOR
      DOUBLE PRECISION, PARAMETER :: OMEGA = 1.0  !1.2
      integer :: iidebug
      parameter( iidebug = 0 )
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Variable
      DOUBLE PRECISION :: Var_tmp(DIMENSION_3)
! Indices
      INTEGER :: I,J,K, IJK
      INTEGER :: ITER

      DOUBLE PRECISION oAm
!-----------------------------------------------

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
         IF(.NOT.IS_ON_myPE_owns(i,j,k)) CYCLE

         OAM = ONE/A_M(IJK,0)
         A_M(IJK,0) = ONE
         A_M(IJK,-2) = A_M(IJK,-2)*OAM
         A_M(IJK,-1) = A_M(IJK,-1)*OAM
         A_M(IJK,1) = A_M(IJK,1)*OAM
         A_M(IJK,2) = A_M(IJK,2)*OAM
         A_M(IJK,-3) = A_M(IJK,-3)*OAM
         A_M(IJK,3) = A_M(IJK,3)*OAM
         B_M(IJK) = B_M(IJK)*OAM
      ENDDO
      ENDDO
      ENDDO

      DO ITER = 1, ITMAX
         IF (DO_K) THEN

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
              IJK = FUNIJK(i,j,k)
              IF (.NOT.IS_ON_myPE_owns(i,j,k)) CYCLE
              VAR_tmp(IJK) = VAR(IJK) + OMEGA*(B_M(IJK)-&
                 A_M(IJK,-1)*VAR(funijk(iminus(i,j,k),j,k))-A_M(IJK,1)*VAR(funijk(iplus(i,j,k),j,k))-&
                 A_M(IJK,-2)*VAR(funijk(i,jminus(i,j,k),k))-A_M(IJK,2)*VAR(funijk(i,jplus(i,j,k),k))-&
                 A_M(IJK,-3)*VAR(funijk(i,j,kminus(i,j,k)))-A_M(IJK,3)*VAR(funijk(i,j,kplus(i,j,k)))-&
                 VAR(IJK))
            ENDDO
            ENDDO
            ENDDO
         ELSE

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
             IJK = FUNIJK(i,j,k)
             IF (.NOT.IS_ON_myPE_owns(i,j,k)) CYCLE
             VAR_tmp(IJK) = VAR(IJK) + OMEGA*(B_M(IJK)-&
                  A_M(IJK,-2)*VAR(funijk(i,jminus(i,j,k),k)) - A_M(IJK,2)*VAR(funijk(i,jplus(i,j,k),k))-&
                  A_M(IJK,-1)*VAR(funijk(iminus(i,j,k),j,k)) - A_M(IJK,1)*VAR(funijk(iplus(i,j,k),j,k))-&
                  VAR(IJK))
           ENDDO
           ENDDO
           ENDDO
         ENDIF

      call send_recv(var,2)
      ENDDO

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)
        VAR(IJK) = VAR_tmp(IJK)
      ENDDO
      ENDDO
      ENDDO

      ITER_TOT(VNO) = ITER


      RETURN
      END SUBROUTINE LEQ_SOR
