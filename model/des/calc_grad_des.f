      SUBROUTINE CALC_GRAD_DES(PHI, DEL_PHI)

! Modules
!-----------------------------------------------

      use geometry, only: oDX, oDY, oDZ
      USE geometry, only: IMIN1, JMIN1, KMIN1
      USE geometry, only: IMAX1, JMAX1, KMAX1
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

      use functions, only: fluid_at
      use functions, only: iplus, iminus, jplus, jminus, kplus, kminus
      use functions, only: AVG

      USE param1, only: ZERO

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: PHI&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(OUT) :: DEL_PHI&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)

! Local variables
!---------------------------------------------------------------------//
! general i, j, k indices
      INTEGER :: I, J, K
!......................................................................!


        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         DEL_PHI(I,J,K,:) = ZERO
         IF(.NOT.fluid_at(i,j,k)) CYCLE

         IF((I>IMIN1).AND.(I<IMAX1)) THEN
            DEL_PHI(I,J,K,1) = oDX*(AVG(PHI(I,J,K),PHI(iplus(i,j,k),j,k)) -     &
               AVG(PHI(iminus(i,j,k),j,k),PHI(I,J,K)))
         ELSEIF(I == IMIN1) THEN
            DEL_PHI(I,J,K,1) = 2.0d0*oDX *                            &
               (AVG(PHI(I,J,K),PHI(iplus(i,j,k),j,k)) -  PHI(I,J,K))
         ELSEIF(I == IMAX1) THEN
            DEL_PHI(I,J,K,1) = 2.0d0*oDX *                            &
               (PHI(I,J,K) - AVG(PHI(iminus(i,j,k),j,k), PHI(I,J,K)))
         ELSE
            DEL_PHI(I,J,K,1) = ZERO
         ENDIF


         IF((J>JMIN1) .AND. (J<JMAX1)) THEN
            DEL_PHI(I,J,K,2) = oDY*(AVG(PHI(I,J,K),PHI(i,jplus(i,j,k),k)) -     &
               AVG(PHI(i,jminus(i,j,k),k),PHI(I,J,K)))
         ELSEIF(J == JMIN1) THEN
            DEL_PHI(I,J,K,2) = 2.0d0*oDY *                            &
               (AVG(PHI(I,J,K),PHI(i,jplus(i,j,k),k)) - PHI(I,J,K))
         ELSEIF(J == JMAX1) THEN
            DEL_PHI(I,J,K,2) = 2.0d0*oDY *                            &
               (PHI(I,J,K)- AVG(PHI(i,jminus(i,j,k),k),PHI(I,J,K)))
         ELSE
            DEL_PHI(I,J,K,2) = ZERO
         ENDIF


         IF((K>KMIN1) .AND. (K<KMAX1)) THEN
            DEL_PHI(I,J,K,3) = oDZ*(AVG(PHI(I,J,K),PHI(i,j,kplus(i,j,k))) -  &
               AVG(PHI(i,j,kminus(i,j,k)),PHI(I,J,K)))
         ELSEIF(K == KMIN1) THEN
            DEL_PHI(I,J,K,3) = 2.0d0*oDZ *                         &
               (AVG(PHI(I,J,K),PHI(i,j,kplus(i,j,k))) - PHI(I,J,K))
         ELSEIF(K == KMAX1) THEN
            DEL_PHI(I,J,K,3) = 2.0d0*oDZ *                         &
               (PHI(I,J,K) - AVG(PHI(i,j,kminus(i,j,k)),PHI(I,J,K)))
         ELSE
            DEL_PHI(I,J,K,3) = ZERO
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CALC_GRAD_DES
