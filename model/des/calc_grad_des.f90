MODULE CALC_GRAD_DES_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
      SUBROUTINE CALC_GRAD_DES(PHI, DEL_PHI, flag, dx, dy, dz)

! Modules
!-----------------------------------------------

      USE geometry, only: IMIN1, JMIN1, KMIN1
      USE geometry, only: IMAX1, JMAX1, KMAX1
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

      use functions, only: iplus, iminus, jplus, jminus, kplus, kminus
      use functions, only: AVG

      USE param1, only: ZERO

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
      real(c_real), INTENT(IN) :: PHI&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_real), INTENT(OUT) :: DEL_PHI&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
      INTEGER, DIMENSION(:,:,:,:), INTENT(IN) :: FLAG
      real(c_real), intent(in   ) :: dx, dy, dz

! Local variables
!---------------------------------------------------------------------//
! general i, j, k indices
      INTEGER :: I, J, K
      real(c_real) :: odx, ody, odz
!......................................................................!

        odx = 1.d0 / dx
        ody = 1.d0 / dy
        odz = 1.d0 / dz

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         DEL_PHI(I,J,K,:) = ZERO
         IF(.NOT.1.eq.flag(i,j,k,1)) CYCLE

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
   END MODULE CALC_GRAD_DES_MODULE
