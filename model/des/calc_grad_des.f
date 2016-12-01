
      SUBROUTINE CALC_GRAD_DES(PHI, DEL_PHI)

! Modules
!-----------------------------------------------

      use geometry, only: DO_K
      use geometry, only: oDX, oDY, oDZ
      USE geometry, only: IMIN1, JMIN1, KMIN1
      USE geometry, only: IMAX1, JMAX1, KMAX1
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

      use functions, only: Funijk
      use functions, only: fluid_at
      use functions, only: iplus, iminus, jplus, jminus, kplus, kminus
      use fun_avg, only: AVG

      USE param, only: DIMENSION_3
      USE param1, only: ZERO

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: PHI(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: DEL_PHI(3,DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! general i, j, k indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM
!......................................................................!


        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
         DEL_PHI(:,IJK) = ZERO
         IF(.NOT.fluid_at(i,j,k)) CYCLE

         IMJK = funijk(iminus(i,j,k),j,k)
         IPJK = funijk(iplus(i,j,k),j,k)

         IF((I>IMIN1).AND.(I<IMAX1)) THEN
            DEL_PHI(1,IJK) = oDX*(AVG(PHI(IJK),PHI(IPJK)) -     &
               AVG(PHI(IMJK),PHI(IJK)))
         ELSEIF(I == IMIN1) THEN
            DEL_PHI(1,IJK) = 2.0d0*oDX *                            &
               (AVG(PHI(IJK),PHI(IPJK)) -  PHI(IJK))
         ELSEIF(I == IMAX1) THEN
            DEL_PHI(1,IJK) = 2.0d0*oDX *                            &
               (PHI(IJK) - AVG(PHI(IMJK), PHI(IJK)))
         ELSE
            DEL_PHI(1,IJK) = ZERO
         ENDIF


         IJMK = funijk(i,jminus(i,j,k),k)
         IJPK = funijk(i,jplus(i,j,k),k)

         IF((J>JMIN1) .AND. (J<JMAX1)) THEN
            DEL_PHI(2,IJK) = oDY*(AVG(PHI(IJK),PHI(IJPK)) -     &
               AVG(PHI(IJMK),PHI(IJK)))
         ELSEIF(J == JMIN1) THEN
            DEL_PHI(2,IJK) = 2.0d0*oDY *                            &
               (AVG(PHI(IJK),PHI(IJPK)) - PHI(IJK))
         ELSEIF(J == JMAX1) THEN
            DEL_PHI(2,IJK) = 2.0d0*oDY *                            &
               (PHI(IJK)- AVG(PHI(IJMK),PHI(IJK)))
         ELSE
            DEL_PHI(2,IJK) = ZERO
         ENDIF

         IF(DO_K) THEN

            IJKM = funijk(i,j,kminus(i,j,k))
            IJKP = funijk(i,j,kplus(i,j,k))

            IF((K>KMIN1) .AND. (K<KMAX1)) THEN
               DEL_PHI(3,IJK) = oDZ*(AVG(PHI(IJK),PHI(IJKP)) -  &
                  AVG(PHI(IJKM),PHI(IJK)))
            ELSEIF(K == KMIN1) THEN
               DEL_PHI(3,IJK) = 2.0d0*oDZ *                         &
                  (AVG(PHI(IJK),PHI(IJKP)) - PHI(IJK))
            ELSEIF(K == KMAX1) THEN
               DEL_PHI(3,IJK) = 2.0d0*oDZ *                         &
                  (PHI(IJK) - AVG(PHI(IJKM),PHI(IJK)))
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CALC_GRAD_DES
