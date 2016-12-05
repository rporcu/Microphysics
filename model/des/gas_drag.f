!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_DRAG_W                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_DRAG_U(A_M, B_M, f_gds, drag_am, drag_bm, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use particle_filter, only: DES_INTERP_SCHEME_ENUM, DES_INTERP_GARG
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Volume of X-momentum cell
      use geometry, only: VOL

! Global Parameters:
!---------------------------------------------------------------------//
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
! Flag: Fluid exists at indexed cell
      use functions, only: flow_at_e
! IJK of cell to east.
      use functions, only: ieast
! Function for averaging to a scalar cell's east face.
      use functions, only: AVG
! Domain index bounds.
      use compar, only: ISTART2, JSTART2, KSTART2
      use compar, only: IEND2, JEND2, KEND2

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: drag_am&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
! Error index
      INTEGER, INTENT(INOUT) :: IER


! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION :: tmp_A, tmp_B
! Averaging factor
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN


! Average the interpoalted drag force from the cell corners to the cell face.
      IF(DES_INTERP_SCHEME_ENUM == DES_INTERP_GARG)THEN

         AVG_FACTOR = 0.25d0

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

            IF(flow_at_e(i,j,k)) then

               IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
               IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
               IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

               tmp_A = -AVG_FACTOR*(DRAG_AM(i,j,k) + DRAG_AM(I, J-1, K))
               tmp_B = -AVG_FACTOR*(DRAG_BM(i,j,k,1) + DRAG_BM(I, J-1, K,1))

               tmp_A = tmp_A - AVG_FACTOR*                             &
                  (DRAG_AM(I, J, K-1) + DRAG_AM(I, J-1, K-1))
               tmp_B = tmp_B - AVG_FACTOR*                             &
                  (DRAG_BM(I, J, K-1,1) + DRAG_BM(I, J-1, K-1,1))

               A_M(I,J,K,0) = A_M(I,J,K,0) + tmp_A*VOL
               B_M(I,J,K) = B_M(I,J,K) + tmp_B*VOL
            endif
         ENDDO
         ENDDO
         ENDDO

      ELSE

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

            IF(flow_at_e(i,j,k)) THEN

               tmp_A = AVG(F_GDS(i,j,k), F_GDS(ieast(i,j,k),j,k))
               tmp_B = AVG(DRAG_BM(i,j,k,1), DRAG_BM(ieast(i,j,k),j,k,1))

               A_M(I,J,K,0) = A_M(I,J,K,0) - VOL * tmp_A
               B_M(I,J,K) = B_M(I,J,K) - VOL * tmp_B
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDIF



      END SUBROUTINE GAS_DRAG_U

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_DRAG_V                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_DRAG_V(A_M, B_M, f_gds, drag_am, drag_bm, IER)


! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use particle_filter, only: DES_INTERP_SCHEME_ENUM, DES_INTERP_GARG
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Volume of Y-momentum cell
      use geometry, only: VOL

! Global Parameters:
!---------------------------------------------------------------------//
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
! Flag: Fluid exists at indexed cell
      use functions, only: flow_at_n
! IJK of cell to north.
      use functions, only: jnorth
! Function for averaging to a scalar cell's north face.
      use functions, only: AVG
! Domain index bounds.
      use compar, only: ISTART2, JSTART2, KSTART2
      use compar, only: IEND2, JEND2, KEND2


      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: drag_am&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION tmp_A, tmp_B
! Averaging factor
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      IF(DES_INTERP_SCHEME_ENUM == DES_INTERP_GARG)THEN

         AVG_FACTOR = 0.25d0

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

            IF(flow_at_n(i,j,k)) then

               IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
               IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
               IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

               tmp_A = -AVG_FACTOR*(DRAG_AM(i,j,k) + DRAG_AM(I-1,J,K))
               tmp_B = -AVG_FACTOR*(DRAG_BM(i,j,k,2) + DRAG_BM(I-1,J,K,2))

               tmp_A = tmp_A - AVG_FACTOR*                             &
                  (DRAG_AM(I,J,K-1) + DRAG_AM(I-1,J,K-1))
               tmp_B = tmp_B - AVG_FACTOR*                             &
                  (DRAG_BM(I,J,K-1,2) + DRAG_BM(I-1,J,K-1,2))

               A_M(I,J,K,0) = A_M(I,J,K,0) + tmp_A*VOL
               B_M(I,J,K) = B_M(I,J,K) + tmp_B*VOL
            endif
         ENDDO
         ENDDO
         ENDDO

      ELSE

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

            IF(flow_at_n(i,j,k)) THEN

               tmp_A = AVG(F_GDS(I,J,K), F_GDS(i,jnorth(i,j,k),k))
               tmp_B = AVG(DRAG_BM(i,j,k,2), DRAG_BM(i,jnorth(i,j,k),k,2))

               A_M(I,J,K,0) = A_M(I,J,K,0) - VOL * tmp_A
               B_M(I,J,K) = B_M(I,J,K) - VOL * tmp_B
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE GAS_DRAG_V

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_DRAG_W                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_DRAG_W(A_M, B_M, f_gds, drag_am, drag_bm, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use particle_filter, only: DES_INTERP_SCHEME_ENUM, DES_INTERP_GARG
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Volume of Z-momentum cell
      use geometry, only: VOL

! Global Parameters:
!---------------------------------------------------------------------//
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
! Flag: Fluid exists at indexed cell
      use functions, only: flow_at_t
! IJK of cell to top.
      use functions, only: ktop
! Function for averaging to a scalar cell's north face.
      use functions, only: AVG
! Domain index bounds.
      use compar, only: ISTART2, JSTART2, KSTART2
      use compar, only: IEND2, JEND2, KEND2

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(INOUT) :: A_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, -3:3)
      DOUBLE PRECISION, INTENT(INOUT) :: B_m&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: f_gds&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: drag_am&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: drag_bm&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,3)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION tmp_A, tmp_B
! Averaging factor
! (=0.25 in 3D and =0.5 in 2D)
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      IF(DES_INTERP_SCHEME_ENUM == DES_INTERP_GARG)THEN

         AVG_FACTOR = 0.25d0

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

            IF(flow_at_t(i,j,k)) then

               IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
               IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
               IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

               tmp_A = -AVG_FACTOR*(DRAG_AM(i,j,k) + DRAG_AM(I-1,J,K) +        &
                  DRAG_AM(I,J-1,K) + DRAG_AM(I-1,J-1,K))

               tmp_B = -AVG_FACTOR*(DRAG_BM(i,j,k,3) + DRAG_BM(I-1,J,K,3) +    &
                  DRAG_BM(I,J-1,K,3) + DRAG_BM(I-1,J-1,K,3))

               A_M(I,J,K,0) = A_M(I,J,K,0) + tmp_A*VOL
               B_M(I,J,K) = B_M(I,J,K) + tmp_B*VOL
            endif
         ENDDO
         ENDDO
         ENDDO

      ELSE

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

            IF(flow_at_t(i,j,k)) THEN
               tmp_A = AVG(F_GDS(I,J,K), F_GDS(i,j,ktop(i,j,k)))
               tmp_B = AVG(DRAG_BM(i,j,k,3), DRAG_BM(i,j,ktop(i,j,k),3))

               A_M(I,J,K,0) = A_M(I,J,K,0) - VOL * tmp_A
               B_M(I,J,K) = B_M(I,J,K) - VOL * tmp_B
            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ENDIF

      RETURN
      END SUBROUTINE GAS_DRAG_W
