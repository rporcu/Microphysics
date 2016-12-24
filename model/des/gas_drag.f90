MODULE GAS_DRAG_MODULE
   CONTAINS

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
      SUBROUTINE GAS_DRAG_U(A_M, B_M, f_gds, drag_am, drag_bm, flag, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values

! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Volume of X-momentum cell
      use geometry, only: VOL

! Global Parameters:
!---------------------------------------------------------------------//
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
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
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
! Error index
      INTEGER, INTENT(INOUT) :: IER


! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN


! Average the interpoalted drag force from the cell corners to the cell face.
      DO K = kstart3, kend3
         DO J = jstart3, jend3
            DO I = istart3, iend3

               IF(flag(i,j,k,2)>= 2000 .and. &
                  flag(i,j,k,2)<=2011) then

                  A_M(I,J,K,0) = A_M(I,J,K,0) - 0.5d0*VOL * &
                     (F_GDS(i,j,k) + F_GDS(ieast(i,j,k),j,k))

                  B_M(I,J,K) = B_M(I,J,K) - 0.5d0* VOL *&
                     (DRAG_BM(i,j,k,1) + DRAG_BM(ieast(i,j,k),j,k,1))

               ENDIF
            ENDDO
         ENDDO
      ENDDO

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
      SUBROUTINE GAS_DRAG_V(A_M, B_M, f_gds, drag_am, drag_bm, flag, IER)


! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values

! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Volume of Y-momentum cell
      use geometry, only: VOL

! Global Parameters:
!---------------------------------------------------------------------//
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
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
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K
! Averaging factor
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      DO K = kstart3, kend3
         DO J = jstart3, jend3
            DO I = istart3, iend3
               IF(flag(i,j,k,3) >= 2000 .and. &
                  flag(i,j,k,3) <= 2011) then
                  A_M(I,J,K,0) = A_M(I,J,K,0) - VOL * 0.5d0*&
                     (F_GDS(I,J,K) + F_GDS(i,jnorth(i,j,k),k))
                  B_M(I,J,K) = B_M(I,J,K) - VOL * 0.5d0*&
                     (DRAG_BM(i,j,k,2)+DRAG_BM(i,jnorth(i,j,k),k,2))
               ENDIF
            ENDDO
         ENDDO
      ENDDO

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
      SUBROUTINE GAS_DRAG_W(A_M, B_M, f_gds, drag_am, drag_bm, flag, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values

! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Volume of Z-momentum cell
      use geometry, only: VOL

! Global Parameters:
!---------------------------------------------------------------------//
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
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
      INTEGER, INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3,4)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      DO K = kstart3, kend3
         DO J = jstart3, jend3
            DO I = istart3, iend3
               IF(flag(i,j,k,4) >= 2000 .and. &
                  flag(i,j,k,4) <= 2011) then
                  A_M(I,J,K,0) = A_M(I,J,K,0) - VOL * 0.5d0*&
                     (F_GDS(I,J,K) + F_GDS(i,j,ktop(i,j,k)))
                  B_M(I,J,K) = B_M(I,J,K) - VOL * 0.5d0*&
                     (DRAG_BM(i,j,k,3) + DRAG_BM(i,j,ktop(i,j,k),3))
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE GAS_DRAG_W
END MODULE GAS_DRAG_MODULE
