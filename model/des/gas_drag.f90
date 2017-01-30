MODULE GAS_DRAG_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

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
      SUBROUTINE GAS_DRAG_U(slo, shi, &
                            A_M, B_M, f_gds, drag_bm, flag, dx, dy, dz)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values

! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED

! Global Parameters:
!---------------------------------------------------------------------//
      ! IJK of cell to east.
      use functions, only: ieast
      ! Function for averaging to a scalar cell's east face.
      use functions, only: AVG

      IMPLICIT NONE

      integer, intent(in   ) :: slo(3),shi(3)

      ! Septadiagonal matrix A_m
      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)

      ! Vector b_m
      real(c_real), INTENT(INOUT) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), INTENT(IN   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)
      INTEGER, INTENT(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      real(c_real), intent(in   ) :: dx, dy, dz

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K
      real(c_real) :: vol
      vol = dx*dy*dz
!......................................................................!

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN


! Average the interpoalted drag force from the cell corners to the cell face.
      DO K = slo(3),shi(3)
         DO J = slo(2),shi(2)
            DO I = slo(1),shi(1)

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
      SUBROUTINE GAS_DRAG_V(slo, shi, &
                            A_M, B_M, f_gds, drag_bm, flag, dx, dy, dz)


! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values

! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED

! Global Parameters:
!---------------------------------------------------------------------//
      ! IJK of cell to north.
      use functions, only: jnorth

      ! Function for averaging to a scalar cell's north face.
      use functions, only: AVG

      IMPLICIT NONE

      integer, intent(in   ) :: slo(3),shi(3)

! Dummy Arguments:
!---------------------------------------------------------------------//
      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
      real(c_real), INTENT(INOUT) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)
      INTEGER, INTENT(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), intent(in   ) :: dx, dy, dz

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K
      real(c_real) :: vol
      vol = dx*dy*dz
!......................................................................!

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      DO K = slo(3),shi(3)
         DO J = slo(2),shi(2)
            DO I = slo(1),shi(1)
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
      SUBROUTINE GAS_DRAG_W(slo, shi, &
                            A_M, B_M, f_gds, drag_bm, flag, dx, dy, dz)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values

! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED

! Global Parameters:
!---------------------------------------------------------------------//
      ! IJK of cell to top.
      use functions, only: ktop

      ! Function for averaging to a scalar cell's north face.
      use functions, only: AVG

      IMPLICIT NONE

      integer, intent(in   ) :: slo(3),shi(3)

      real(c_real), INTENT(INOUT) :: A_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),-3:3)
      real(c_real), INTENT(INOUT) :: B_m&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)
      INTEGER, INTENT(IN   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
      real(c_real), intent(in   ) :: dx, dy, dz

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K
      real(c_real) :: vol
      vol = dx*dy*dz
!......................................................................!

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      DO K = slo(3),shi(3)
         DO J = slo(2),shi(2)
            DO I = slo(1),shi(1)
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
