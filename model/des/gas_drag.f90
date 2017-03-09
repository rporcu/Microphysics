module gas_drag_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   implicit none

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: gas_drag_u                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine gas_drag_u(slo, shi, alo, ahi, &
                            A_m, b_m, f_gds, drag_bm, vol)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values

      ! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED

      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: alo(3),ahi(3)

      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)
      real(c_real), intent(in   ) :: vol

      integer :: I, J, K
!......................................................................!

      ! Skip this routine if the gas/solids are only one-way coupled.
      if (DES_ONEWAY_COUPLED) return

      ! Average the interpolated drag force from the cell corners to the cell face.
      DO K = alo(3),ahi(3)
         DO J = alo(2),ahi(2)
            DO I = alo(1),ahi(1)

                  A_m(I,J,K,0) = A_m(I,J,K,0) - 0.5d0*VOL * &
                     (f_gds(i,j,k) + f_gds(i+1,j,k))

                  b_m(I,J,K) = b_m(I,J,K) - 0.5d0* VOL *&
                     (drag_bm(i,j,k,1) + drag_bm(i+1,j,k,1))

!              end if
            end do
         end do
      end do

      end subroutine gas_drag_u

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: gas_drag_v                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine gas_drag_v(slo, shi, alo, ahi, &
                            A_m, b_m, f_gds, drag_bm, vol)


      ! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED

      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: alo(3),ahi(3)

      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(inout) :: B_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)
      real(c_real), intent(in   ) :: vol

      integer :: I, J, K

      ! Skip this routine if the gas/solids are only one-way coupled.
      if (DES_ONEWAY_COUPLED) RETURN

      DO K = alo(3),ahi(3)
         DO J = alo(2),ahi(2)
            DO I = alo(1),ahi(1)
                  A_m(I,J,K,0) = A_m(I,J,K,0) - VOL * 0.5d0*&
                     (f_gds(I,J,K) + f_gds(i,j+1,k))
                  b_m(I,J,K) = b_m(I,J,K) - VOL * 0.5d0*&
                     (drag_bm(i,j,k,2)+drag_bm(i,j+1,k,2))
            end do
         end do
      end do

      end subroutine gas_drag_v

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: gas_drag_w                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine gas_drag_w(slo, shi, alo, ahi, &
                            A_m, b_m, f_gds, drag_bm, vol)

      ! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED

      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: alo(3),ahi(3)

      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(inout) :: B_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real), intent(in   ) :: f_gds&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: drag_bm&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),3)
      real(c_real), intent(in   ) :: vol

      integer :: I, J, K

      ! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      DO K = alo(3),ahi(3)
         DO J = alo(2),ahi(2)
            DO I = alo(1),ahi(1)
               A_m(I,J,K,0) = A_m(I,J,K,0) - VOL * 0.5d0*&
                  (f_gds(I,J,K) + f_gds(i,j,k+1))
               b_m(I,J,K) = b_m(I,J,K) - VOL * 0.5d0*&
                  (drag_bm(i,j,k,3) + drag_bm(i,j,k+1,3))
            end do
         end do
      end do

      end subroutine gas_drag_w
end module gas_drag_module
