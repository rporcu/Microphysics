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
      subroutine gas_drag_u(lo, hi, slo, shi, alo, ahi, dlo, dhi, &
                            A_m, b_m, f_gds_u, drag_u, vol, domlo, domhi)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values

      ! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED

      integer, intent(in   ) ::  lo(3), hi(3)
      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: dlo(3),dhi(3)
      integer, intent(in   ) :: domlo(3), domhi(3)

      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real), intent(in   ) :: f_gds_u&
         (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
      real(c_real), intent(in   ) :: drag_u&
         (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
      real(c_real), intent(in   ) :: vol

      integer :: I, J, K
!......................................................................!

      ! Skip this routine if the gas/solids are only one-way coupled.
      if (DES_ONEWAY_COUPLED) return

      ! Average the interpolated drag force from the cell corners to the cell face.
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = max(domlo(1)+1,lo(1)), min(domhi(1), hi(1))
               A_m(i,j,k,0) = A_m(i,j,k,0) - vol * f_gds_u(i,j,k)
               b_m(i,j,k  ) = b_m(i,j,k  ) - vol *  drag_u(i,j,k)
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
      subroutine gas_drag_v(lo, hi, slo, shi, alo, ahi, dlo, dhi, &
                            A_m, b_m, f_gds_v, drag_v, vol, domlo, domhi)


      ! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED

      integer, intent(in   ) ::  lo(3), hi(3)
      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: dlo(3),dhi(3)
      integer, intent(in   ) :: domlo(3), domhi(3)

      real(c_real), intent(inout) :: A_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real), intent(in   ) :: f_gds_v&
         (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
      real(c_real), intent(in   ) :: drag_v&
         (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
      real(c_real), intent(in   ) :: vol

      integer :: i, j, k

      ! Skip this routine if the gas/solids are only one-way coupled.
      if (DES_ONEWAY_COUPLED) return

      do k = lo(3), hi(3)
         do j = max(domlo(2)+1,lo(2)), min(domhi(2), hi(2))
            do i = lo(1), hi(1)
               A_m(i,j,k,0) = A_m(i,j,k,0) - vol * f_gds_v(i,j,k)
               b_m(i,j,k  ) = b_m(i,j,k  ) - vol *  drag_v(i,j,k)
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
      subroutine gas_drag_w(lo, hi, slo, shi, alo, ahi, dlo, dhi, &
                            A_m, b_m, f_gds_w, drag_w, vol, domlo, domhi)

      ! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED

      integer, intent(in   ) ::  lo(3), hi(3)
      integer, intent(in   ) :: slo(3),shi(3)
      integer, intent(in   ) :: alo(3),ahi(3)
      integer, intent(in   ) :: dlo(3),dhi(3)
      integer, intent(in   ) :: domlo(3), domhi(3)

      real(c_real), intent(inout) :: a_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),-3:3)
      real(c_real), intent(inout) :: b_m&
         (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
      real(c_real), intent(in   ) :: f_gds_w&
         (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
      real(c_real), intent(in   ) :: drag_w&
         (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
      real(c_real), intent(in   ) :: vol

      integer :: i, j, k

      ! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) return

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               A_m(i,j,k,0) = A_m(i,j,k,0) - vol * f_gds_w(i,j,k)
               b_m(i,j,k  ) = b_m(i,j,k  ) - vol *  drag_w(i,j,k)
            end do
         end do
      end do

      end subroutine gas_drag_w
end module gas_drag_module
