module set_max2_module

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: set_max2                                               !
!  Purpose: calculate domain bounds.                                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine set_max2

      use geometry, only: imax, jmax, kmax
      use geometry, only: domlo, domhi

      implicit none

      domlo(1) = 0 + 2
      domlo(2) = 0 + 2
      domlo(3) = 0 + 2

      domhi(1) = imax-1 + 2
      domhi(2) = jmax-1 + 2
      domhi(3) = kmax-1 + 2

      end subroutine set_max2

end module set_max2_module
