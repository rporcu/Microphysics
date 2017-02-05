!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Author: Avinash Vaidheeswaran                      Date: July, 2016 C
!  Reviewer: J.Musser                                                  C
!                                                                      C
!  Purpose: Calculates the exact solution to the Couette flow case and C
!  compares with the MFIX solution.                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
SUBROUTINE USR3(slo, shi, u_g, v_g, w_g, p_g, dx, dy, dz)

      use geometry, only: domlo, domhi

      use param1, only: small_number, half
      use bl_fort_module, only : c_real

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      double precision, intent(in) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      double precision, intent(in) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      double precision, intent(in) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      double precision, intent(in) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in) :: dx, dy, dz

      ! looping indices
      integer :: i, j, k

      ! Calculated height of cell center
      double precision  :: xt

      ! Exact and numerical solutions
      double precision  :: lwg, wg_MFIX

      ! Absolute and absolute relative errors
      double precision :: absErr, relErr

      ! file unit for output data
      integer, parameter :: fUNIT = 2030

      ! Open file for output.
      OPEN(UNIT=fUNIT, FILE='POST_VEL.dat', &
         POSITION="APPEND",STATUS='OLD')

      ! Calculate then initial height
      xt = -0.5d0*dx

      ! generate grid locations for exact solution calculation
      j = domlo(2) + (domhi(2)-domlo(2))/2
      k = domlo(3) + (domhi(3)-domlo(3))/2
      do i = domlo(1), domhi(1)

         ! Calculate cell height (assumed uniform spacing)
         xt = xt + dx

         ! Calculate exact solution
         lwg = wg(xt)

         ! Get the MFIX solution
         wg_MFIX = w_g(i,j,k)

         absErr = abs(lwg - wg_MFIX)
         relErr = abs(absErr/max(abs(lwg), SMALL_NUMBER))

         write(fUnit, 1200) xt, lwg, wg_MFIX, absErr, relErr

      end do

      write(fUnit,"(5/)")
      close(fUnit)

 1200 FORMAT(5(3x,es13.6))

      RETURN

      contains


!----------------------------------------------------------------------//
      double precision function wg(x)

      Use bc, only: delP_z
      Use bc, only: ww => BC_ww_g
      Use bc, only: BC_P_g, BC_TYPE
      use fld_const, only: Mu_g0
      Use geometry, only: HEIGHT => XLENGTH
      Use geometry, only: ZLENGTH
      Use geometry, only: CYCLIC_Z_PD

      implicit none

      double precision, intent(in) :: x
      double precision :: dPdz

      if(CYCLIC_Z_PD) then
         dPdz = -delP_z/zLength
      elseif(bc_type(3) == 'P_OUTFLOW' .and. &
             bc_type(4) == 'P_INFLOW') then
         dPdz = (BC_P_g(4) - BC_P_g(3))
      elseif(bc_type(4) == 'P_OUTFLOW' .and. &
             bc_type(3) == 'P_INFLOW') then
         dPdz = (BC_P_g(3) - BC_P_g(4))
      endif

      wg = ww(1)*(x/height) + &
            (HALF/Mu_g0)*(dPdz)*(x**2 - height*x)

      end function wg

END SUBROUTINE USR3
