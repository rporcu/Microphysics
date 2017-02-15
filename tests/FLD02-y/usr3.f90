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
   SUBROUTINE USR3(u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
      p_g, slo, shi, dx, dy, dz)

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use geometry, only: domlo, domhi

      use param1, only: small_number, half

      IMPLICIT NONE

      integer, intent(in   ) :: ulo(3),uhi(3)
      integer, intent(in   ) :: vlo(3),vhi(3)
      integer, intent(in   ) :: wlo(3),whi(3)
      integer, intent(in   ) :: slo(3),shi(3)

      real(c_real), intent(inout) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(c_real), intent(inout) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in   ) :: dx, dy, dz


      ! looping indices
      integer :: i, j, k

      ! Calculated height of cell center
      double precision  :: zt

      ! Exact and numerical solutions
      double precision  :: lvg, vg_MFIX

      ! Absolute and absolute relative errors
      double precision :: absErr, relErr

      ! file unit for output data
      integer, parameter :: fUNIT = 2030

      ! Open file for output.
      OPEN(UNIT=fUNIT, FILE='POST_VEL.dat', &
         POSITION="APPEND",STATUS='OLD')

      ! Calculate then initial height
      zt = -0.5d0*dz

      ! generate grid locations for exact solution calculation
      i = domlo(1) + (domhi(1)-domlo(1))/2
      j = domlo(2) + (domhi(2)-domlo(2))/2
      do k = domlo(3), domhi(3)

         ! Calculate cell height (assumed uniform spacing)
         zt = zt + dz

         ! Calculate exact solution
         lvg = vg(zt)

         ! Get the MFIX solution
         vg_MFIX = v_g(i,j,k)

         absErr = abs(lvg - vg_MFIX)
         relErr = abs(absErr/max(abs(lvg), SMALL_NUMBER))

         write(fUnit, 1200) zt, lvg, vg_MFIX, absErr, relErr

      end do

      write(fUnit,"(5/)")
      close(fUnit)

 1200 FORMAT(5(3x,es13.6))

      RETURN

      contains


!----------------------------------------------------------------------//
      double precision function vg(z)

      Use bc, only: delP_y
      Use bc, only: vw => BC_vw_g
      Use bc, only: BC_P_g, BC_TYPE
      use fld_const, only: Mu_g0
      Use geometry, only: HEIGHT => ZLENGTH
      Use geometry, only: YLENGTH
      Use geometry, only: CYCLIC_Y_PD

      implicit none

      double precision, intent(in) :: z
      double precision :: dPdy

      if(CYCLIC_Y_PD) then
         dPdy = -delP_y/yLength
      elseif(bc_type(3) == 'P_OUTFLOW' .and. &
             bc_type(4) == 'P_INFLOW') then
         dPdy = (BC_P_g(4) - BC_P_g(3))
      elseif(bc_type(4) == 'P_OUTFLOW' .and. &
             bc_type(3) == 'P_INFLOW') then
         dPdy = (BC_P_g(3) - BC_P_g(4))
      endif

      vg = vw(1)*(z/height) + &
            (HALF/Mu_g0)*(dPdy)*(z**2 - height*z)

      end function vg

END SUBROUTINE USR3
