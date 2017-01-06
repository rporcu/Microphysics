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

      use geometry, only: imin1, jmin1, kmin1
      use geometry, only: imax1, jmax1, kmax1

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
      double precision  :: yt
! Exact and numerical solutions
      double precision  :: lUg, Ug_MFIX
! Absolute and absolute relative errors
      double precision :: absErr, relErr
! file unit for output data
      integer, parameter :: fUNIT = 2030

! Open file for output.
      OPEN(UNIT=fUNIT, FILE='POST_VEL.dat', &
         POSITION="APPEND",STATUS='OLD')

! Calculate then initial height
      yt = -0.5d0*dy

! generate grid locations for exact solution calculation
      k = kmin1 + (kmax1-kmin1)/2
      i = imin1 + (imax1-imin1)/2
      do j = jmin1, jmax1
! Calculate cell height (assumed uniform spacing)
         yt = yt + dy
! Calculate exact solution
         lUg = Ug(yt)

! Get the MFIX solution
         Ug_MFIX = U_G(i,j,k)

         absErr = abs(lUg - Ug_MFIX)
         relErr = abs(absErr/max(abs(lUg), SMALL_NUMBER))

         write(fUnit, 1200) yt, lUg, Ug_MFIX, absErr, relErr

      end do

      write(fUnit,"(5/)")
      close(fUnit)

 1200 FORMAT(5(3x,es13.6))

      RETURN

      contains


!----------------------------------------------------------------------//
      double precision function Ug(y)

      Use bc, only: delP_x
      Use bc, only: Uw => BC_Uw_g
      Use bc, only: BC_P_g, BC_TYPE
      use fld_const, only: Mu_g0
      Use geometry, only: HEIGHT => YLENGTH
      Use geometry, only: XLENGTH
      Use geometry, only: CYCLIC_X_PD

      implicit none

      double precision, intent(in) :: y
      double precision :: dPdX

      if(CYCLIC_X_PD) then
         dPdX = -delP_x/xLength
      elseif(bc_type(3) == 'P_OUTFLOW' .and. &
         bc_type(4) == 'P_INFLOW') then
         dPdX = (BC_P_g(4) - BC_P_g(3))
      elseif(bc_type(4) == 'P_OUTFLOW' .and. &
         bc_type(3) == 'P_INFLOW') then
         dPdX = (BC_P_g(3) - BC_P_g(4))
      endif

      Ug = Uw(1)*(y/height) + &
            (HALF/Mu_g0)*(dPdX)*(y**2 - height*y)

      end function Ug
      END SUBROUTINE USR3
