!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Purpose: This routine is called after the time loop ends and is     C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.              C
!           This routine is not called from an IJK loop, hence         C
!           all indices are undefined.                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE USR3(u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
      p_g, slo, shi, dx, dy, dz)

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use geometry, only: domlo,domhi

      implicit none

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

! Calculated height of cell
      double precision  :: xt, yt
! Exact and Numerical solutions
      double precision :: lUg, Ug_MFIX
      double precision :: lPg, Pg_MFIX
! Absolute and absolute relative errors
      double precision :: absErr, relErr
! File unit for output data
      integer, parameter :: fUnit= 2030
! Norm Errors
      double precision :: L1, L2, LI

! Open file for output
      open(unit=fUnit, file='POST_UG.dat', &
         position='append', status='old')

      L1 = 0.0d0
      L2 = 0.0d0
      LI = 0.0d0

! Calculate the initial height
      yt = -0.5d0*dy

! Calculate the U velocity solution at center of domain
      k = domlo(3) + (domhi(3)-domlo(3))/2
      i = domlo(1) + (domhi(1)-domlo(1))/2
      do j = domlo(2), domhi(2)
         yt = yt + dy
! Calculate the exact solution
         lUg = Ug(yt)
         Ug_MFIX = U_G(i,j,k)

         absErr = abs(lUg - Ug_MFIX)
         relErr = abs(absErr/max(abs(lUg), 1.0d-15))

         L1 = L1 + absErr
         L2 = L2 + absErr*absErr
         LI = max(LI, absErr)

         write(fUnit,1100) yt, lUg, Ug_MFIX, absErr, relErr
      end do

      close(fUnit)

      open(unit=fUnit, file='POST_PG.dat', &
         position='append', status='old')

      L1 = 0.0d0
      L2 = 0.0d0
      LI = 0.0d0

! Calculate the initial height
      xt = -0.5d0*dx

! Calculate the U velocity solution at center of domain
      k = domlo(3) + (domhi(3)-domlo(3))/2
      j = domlo(2) + (domhi(2)-domlo(2))/2

      do i = domlo(1), domhi(1)
         xt = xt + dx
! Calculate the exact solution
         lPg = Pg(xt)
         Pg_MFIX = P_G(i,j,k)

         write(6,*)xt,lPg,pg_mfix

         absErr = abs(lPg - Pg_MFIX)
         relErr = abs(absErr/max(abs(lPg), 1.0d-15))

         L1 = L1 + absErr
         L2 = L2 + absErr*absErr
         LI = max(LI, absErr)

         write(fUnit,1100) xt, lPg, Pg_MFIX, absErr, relErr
      end do
      close(fUnit)


      RETURN
 1100 Format(5(3x,es13.6))
 1200 Format(3x,I3,3(3x,es13.6))

      contains

!----------------------------------------------------------------------!
! Function: Calculate the exact solution for pressure.                 !
!----------------------------------------------------------------------!
      double precision function Pg(x)

      use ic, only: Pg0 => IC_P_g
      use bc, only: delP_x
      use geometry, only: xLength

      double precision, intent(in) :: x
      double precision :: dPdX

      dPdx = -delP_x/xLength

      Pg = Pg0(1) + dPdX*(x - xLength)

      end function Pg

!----------------------------------------------------------------------!
! Function: Calculate the exact solution for pressure.                 !
!----------------------------------------------------------------------!
      double precision function Ug(y)

      use bc, only: delP_x
      use geometry, only: xLength
      use geometry, only: yLength

      use fld_const, only: Mu_g0

      double precision, intent(in) :: y
      double precision :: dPdX

      dPdx = -delP_x/xLength

      Ug = (1.0d0/(2.0d0*mu_g0))*dPdX*(y**2 -yLength*y)

      end function Ug

      end subroutine usr3
