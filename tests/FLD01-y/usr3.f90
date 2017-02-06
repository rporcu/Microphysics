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
      SUBROUTINE USR3(slo, shi, u_g, v_g, w_g, p_g, dx, dy, dz)

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use geometry, only: jmax,kmax,domlo,domhi

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      real(c_real), intent(in) :: u_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in) :: v_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in) :: w_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in) :: dx, dy, dz

      ! looping indices
      integer :: i, j, k

      ! Calculated height of cell
      double precision  :: yt, zt

      ! Exact and Numerical solutions
      double precision :: lvg, vg_MFIX
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
      zt = -0.5d0*dz

! Calculate the U velocity solution at center of domain
      i = domlo(1) + (domhi(1)-domlo(1))/2
      j = domlo(2) + (domhi(2)-domlo(2))/2
      do k = domlo(3), domhi(3)
         zt = zt + dz

         ! Calculate the exact solution
         lvg = vg(zt)
         vg_MFIX = v_g(i,j,k)

         absErr = abs(lvg - vg_MFIX)
         relErr = abs(absErr/max(abs(lvg), 1.0d-15))

         L1 = L1 + absErr
         L2 = L2 + absErr*absErr
         LI = max(LI, absErr)

         write(fUnit,1100) zt, lvg, vg_MFIX, absErr, relErr
      end do

      close(fUnit)


      open(unit=fUnit, file='POST_PG.dat', &
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
         lPg = Pg(yt)
         Pg_MFIX = P_G(i,j,k)

         write(6,*)yt,lPg,pg_mfix

         absErr = abs(lPg - Pg_MFIX)
         relErr = abs(absErr/max(abs(lPg), 1.0d-15))

         L1 = L1 + absErr
         L2 = L2 + absErr*absErr
         LI = max(LI, absErr)

         write(fUnit,1100) yt, lPg, Pg_MFIX, absErr, relErr
      end do
      close(fUnit)

      RETURN
 1100 Format(5(3x,es13.6))
 1200 Format(3x,I3,3(3x,es13.6))

      contains

!----------------------------------------------------------------------!
! Function: Calculate the exact solution for pressure.                 !
!----------------------------------------------------------------------!
      double precision function Pg(y)

      use ic, only: Pg0 => IC_P_g
      use bc, only: delP_y
      use geometry, only: yLength

      double precision, intent(in) :: y
      double precision :: dPdy

      dPdy = -delP_y/yLength

      Pg = Pg0(1) + dPdY*(y - yLength)

      end function Pg

!----------------------------------------------------------------------!
! Function: Calculate the exact solution for pressure.                 !
!----------------------------------------------------------------------!
      double precision function vg(z)

      use bc, only: delP_y
      use geometry, only: yLength
      use geometry, only: zLength

      use fld_const, only: Mu_g0

      double precision, intent(in) :: z
      double precision :: dPdY

      dPdy = -delP_y/yLength

      vg = (1.0d0/(2.0d0*mu_g0))*dPdY*(z**2 -zLength*z)

      end function vg

      end subroutine usr3
