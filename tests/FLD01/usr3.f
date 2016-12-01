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
      SUBROUTINE USR3(u_g, v_g, w_g, p_g)

      use geometry, only: dx, dy

      use geometry, only: imin1, imax1, imax
      use geometry, only: jmin1, jmax1, jmax
      use geometry, only: kmin1, kmax1

      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3

      IMPLICIT NONE

      double precision, intent(in) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

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
      k = kmin1 + (kmax1-kmin1)/2
      i = imin1 + (imax1-imin1)/2
      do j = jmin1, jmax1
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

! Store Norms for U velocity
      open(unit=fUnit, file='POST_UG_NORMS.dat', &
         position='append', status='old')

      L1 = L1/dble(jmax1-jmin1+1)
      L2 = L2/dble(jmax1-jmin1+1)

      write(fUnit,1200) jmax, L1, L2, LI
      close(fUnit)

      open(unit=fUnit, file='POST_PG.dat', &
         position='append', status='old')

      L1 = 0.0d0
      L2 = 0.0d0
      LI = 0.0d0

! Calculate the initial height
      xt = -0.5d0*dx

! Calculate the U velocity solution at center of domain
      k = kmin1 + (kmax1-kmin1)/2
      j = jmin1 + (jmax1-jmin1)/2
      do i = imin1, imax1
         xt = xt + dx
! Calculate the exact solution
         lPg = Pg(xt)
         Pg_MFIX = P_G(i,j,k)

         absErr = abs(lPg - Pg_MFIX)
         relErr = abs(absErr/max(abs(lPg), 1.0d-15))

         L1 = L1 + absErr
         L2 = L2 + absErr*absErr
         LI = max(LI, absErr)

         write(fUnit,1100) xt, lPg, Pg_MFIX, absErr, relErr
      end do
      close(fUnit)

      open(unit=fUnit, file='POST_PG_NORMS.dat', &
         position='append', status='old')

      L1 = L1/dble(imax1-imin1+1)
      L2 = L2/dble(imax1-imin1+1)

      write(fUnit,1200) imax, L1, L2, LI
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
