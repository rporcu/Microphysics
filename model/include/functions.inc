! This file is included by functions_mod.f

! These function are are defined in this file so it can be included by
! individual subprograms in order to inline the (short) functions
! defined in this file.

! All functions that may benefits from inlining should be defined in this file.

! For code that is not performance-critical, do:

! SUBROUTINE subname
! USE functions
! ...
! END SUBROUTINE subname

! For code that is performance-critical, do:

! SUBROUTINE subname
! ...
! CONTAINS
! INCLUDE 'functions.inc'
! END SUBROUTINE subname

!---------------------------------------------------------------------//
! Functions for generating the LOCAL 3-D array index IJK from the
! 1-D indices I, J, and K.
      INTEGER FUNCTION funijk(li,lj,lk)
      USE compar, only: c0, c1, c2
      implicit none
      INTEGER, INTENT(IN) :: LI, LJ, LK
      funijk = c0 + li + lj*c1 + lk*c2
      END FUNCTION funijk

!---------------------------------------------------------------------//
      ! logical function to identify a fluid cell
      logical function fluid_at(i,j,k)
         use geometry, only: flag
         implicit none
         integer, intent(in) :: i,j,k
         fluid_at = (flag(i,j,k,1).eq.1)
      end function fluid_at


! Flow boundaries
!---------------------------------------------------------------------//
! identify flow at east boundary
      LOGICAL FUNCTION FLOW_AT_E(i,j,k)
      use geometry, only: flag
      implicit none
      INTEGER, INTENT(IN) :: i,j,k
      FLOW_AT_E   = FLAG(i,j,k,2) .GE. 2000 .AND.&
                    FLAG(i,j,k,2) .LE. 2011
      END FUNCTION FLOW_AT_E

! identify specified flow north boundary
      LOGICAL FUNCTION FLOW_AT_N(i,j,k)
      use geometry, only: flag
      implicit none
      INTEGER, INTENT(IN) :: i,j,k
      FLOW_AT_N   = FLAG(i,j,k,3) .GE. 2000 .AND.&
                    FLAG(i,j,k,3) .LE. 2011
      END FUNCTION FLOW_AT_N

! identify specified flow top boundary
      LOGICAL FUNCTION FLOW_AT_T(i,j,k)
      use geometry, only: flag
      implicit none
      INTEGER, INTENT(IN) :: i,j,k
      FLOW_AT_T   = FLAG(i,j,k,4) .GE. 2000 .AND.&
                    FLAG(i,j,k,4) .LE. 2011
      END FUNCTION FLOW_AT_T

! identify specified flow east boundary
      LOGICAL FUNCTION MFLOW_AT_E(i,j,k)
      use geometry, only: flag
      implicit none
      INTEGER i,j,k
      MFLOW_AT_E  = FLAG(i,j,k,2) .EQ. 2020 .OR. &
                    FLAG(i,j,k,2) .EQ. 2021 .OR. &
                    FLAG(i,j,k,2) .EQ. 2031
      END FUNCTION MFLOW_AT_E

! identify specified flow north boundary
      LOGICAL FUNCTION MFLOW_AT_N(i,j,k)
      use geometry, only: flag
      implicit none
      INTEGER i,j,k
      MFLOW_AT_N  = FLAG(i,j,k,3) .EQ. 2020 .OR. &
                    FLAG(i,j,k,3) .EQ. 2021 .OR. &
                    FLAG(i,j,k,3) .EQ. 2031
      END FUNCTION MFLOW_AT_N

! identify specified flow top boundary
      LOGICAL FUNCTION MFLOW_AT_T(i,j,k)
      use geometry, only: flag
      implicit none
      INTEGER i,j,k
      MFLOW_AT_T  = FLAG(i,j,k,4) .EQ. 2020 .OR. &
                    FLAG(i,j,k,4) .EQ. 2021 .OR. &
                    FLAG(i,j,k,4) .EQ. 2031
      END FUNCTION MFLOW_AT_T


! Functions to identify a impermeable and/or semi-permeable surface at
! indicated boundary (specific type of internal surface)
!---------------------------------------------------------------------//
! Logical function to identify IP (impermeable surface) East of the cell
      logical function ip_at_e(i,j,k)
      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k
      ip_at_e = flag(i,j,k,2) .LT. 1000
      end function ip_at_e

! Logical function to identify IP (impermeable surface) North of the cell
      logical function ip_at_n(i,j,k)
      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k
      ip_at_n = flag(i,j,k,3) .LT. 1000
      end function ip_at_n

! Logical function to identify IP (impermeable surface) Top of the cell
      logical function ip_at_t(i,j,k)
      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k
      ip_at_t = flag(i,j,k,4) .LT. 1000
      end function ip_at_t

! Misc
!---------------------------------------------------------------------//


      integer function iminus(i,j,k)
      use compar, only: istart3
      use geometry, only  : imin1,imax1,cyclic_x
      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k

      if (i.eq.imin1 .and.  flag(i,j,k,1) > 1) then
         iminus = i
      else
         iminus = max(istart3, i-1)
      end if

      if (cyclic_x .and. (iminus < imin1) .and. &
         flag(i,j,k,1)<100) iminus = imax1

      end function iminus


      ! ************************************************

      integer function iplus(i,j,k)
      use compar, only: iend3
      use geometry, only  : imin1,imax1,cyclic_x
      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k

      if (i.eq.imax1 .and. flag(i,j,k,1) > 1) then
         iplus = i
      else
         iplus = min(iend3, i+1)
      end if

      if (cyclic_x .and. (iplus > imax1) .and. &
         flag(i,j,k,1)<100) iplus = imin1

      end function iplus

      ! ************************************************

      integer function jminus(i,j,k)
      use compar, only: jstart3
      use geometry, only  : jmin1,jmax1,cyclic_y
      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k

      if (j.eq.jmin1 .and. flag(i,j,k,1) > 1) then
         jminus = j
      else
         jminus = max(jstart3, j-1)
      end if

      if (cyclic_y .and. (jminus < jmin1) .and. &
         flag(i,j,k,1)<100) jminus = jmax1

      end function jminus

      ! ************************************************

      integer function jplus(i,j,k)
      use compar, only: jend3
      use geometry, only  : jmin1,jmax1,cyclic_y
      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k

      if (j.eq.jmax1 .and. flag(i,j,k,1) > 1) then
         jplus = j
      else
         jplus = min(jend3, j+1)
      end if

      if (cyclic_y .and. (jplus > jmax1) .and. &
         flag(i,j,k,1)<100) jplus = jmin1

      end function jplus

      ! ************************************************

      integer function kminus(i,j,k)

      use compar, only: kstart3
      use geometry, only  : kmin1,kmax1,cyclic_z
      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k

      if (k.eq.kmin1 .and. flag(i,j,k,1) > 1) then
         kminus = k
      else
         kminus = max(kstart3, k-1)
      end if

      if (cyclic_z .and. (kminus < kmin1) .and. &
         flag(i,j,k,1)<100) kminus = kmax1

      end function kminus

      ! ************************************************

      integer function kplus(i,j,k)

      use compar, only: kend3
      use geometry, only: kmin1,kmax1,cyclic_z
      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k

      if (k.eq.kmax1 .and. flag(i,j,k,1) > 1) then
         kplus = k
      else
         kplus = min(kend3, k+1)
      end if

      if (cyclic_z .and. (kplus > kmax1) .and. &
         flag(i,j,k,1)<100) kplus = kmin1

      end function kplus

      ! ************************************************

      integer function ieast(i,j,k)
      use geometry, only: flag

      implicit none
      integer, intent(in) :: i,j,k

      if (flag(iplus(i,j,k),j,k,1)>=100) then
         ieast = i
      else
         ieast = iplus(i,j,k)
      end if

      end function ieast

      ! ************************************************

      integer function iwest(i,j,k)

      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k

      if (flag(iminus(i,j,k),j,k,1) >= 100) then
         iwest = i
      else
         iwest = iminus(i,j,k)
      end if

      end function iwest

      ! ************************************************

      integer function jnorth(i,j,k)

      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k

      if (flag(i,jplus(i,j,k),k,1) >= 100) then
         jnorth = j
      else
         jnorth = jplus(i,j,k)
      end if

      end function jnorth

      ! ************************************************

      integer function jsouth(i,j,k)

      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k

      if (flag(i,jminus(i,j,k),k,1) >= 100) then
         jsouth = j
      else
         jsouth = jminus(i,j,k)
      end if

      end function jsouth

      ! ************************************************

      integer function ktop(i,j,k)

      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k
      integer             :: kp

      kp = kplus(i,j,k)

      if (flag(i,j,kp,1)>=100) then
         ktop = k
      else
         ktop = kp
      end if

      end function ktop

      ! ************************************************

      integer function kbot(i,j,k)

      use geometry, only: flag
      implicit none
      integer, intent(in) :: i,j,k
      integer             :: km

      km = kminus(i,j,k)

      if (flag(i,j,km,1)>=100) then
         kbot = k
      else
         kbot = km
      end if

      end function kbot

! ************************************************

      integer function im1(i)

      use compar, only: nodesi, istart3
      use geometry, only  : imin1, imin2, imax1, cyclic_x
      implicit none
      integer, intent(in) :: i

      if(i==imin2) then
         im1 = imin2
      elseif(cyclic_x .and. nodesi == 1 .and. i==imin1 ) then
         im1 = imax1
      else
         im1 = max(istart3, i-1)
      endif

      end function im1

! ************************************************

      integer function ip1(i)

      use compar, only: nodesi, iend3
      use geometry, only  : imin1, imax2, imax1, cyclic_x
      implicit none
      integer, intent(in) :: i

      if(i==imax2) then
         ip1 = imax2
      elseif(cyclic_x .and. nodesi == 1 .and. i==imax1 ) then
         ip1 = imin1
      else
         ip1 = min(iend3, i+1)
      endif

      end function ip1

! ************************************************

      integer function jm1(j)

      use compar, only: nodesj, jstart3
      use geometry, only  : jmin1, jmin2, jmax1, cyclic_y
      implicit none
      integer, intent(in) :: j

      if(j==jmin2) then
         jm1 = jmin2
      elseif(cyclic_y .and. nodesj == 1 .and. j==jmin1 ) then
         jm1 = jmax1
      else
         jm1 = max(jstart3, j-1)
      endif

      end function jm1

! ************************************************

      integer function jp1(j)

      use compar, only: nodesj, jend3
      use geometry, only  : jmin1, jmax2, jmax1, cyclic_y
      implicit none
      integer, intent(in) :: j

      if(j==jmax2) then
         jp1 = jmax2
      elseif(cyclic_y .and. nodesj == 1 .and. j==jmax1 ) then
         jp1 = jmin1
      else
         jp1 = min(jend3, j+1)
      endif

      end function jp1

! ************************************************

      integer function km1(k)

      use compar, only: nodesk, kstart3
      use geometry, only: kmin1, kmin2, kmax1, cyclic_z
      implicit none
      integer, intent(in) :: k

      if(k==kmin2) then
         km1 = kmin2
      elseif(cyclic_z .and. nodesk == 1 .and. k==kmin1 ) then
         km1 = kmax1
      else
         km1 = max(kstart3, k-1)
      endif

      end function km1

! ************************************************

      integer function kp1(k)

      use compar, only: nodesk, kend3
      use geometry, only: kmin1, kmax2, kmax1, cyclic_z
      implicit none
      integer, intent(in) :: k

      if(k==kmax2) then
         kp1 = kmax2
      elseif(cyclic_z .and. nodesk == 1 .and. k==kmax1 ) then
         kp1 = kmin1
      else
         kp1 = min(kend3, k+1)
      endif

      end function kp1

! Misc
!---------------------------------------------------------------------//

! Function that returns the maximum of zero or input
      DOUBLE PRECISION FUNCTION ZMAX(XXX)
      USE param1, only: zero
      implicit none
      DOUBLE PRECISION XXX
      ZMAX       = MAX(XXX, ZERO)
      END FUNCTION ZMAX

! Arithmetic average
      DOUBLE PRECISION FUNCTION AVG(XXXm, XXXp)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) ::XXXp, XXXm
      AVG = 0.5d0 *(XXXm + XXXp)
      END FUNCTION AVG

! Harmonic average
      DOUBLE PRECISION FUNCTION AVG_H(XXXm, XXXp)
      USE param1, only: SMALL_NUMBER
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XXXp, XXXm
      AVG_H = XXXm * XXXp / MAX(SMALL_NUMBER, 0.5d0*(XXXm + XXXp) )
      END FUNCTION AVG_H
