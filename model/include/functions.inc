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

! Misc
!---------------------------------------------------------------------//


      integer function iminus(i,j,k)
      use geometry, only  : domlo,domhi,cyclic_x,flag_mod
      implicit none
      integer, intent(in) :: i,j,k

      if (i.eq.domlo(1) .and.  flag_mod(i,j,k,1) > 1) then
         iminus = i
      else
         iminus = max(1, i-1)
      end if

      if (cyclic_x .and. (iminus < domlo(1)) .and. &
         flag_mod(i,j,k,1)<100) iminus = domhi(1)

      end function iminus


      ! ************************************************

      integer function iplus(i,j,k)
      use compar, only: iend3
      use geometry, only  : domlo,domhi,cyclic_x,flag_mod
      implicit none
      integer, intent(in) :: i,j,k

      if (i.eq.domhi(1) .and. flag_mod(i,j,k,1) > 1) then
         iplus = i
      else
         iplus = min(iend3, i+1)
      end if

      if (cyclic_x .and. (iplus > domhi(1)) .and. &
         flag_mod(i,j,k,1)<100) iplus = domlo(1)

      end function iplus

      ! ************************************************

      integer function jminus(i,j,k)
      use geometry, only  : domlo,domhi,cyclic_y,flag_mod
      implicit none
      integer, intent(in) :: i,j,k

      if (j.eq.domlo(2) .and. flag_mod(i,j,k,1) > 1) then
         jminus = j
      else
         jminus = max(1, j-1)
      end if

      if (cyclic_y .and. (jminus < domlo(2)) .and. &
         flag_mod(i,j,k,1)<100) jminus = domhi(2)

      end function jminus

      ! ************************************************

      integer function jplus(i,j,k)
      use compar, only: jend3
      use geometry, only  : domlo,domhi,cyclic_y,flag_mod
      implicit none
      integer, intent(in) :: i,j,k

      if (j.eq.domhi(2) .and. flag_mod(i,j,k,1) > 1) then
         jplus = j
      else
         jplus = min(jend3, j+1)
      end if

      if (cyclic_y .and. (jplus > domhi(2)) .and. &
         flag_mod(i,j,k,1)<100) jplus = domlo(2)

      end function jplus

      ! ************************************************

      integer function kminus(i,j,k)

      use geometry, only  : domlo,domhi,cyclic_z,flag_mod
      implicit none
      integer, intent(in) :: i,j,k

      if (k.eq.domlo(3) .and. flag_mod(i,j,k,1) > 1) then
         kminus = k
      else
         kminus = max(1, k-1)
      end if

      if (cyclic_z .and. (kminus < domlo(3)) .and. &
         flag_mod(i,j,k,1)<100) kminus = domhi(3)

      end function kminus

      ! ************************************************

      integer function kplus(i,j,k)

      use compar, only: kend3
      use geometry, only: domlo,domhi,cyclic_z,flag_mod
      implicit none
      integer, intent(in) :: i,j,k

      if (k.eq.domhi(3) .and. flag_mod(i,j,k,1) > 1) then
         kplus = k
      else
         kplus = min(kend3, k+1)
      end if

      if (cyclic_z .and. (kplus > domhi(3)) .and. &
         flag_mod(i,j,k,1)<100) kplus = domlo(3)

      end function kplus

      ! ************************************************

      integer function ieast(i,j,k)
      use geometry, only: flag_mod

      implicit none
      integer, intent(in) :: i,j,k

      if (flag_mod(iplus(i,j,k),j,k,1)>=100) then
         ieast = i
      else
         ieast = iplus(i,j,k)
      end if

      end function ieast

      ! ************************************************

      integer function iwest(i,j,k)

      use geometry, only: flag_mod
      implicit none
      integer, intent(in) :: i,j,k

      if (flag_mod(iminus(i,j,k),j,k,1) >= 100) then
         iwest = i
      else
         iwest = iminus(i,j,k)
      end if

      end function iwest

      ! ************************************************

      integer function jnorth(i,j,k)

      use geometry, only: flag_mod
      implicit none
      integer, intent(in) :: i,j,k

      if (flag_mod(i,jplus(i,j,k),k,1) >= 100) then
         jnorth = j
      else
         jnorth = jplus(i,j,k)
      end if

      end function jnorth

      ! ************************************************

      integer function jsouth(i,j,k)

      use geometry, only: flag_mod
      implicit none
      integer, intent(in) :: i,j,k

      if (flag_mod(i,jminus(i,j,k),k,1) >= 100) then
         jsouth = j
      else
         jsouth = jminus(i,j,k)
      end if

      end function jsouth

      ! ************************************************

      integer function ktop(i,j,k)

      use geometry, only: flag_mod
      implicit none
      integer, intent(in) :: i,j,k
      integer             :: kp

      kp = kplus(i,j,k)

      if (flag_mod(i,j,kp,1)>=100) then
         ktop = k
      else
         ktop = kp
      end if

      end function ktop

      ! ************************************************

      integer function kbot(i,j,k)

      use geometry, only: flag_mod
      implicit none
      integer, intent(in) :: i,j,k
      integer             :: km

      km = kminus(i,j,k)

      if (flag_mod(i,j,km,1)>=100) then
         kbot = k
      else
         kbot = km
      end if

      end function kbot

! ************************************************

      integer function im1(i)

      use compar, only: nodesi
      use geometry, only  : domlo,domhi,cyclic_x
      implicit none
      integer, intent(in) :: i

      if(i .eq. domlo(1)-1) then
         im1 = domlo(1)-1
      elseif(cyclic_x .and. nodesi == 1 .and. i==domlo(1) ) then
         im1 = domhi(1)
      else
         im1 = max(1, i-1)
      endif

      end function im1

! ************************************************

      integer function ip1(i)

      use compar, only: nodesi, iend3
      use geometry, only  : domlo, domhi, cyclic_x
      implicit none
      integer, intent(in) :: i

      if(i .eq. domhi(1)+1) then
         ip1 = domhi(1)+1
      elseif(cyclic_x .and. nodesi == 1 .and. i==domhi(1) ) then
         ip1 = domlo(1)
      else
         ip1 = min(iend3, i+1)
      endif

      end function ip1

! ************************************************

      integer function jm1(j)

      use compar  , only: nodesj
      use geometry, only  : domlo, domhi, cyclic_y
      implicit none
      integer, intent(in) :: j

      if(j .eq. domlo(2)-1) then
         jm1 = domlo(2)-1
      elseif(cyclic_y .and. nodesj == 1 .and. j==domlo(2) ) then
         jm1 = domhi(2)
      else
         jm1 = max(1, j-1)
      endif

      end function jm1

! ************************************************

      integer function jp1(j)

      use compar, only: nodesj, jend3
      use geometry, only  : domlo, domhi, cyclic_y
      implicit none
      integer, intent(in) :: j

      if(j .eq. domhi(2)+1) then
         jp1 = domhi(2)+1
      elseif(cyclic_y .and. nodesj == 1 .and. j==domhi(2) ) then
         jp1 = domlo(2)
      else
         jp1 = min(jend3, j+1)
      endif

      end function jp1

! ************************************************

      integer function km1(k)

      use compar  , only: nodesk
      use geometry, only: domlo,domhi,cyclic_z
      implicit none
      integer, intent(in) :: k

      if(k .eq. domlo(3)-1) then
         km1 = domlo(3)-1
      elseif(cyclic_z .and. nodesk == 1 .and. k==domlo(3) ) then
         km1 = domhi(3)
      else
         km1 = max(1, k-1)
      endif

      end function km1

! ************************************************

      integer function kp1(k)

      use compar, only: nodesk, kend3
      use geometry, only: domlo,domhi,cyclic_z
      implicit none
      integer, intent(in) :: k

      if(k .eq. domhi(3)+1) then
         kp1 = domhi(3)+1
      elseif(cyclic_z .and. nodesk == 1 .and. k==domhi(3) ) then
         kp1 = domlo(3)
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
