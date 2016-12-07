!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ZERO_NORM_VEL                                           C
!  Purpose: Set the velocity component normal to a wall to zero        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ZERO_NORM_VEL(u_g,v_g,w_g,flag)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1   , only: zero
      USE compar   , only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE geometry , only: imax2, jmax2, kmax2
      USE functions, only: iminus, jminus, kminus
      USE functions, only: ip_at_e, ip_at_n, ip_at_t

      IMPLICIT NONE

      double precision, intent(inout) ::  u_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) ::  v_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) ::  w_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      integer, intent(in) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,0:4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      integer :: i,j,k

      do k = kstart3, kend3
        do j = jstart3, jend3
          do i = istart3, iend3


            if (flag(i,j,k,1)<100) then
               if (ip_at_e(i,j,k)) u_g(i,j,k) = zero
               if (ip_at_n(i,j,k)) v_g(i,j,k) = zero
               if (ip_at_t(i,j,k)) w_g(i,j,k) = zero
            else

               U_G(I,J,K) = ZERO
               V_G(I,J,K) = ZERO
               W_G(I,J,K) = ZERO

               if(flag(i,j,k,1) /= 106 .and. flag(i,j,k,1) /= 107) then
                  if(i < imax2) u_g(iminus(i,j,k),j,k) = zero
                  if(j < jmax2) v_g(i,jminus(i,j,k),k) = zero
                  if(k < kmax2) w_g(i,j,kminus(i,j,k)) = zero
               endif
            endif

          end do
        end do
      end do

      END SUBROUTINE ZERO_NORM_VEL
