!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Purpose: This routine is called after the time loop ends and is
!           user-definable.  The user may insert code in this routine
!           or call appropriate user defined subroutines.
!           This routine is not called from an IJK loop, hence
!           all indices are undefined.                                 C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR3(u_g, v_g, w_g, p_g)

      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3

      Use usr

      IMPLICIT NONE

      double precision, intent(in) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision, intent(in) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)


      RETURN
      END SUBROUTINE USR3
