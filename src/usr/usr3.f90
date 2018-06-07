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
      subroutine usr3(u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
                      p_g, slo, shi, dx, dy, dz)

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      ! NOTE -- this routine may be called with staggered or cell-centered
      !         velocities so be sure to use the components appropriately

      implicit none 

      integer     , intent(in   ) :: ulo(3),uhi(3)
      integer     , intent(in   ) :: vlo(3),vhi(3)
      integer     , intent(in   ) :: wlo(3),whi(3)
      integer     , intent(in   ) :: slo(3),shi(3)

      real(rt), intent(in) :: u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(rt), intent(in) :: v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(rt), intent(in) :: w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
      real(rt), intent(in) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(rt), intent(in) :: dx, dy, dz

      end subroutine usr3
