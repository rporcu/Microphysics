!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Module name: ps_mod.f                                               C
!                                                                      C
!  Purpose: Common block containing point source data.                 C
!                                                                      C
!  Author: J. Musser                                  Date: 10-Jun-13  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE ps

      use amrex_fort_module, only : rt => amrex_real
      use iso_c_binding , only: c_int

      use param, only: dim_ps, dim_n_g, dim_m, dim_n_s

! Run-time logical indicating that point sources are present.
      logical :: POINT_SOURCE

! Physical location of point sources.
      real(rt) :: PS_X_w(DIM_PS)  ! West
      real(rt) :: PS_X_e(DIM_PS)  ! East
      real(rt) :: PS_Y_s(DIM_PS)  ! South
      real(rt) :: PS_Y_n(DIM_PS)  ! North
      real(rt) :: PS_Z_b(DIM_PS)  ! Bottom
      real(rt) :: PS_Z_t(DIM_PS)  ! Top

! Gas mass flow rate through the point source:
      real(rt) :: PS_MASSFLOW_g (DIM_PS)

! Velocity vector for gas point source: (normalized)
      real(rt) :: PS_U_g(DIM_PS) ! X-axis
      real(rt) :: PS_V_g(DIM_PS) ! Y-axis
      real(rt) :: PS_W_g(DIM_PS) ! Z-axis

! Gas phase velocity magnitude: (calculated)
      real(rt) :: PS_VEL_MAG_G(DIM_PS)

! Gas phase species mass fractions
      real(rt) :: PS_X_g(DIM_PS, DIM_N_g)

! Gas phase temperature.
      real(rt) :: PS_T_g(DIM_PS)
      real(rt) :: PS_CpxMFLOW_g(DIM_PS)

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: ps_defined                                               !
!                                                                      !
! Purpose: Return if a PS region has been defined based on coordinates !
! defined in the input deck.                                           !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   logical function ps_defined(icv)

      use param, only: is_defined

      integer, intent(in) :: icv

      ps_defined = is_defined(ps_x_w(icv)) .or. is_defined(ps_x_e(icv)) .or. &
                   is_defined(ps_y_s(icv)) .or. is_defined(ps_y_n(icv)) .or. &
                   is_defined(ps_z_b(icv)) .or. is_defined(ps_z_t(icv))

   end function ps_defined

end module ps
