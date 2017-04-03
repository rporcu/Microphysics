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

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use param, only: dimension_ps, dim_n_g, dim_m, dim_n_s

! Run-time logical indicating that point sources are present.
      LOGICAL :: POINT_SOURCE

      LOGICAL :: PS_DEFINED(DIMENSION_PS)

! Physical location of point sources.
      real(c_real) :: PS_X_w(DIMENSION_PS)  ! West
      real(c_real) :: PS_X_e(DIMENSION_PS)  ! East
      real(c_real) :: PS_Y_s(DIMENSION_PS)  ! South
      real(c_real) :: PS_Y_n(DIMENSION_PS)  ! North
      real(c_real) :: PS_Z_b(DIMENSION_PS)  ! Bottom
      real(c_real) :: PS_Z_t(DIMENSION_PS)  ! Top

! Cell indices delineating point source region:
      integer :: PS_I_w(DIMENSION_PS)  ! West
      integer :: PS_I_e(DIMENSION_PS)  ! East
      integer :: PS_J_s(DIMENSION_PS)  ! South
      integer :: PS_J_n(DIMENSION_PS)  ! North
      integer :: PS_K_b(DIMENSION_PS)  ! Bottom
      integer :: PS_K_t(DIMENSION_PS)  ! Top

! Gas mass flow rate through the point source:
      real(c_real) PS_MASSFLOW_g (DIMENSION_PS)

! Velocity vector for gas point source: (normalized)
      real(c_real) :: PS_U_g(DIMENSION_PS) ! X-axis
      real(c_real) :: PS_V_g(DIMENSION_PS) ! Y-axis
      real(c_real) :: PS_W_g(DIMENSION_PS) ! Z-axis

! Gas phase velocity magnitude: (calculated)
      real(c_real) :: PS_VEL_MAG_G(DIMENSION_PS)

! Gas phase species mass fractions
      real(c_real) :: PS_X_g(DIMENSION_PS, DIM_N_g)

! Gas phase temperature.
      real(c_real) :: PS_T_g(DIMENSION_PS)
      real(c_real) :: PS_CpxMFLOW_g(DIMENSION_PS)

! Solids mass flow rate through the point source:
      real(c_real) PS_MASSFLOW_s (DIMENSION_PS, DIM_M)

! Velocity vector for solids point sources: (normalized)
      real(c_real) :: PS_U_s(DIMENSION_PS, DIM_M) ! X-axis
      real(c_real) :: PS_V_s(DIMENSION_PS, DIM_M) ! Y-axis
      real(c_real) :: PS_W_s(DIMENSION_PS, DIM_M) ! Z-axis

! Solids phase velocity magnitude: (calculated)
      real(c_real) :: PS_VEL_MAG_S(DIMENSION_PS, DIM_M)

! Solids phase species mass fractions
      real(c_real) :: PS_X_s(DIMENSION_PS, DIM_M, DIM_N_s)

! Solids phase temperature.
      real(c_real) :: PS_T_s(DIMENSION_PS, DIM_M)
      real(c_real) :: PS_CpxMFLOW_s(DIMENSION_PS, DIM_M)

! Total volume of cells comprising point source cells (calculated)
      real(c_real) :: PS_VOLUME(DIMENSION_PS)

! Legacy variable... to be deleated
      integer, DIMENSION(:), ALLOCATABLE :: POINT_SOURCES


      END MODULE ps
