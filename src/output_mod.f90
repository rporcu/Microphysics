!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: output                                                      !
!  Author: M. Syamlal                                 Date: dd-mmm-yy  !
!                                                                      !
!  Purpose: Contain data for output control.                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE output

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use param, only: DIMENSION_USR

! Interval at which user-defined output files are updated.
      real(c_real) :: USR_DT (DIMENSION_USR), USR_TIME(DIMENSION_USR)
! Interval in number of time steps at which LOG file is written
      integer :: NLOG
! Flag to display messages and residuals on the screen
      logical :: FULL_LOG
! The approximated total disk space (in MB)
      real(c_real) :: DISK_TOT = 0.0d0
! One megabite (MB)
      real(c_real), PARAMETER :: ONEMEG = 1048576

! The following have no direct usage in the code. They are generic
! hooks wereby a user can specify some information. It may be useful
! to have the USR_I/J/K variables calculated in a check routine
! the same way that coordinates for BCs, ICs, PSs, are calculated.
!--------------------------------------------------------------------//
! X coordinate of the west face of user output region
      real(c_real) :: USR_X_w (DIMENSION_USR)
! X coordinate of the east face of user output region
      real(c_real) :: USR_X_e (DIMENSION_USR)
! Y coordinate of the south face of user output region
      real(c_real) :: USR_Y_s (DIMENSION_USR)
! Y coordinate of the north face of user output region
      real(c_real) :: USR_Y_n (DIMENSION_USR)
! Z coordinate of the bottom face of user output region
      real(c_real) :: USR_Z_b (DIMENSION_USR)
! Z coordinate of the top face of user output region
      real(c_real) :: USR_Z_t (DIMENSION_USR)
! I index of the west face of user output region
      integer :: USR_I_w (DIMENSION_USR)
! I index of the east face of user output region
      integer :: USR_I_e (DIMENSION_USR)
! J index of the south face of user output region
      integer :: USR_J_s (DIMENSION_USR)
! J index of the north face of user output region
      integer :: USR_J_n (DIMENSION_USR)
! K index of the bottom face of user output region
      integer :: USR_K_b (DIMENSION_USR)
! K index of the top face of user output region
      integer :: USR_K_t (DIMENSION_USR)
! Type of user-defined output: BINARY or ASCII.
      CHARACTER(LEN=6) :: USR_TYPE (DIMENSION_USR)
! Variables to be written in the user-defined output file.
      CHARACTER(LEN=60) :: USR_VAR (DIMENSION_USR)
! Format for writing user-defined (ASCII) output file.
      CHARACTER(LEN=60) :: USR_FORMAT (DIMENSION_USR)
! Extension for the user-defined output file.
      CHARACTER(LEN=16) :: USR_EXT (DIMENSION_USR)


      END MODULE output
