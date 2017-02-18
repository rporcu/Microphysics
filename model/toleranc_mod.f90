!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: tolerance                                              C
!  Purpose: Specify all tolerance parameters                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUL-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE toleranc

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use param1, only: one

      ! Tolerance used for comparing two numbers for equality in function
      ! compare(a, b)
      real(c_real), parameter :: tol_com = 1.0D-4

      ! Minimum value of solids volume fraction tracked
      real(c_real), parameter :: zero_ep_s = 1.0d-8

      ! Maximum value of velocities set to avoid divergence problems.
      real(c_real) :: max_inlet_vel

      ! User definable factor used to scale max_inlet_vemax_inlet_vel. Default value is 1.
      real(c_real) :: max_inlet_vel_fac

      ! Maximum allowed velocity of gas or solids in case no inlet velocities
      ! (or zero velocities) are defined at inlet (see function check_vel_bound)
      real(c_real), parameter :: max_allowed_vel = 500.0d+2

      ! The following quantities can be specified through the input data
      ! file, with namelist inputs of the same name.
      ! ------------------------------------------------------------------->>

      ! Tolerance in residuals allowed for convergence
      real(c_real) :: tol_resid

      ! Minimum residual for declaring divergence
      real(c_real) :: tol_diverge

      ! Factor for normalizing the residual of gas cont. eq.
      real(c_real) :: norm_g

      ! -------------------------------------------------------------------<<

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Purpose: Test if two small values are nearly equal                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      logical function compare (v1, v2)

      use param1, only: small_number, one

      implicit none

      ! Values to be compared
      real(c_real), intent(in) :: v1, v2

      if (abs(v1) <= small_number) then
         if (abs(v2) <= small_number) then
            compare = .true.
         else
            compare = .false.
         end if
      else
         if (abs(v2/v1 - one) <= tol_com) then
            compare = .true.
         else
            compare = .false.
         end if

      end if
      end function compare

      end module toleranc
