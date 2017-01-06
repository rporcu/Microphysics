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

! Minimum value of solids volume fraction tracked
      real(c_real), PARAMETER :: ZERO_EP_s = 1.0D-8

! Dilute flow threshold. When the volume fraction of a certain phase
! in a cell is smaller than this value the momentum equation for that
! phase is not solved in the cell.
      real(c_real), PARAMETER :: DIL_EP_s = 1.0D-4

! Tolerance used for comparing two numbers for equality in function
! compare(a, b)
      real(c_real), PARAMETER :: TOL_COM = 1.0D-4

! Upper bound for temperatures
      real(c_real), PARAMETER :: TMAX = 4000.D0

! Lower bound for temperatures
      real(c_real), PARAMETER :: TMIN = 250.D0

! Reciprocal of a maximum molecular weight
      real(c_real), PARAMETER :: oMW_MAX = (ONE/500.D0)

! Maximum value of velocities set to avoid divergence problems.
      real(c_real) :: MAX_INLET_VEL

! User definable factor used to scale MAX_INLET_VEL. Default value is 1.
      real(c_real) :: MAX_INLET_VEL_FAC

! Maximum allowed velocity of gas or solids in case no inlet velocities
! (or zero velocities) are defined at inlet (see function
! check_vel_bound)
      real(c_real), PARAMETER :: MAX_ALLOWED_VEL = 500.0D+2

! The following quantities can be specified through the input data
! file, with namelist inputs of the same name.
! ------------------------------------------------------------------->>
! Tolerance in residuals allowed for convergence
      real(c_real) :: TOL_RESID

! Minimum residual for declaring divergence
      real(c_real) :: TOL_DIVERGE

! Factor for normalizing the residual of gas cont. eq.
      real(c_real) :: NORM_g

! -------------------------------------------------------------------<<

! Detect negative Rho_g in physical_prop to reduce DT in iterate
      LOGICAL :: Neg_RHO_G = .FALSE.

      CONTAINS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Purpose: Test if two small values are nearly equal                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      LOGICAL FUNCTION COMPARE (V1, V2)

      use param1, only: small_number, one
      IMPLICIT NONE

! Dummy arguments
! -------------------------------------------------------------------//
! Values to be compared
      real(c_real), INTENT(IN) :: V1, V2

      IF (ABS(V1) <= SMALL_NUMBER) THEN
         IF (ABS(V2) <= SMALL_NUMBER) THEN
            COMPARE = .TRUE.
         ELSE
            COMPARE = .FALSE.
         ENDIF
      ELSE
         IF (ABS(V2/V1 - ONE) <= TOL_COM) THEN
            COMPARE = .TRUE.
         ELSE
            COMPARE = .FALSE.
         ENDIF
      ENDIF
      RETURN
      END FUNCTION COMPARE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Purpose: Test if given variable is smaller than specified tolerance !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      LOGICAL FUNCTION IS_SMALL (V, tol, flag, slo, shi)

      IMPLICIT NONE

      integer     , intent(in   ) :: slo(3),shi(3)

      ! Tolerance value for small
      real(c_real), INTENT(IN) :: tol

      ! Field variable array
      real(c_real), INTENT(IN) :: V&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      integer, intent(in   ) :: flag &
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

! Local variables
!--------------------------------------------------------------------//
      INTEGER :: i, j, k

      IS_SMALL = .FALSE.
      do k = slo(3),shi(3)
        do j = slo(1),shi(2)
          do i = slo(1),shi(1)

             if (1.eq.flag(i,j,k,1)) then
               if (abs(V(i,j,k)) > tol) return
             end if
          end do
        end do
      end do

      IS_SMALL = .TRUE.

      RETURN
      END FUNCTION IS_SMALL

      END MODULE toleranc
