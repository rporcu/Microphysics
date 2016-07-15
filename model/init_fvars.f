!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: INIT_FVARS                                             !
!  Purpose: Initialize all field variables.                            !
!                                                                      !
!  Author: M. Syamlal                                 Date: 23-JAN-94  !
!  Reviewer: J.Musser                                 Date:  8-Oct-13  !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_FVARS

! Global Variables:
!---------------------------------------------------------------------//
! Gas phase volume farction
      USE fldvar, only: EP_G
! Pressures
      USE fldvar, only: P_G    ! Gas
! Densities
      USE fldvar, only: RO_G  ! Gas
! Bulk densities: RO*EP
      USE fldvar, only: ROP_G ! Gas
! Gas velocity components:
      USE fldvar, only: U_G  ! x-axis
      USE fldvar, only: V_G  ! y-axis
      USE fldvar, only: W_G  ! z-axis

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: UNDEFINED
      USE param1, only: ZERO


      IMPLICIT NONE

! Passed Variables:
!---------------------------------------------------------------------//
! NONE

! Local Variables:
!---------------------------------------------------------------------//
! NONE

      IF(allocated(EP_G)) EP_G = UNDEFINED

      IF(allocated(P_G)) P_G  = UNDEFINED

      IF(allocated(RO_G)) RO_G = UNDEFINED
      IF(allocated(ROP_G)) ROP_G = UNDEFINED

      IF(allocated(U_G)) U_G = UNDEFINED
      IF(allocated(V_G)) V_G = UNDEFINED
      IF(allocated(W_G)) W_G = UNDEFINED


      RETURN
      END SUBROUTINE INIT_FVARS
