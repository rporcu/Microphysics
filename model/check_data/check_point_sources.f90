MODULE CHECK_POINT_SOURCES_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

! Parameter constants
      use param1, only: UNDEFINED, UNDEFINED_I, IS_DEFINED, IS_UNDEFINED, ZERO

      use get_ps_module, only: get_ps

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_POINT_SOURCES                                     !
!  Author: J. Musser                                  Date: 10-JUN-13  !
!                                                                      !
!  Purpose: Check point source specifications.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_POINT_SOURCES(dx,dy,dz)

! Global Variables:
!---------------------------------------------------------------------//
! Flag: PS geometry was detected.
      use ps, only: PS_DEFINED

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of PS.
      use param, only: DIMENSION_PS

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter for BCs
      INTEGER :: PSV
      real(c_real), intent(in) :: dx,dy,dz
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_POINT_SOURCES")

! Determine which PSs are DEFINED
      CALL CHECK_PS_GEOMETRY

! Loop over all PS arrays.
      DO PSV = 1, DIMENSION_PS

! Verify user input for defined defined PS.
         IF(PS_DEFINED(PSV)) THEN
            call get_ps(PSV,dx,dy,dz)
            CALL CHECK_PS_GAS_PHASE(PSV)
         ELSE
! Verify that no data was defined for unspecified PS.
            CALL CHECK_PS_OVERFLOW(PSV)
         ENDIF
      ENDDO

! Clear the error manager.
      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_POINT_SOURCES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_PS_GEOMETRY                                       !
!  Author: J. Musser                                  Date: 10-JUN-13  !
!                                                                      !
!  Purpose: Check point source specifications.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_PS_GEOMETRY


! Global Variables:
!---------------------------------------------------------------------//
! Flag: PS contains geometric data and/or specified type
      use ps, only: PS_DEFINED, POINT_SOURCE
! User specified: PS geometry
      use ps, only: PS_X_e, PS_X_w, PS_I_e, PS_I_w
      use ps, only: PS_Y_n, PS_Y_s, PS_J_n, PS_J_s
      use ps, only: PS_Z_t, PS_Z_b, PS_K_t, PS_K_b

! Global Parameters:
!---------------------------------------------------------------------//
! The max number of BCs.
      use param, only: DIMENSION_PS

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival


      implicit none


! Local Variables:
!---------------------------------------------------------------------//
! PS loop counter.
      INTEGER :: PSV
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_PS_GEOMETRY")

! Initialize the PS runtime flag.
      POINT_SOURCE = .FALSE.

! Determine which point source indices have values.
      PSV_LP: do PSV = 1, DIMENSION_PS

         IF (IS_DEFINED(PS_X_W(PSV)))   PS_DEFINED(PSV) = .TRUE.
         IF (IS_DEFINED(PS_X_E(PSV)))   PS_DEFINED(PSV) = .TRUE.
         IF (IS_DEFINED(PS_Y_S(PSV)))   PS_DEFINED(PSV) = .TRUE.
         IF (IS_DEFINED(PS_Y_N(PSV)))   PS_DEFINED(PSV) = .TRUE.
         IF (IS_DEFINED(PS_Z_B(PSV)))   PS_DEFINED(PSV) = .TRUE.
         IF (IS_DEFINED(PS_Z_T(PSV)))   PS_DEFINED(PSV) = .TRUE.
         IF (IS_DEFINED(PS_I_W(PSV))) PS_DEFINED(PSV) = .TRUE.
         IF (IS_DEFINED(PS_I_E(PSV))) PS_DEFINED(PSV) = .TRUE.
         IF (IS_DEFINED(PS_J_S(PSV))) PS_DEFINED(PSV) = .TRUE.
         IF (IS_DEFINED(PS_J_N(PSV))) PS_DEFINED(PSV) = .TRUE.
         IF (IS_DEFINED(PS_K_B(PSV))) PS_DEFINED(PSV) = .TRUE.
         IF (IS_DEFINED(PS_K_T(PSV))) PS_DEFINED(PSV) = .TRUE.

! Skip consistency checks if nothing was defined.
         IF (.NOT.PS_DEFINED(PSV)) cycle PSV_LP

! Flag that one or more point sources has been detected.
         POINT_SOURCE = .TRUE.

         IF(IS_UNDEFINED(PS_X_W(PSV)) .AND. IS_UNDEFINED(PS_I_W(PSV))) THEN
            WRITE(ERR_MSG,1101) PSV, 'PS_X_w and PS_I_w '
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(IS_UNDEFINED(PS_X_E(PSV)) .AND. IS_UNDEFINED(PS_I_E(PSV))) THEN
            WRITE(ERR_MSG,1101) PSV, 'PS_X_e and PS_I_e '
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(IS_UNDEFINED(PS_Y_S(PSV)) .AND. IS_UNDEFINED(PS_J_S(PSV))) THEN
            WRITE(ERR_MSG,1101) PSV, 'PS_Y_s and PS_J_s '
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(IS_UNDEFINED(PS_Y_N(PSV)) .AND. IS_UNDEFINED(PS_J_N(PSV))) THEN
            WRITE(ERR_MSG,1101) PSV, 'PS_Y_n and PS_J_n '
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(IS_UNDEFINED(PS_Z_B(PSV)) .AND. IS_UNDEFINED(PS_K_B(PSV))) THEN
            WRITE(ERR_MSG,1101) PSV, 'PS_Z_b and PS_K_b '
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(IS_UNDEFINED(PS_Z_T(PSV)) .AND. IS_UNDEFINED(PS_K_T(PSV))) THEN
            WRITE(ERR_MSG,1101) PSV, 'PS_Z_t and PS_K_t '
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1101 FORMAT('Error 1101: Point source ',I3,' is ill-defined.',/A,     &
         ' are not specified.',/'Please correct the mfix.dat file.')

      ENDDO PSV_LP

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_PS_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_PS_GAS_PHASE                                      !
!  Author: J. Musser                                  Date: 10-JUN-13  !
!                                                                      !
!  Purpose: Check point source specifications.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_PS_GAS_PHASE(PSV)

! Global Variables:
!---------------------------------------------------------------------//
! Gas phase mass flowrate for PS
      use ps, only: PS_MASSFLOW_G
! Gas phase velocity for PS
      use ps, only: PS_U_g, PS_V_g, PS_W_g

! Global Parameters:
!---------------------------------------------------------------------//

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival

      implicit none

! Dummy Arguments:
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: PSV

!......................................................................!

! Initialze the error manager.
      CALL INIT_ERR_MSG("CHECK_PS_GAS_PHASE")


! Check mass flow and velocity
      IF(IS_UNDEFINED(PS_MASSFLOW_G(PSV))) THEN
         IF(IS_DEFINED(PS_U_g(PSV)) .OR. &
            IS_DEFINED(PS_V_g(PSV)) .OR. &
            IS_DEFINED(PS_W_g(PSV))) THEN

            WRITE(ERR_MSG,1100) PSV, trim(iVar('PS_MASSFLOW_G',PSV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error 1100: Invalid specification for point source ',I3,&
         '.',/A,' is undefined but velocity is given.',/'Please ',    &
         'correct the mfix.dat file.')

         ELSE
            PS_MASSFLOW_G(PSV) = ZERO
            PS_U_g(PSV) = ZERO
            PS_V_g(PSV) = ZERO
            PS_W_g(PSV) = ZERO
         ENDIF

      ELSEIF(ABS(PS_MASSFLOW_G(PSV)) < EPSILON(ZERO)) THEN
         IF(ABS(PS_U_g(PSV)) < EPSILON(ZERO) .OR. &
            ABS(PS_V_g(PSV)) < EPSILON(ZERO) .OR. &
            ABS(PS_W_g(PSV)) < EPSILON(ZERO)) THEN

            WRITE(ERR_MSG,1101) PSV, trim(iVar('PS_MASSFLOW_G',PSV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1101 FORMAT('Error 1101: Invalid specification for point source ',I3,&
         '.',/A,' is zero but velocity is given.',/'Please correct ', &
         'the mfix.dat file.')

! Verify a physical mass flow
      ELSEIF(PS_MASSFLOW_G(PSV) < ZERO) THEN
         WRITE(ERR_MSG,1102) PSV, trim(iVar('PS_MASSFLOW_G',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1102 FORMAT('Error 1102: Invalid specifications for point source ',I3,&
         '.',/A,' < 0.0. Point sources can only add mass to a system',/&
         'Please correct the mfix.dat file.')


! Mass flow is specified:
      ELSE

! Velocity does not have to be defined (no momentum source). If the
! components are UNDEFINED, zero them out.
         IF(IS_UNDEFINED(PS_U_g(PSV))) PS_U_g(PSV) = ZERO
         IF(IS_UNDEFINED(PS_V_g(PSV))) PS_V_g(PSV) = ZERO
         IF(IS_UNDEFINED(PS_W_g(PSV))) PS_W_g(PSV) = ZERO

      ENDIF

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_PS_GAS_PHASE



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_PS_OVERFLOW                                       !
!  Author: J. Musser                                  Date: 10-JUN-13  !
!                                                                      !
!  Purpose: Check point source specifications.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_PS_OVERFLOW(PSV)

! Global Variables:
!---------------------------------------------------------------------//
! Gas phase mass flowrate for PS and velocities
      use ps, only: PS_MASSFLOW_G, PS_U_g, PS_V_g, PS_W_g

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival

      implicit none

! Dummy Arguments:
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: PSV
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_PS_OVERFLOW")


      IF(IS_DEFINED(PS_MASSFLOW_G(PSV))) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_MASSFLOW_G',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IS_DEFINED(PS_U_g(PSV))) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_U_g',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IS_DEFINED(PS_V_g(PSV))) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_V_g',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IS_DEFINED(PS_W_g(PSV))) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_W_g',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1010 FORMAT('Error 1010: ',A,' specified in an undefined PS region.')

      END SUBROUTINE CHECK_PS_OVERFLOW
END MODULE CHECK_POINT_SOURCES_MODULE
