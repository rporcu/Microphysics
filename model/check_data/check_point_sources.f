!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_POINT_SOURCES                                     !
!  Author: J. Musser                                  Date: 10-JUN-13  !
!                                                                      !
!  Purpose: Check point source specifications.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_POINT_SOURCES

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
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter for BCs
      INTEGER :: PSV
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_POINT_SOURCES")

! Determine which PSs are DEFINED
      CALL CHECK_PS_GEOMETRY

! Loop over all PS arrays.
      DO PSV = 1, DIMENSION_PS

! Verify user input for defined defined PS.
         IF(PS_DEFINED(PSV)) THEN
            CALL GET_PS(PSV)
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
! User specified: System geometry
      use geometry, only: NO_I, XLENGTH
      use geometry, only: NO_J, YLENGTH
      use geometry, only: NO_K, ZLENGTH

! Global Parameters:
!---------------------------------------------------------------------//
! The max number of BCs.
      use param, only: DIMENSION_PS
! Parameter constants
      use param1, only: ZERO, UNDEFINED, UNDEFINED_I

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


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

         IF (PS_X_W(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE.
         IF (PS_X_E(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE.
         IF (PS_Y_S(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE.
         IF (PS_Y_N(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE.
         IF (PS_Z_B(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE.
         IF (PS_Z_T(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE.
         IF (PS_I_W(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE.
         IF (PS_I_E(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE.
         IF (PS_J_S(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE.
         IF (PS_J_N(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE.
         IF (PS_K_B(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE.
         IF (PS_K_T(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE.

! Skip consistency checks if nothing was defined.
         IF (.NOT.PS_DEFINED(PSV)) cycle PSV_LP

! Flag that one or more point sources has been detected.
         POINT_SOURCE = .TRUE.

         IF(PS_X_W(PSV)==UNDEFINED .AND. PS_I_W(PSV)==UNDEFINED_I) THEN
            IF(NO_I) THEN
               PS_X_W(PSV) = ZERO
            ELSE
               WRITE(ERR_MSG,1101) PSV, 'PS_X_w and PS_I_w '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
         IF(PS_X_E(PSV)==UNDEFINED .AND. PS_I_E(PSV)==UNDEFINED_I) THEN
            IF(NO_I) THEN
               PS_X_E(PSV) = XLENGTH
            ELSE
               WRITE(ERR_MSG,1101) PSV, 'PS_X_e and PS_I_e '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
         IF(PS_Y_S(PSV)==UNDEFINED .AND. PS_J_S(PSV)==UNDEFINED_I) THEN
            IF(NO_J) THEN
               PS_Y_S(PSV) = ZERO
            ELSE
               WRITE(ERR_MSG,1101) PSV, 'PS_Y_s and PS_J_s '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
         IF(PS_Y_N(PSV)==UNDEFINED .AND. PS_J_N(PSV)==UNDEFINED_I) THEN
            IF(NO_J) THEN
               PS_Y_N(PSV) = YLENGTH
            ELSE
               WRITE(ERR_MSG,1101) PSV, 'PS_Y_n and PS_J_n '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
         IF(PS_Z_B(PSV)==UNDEFINED .AND. PS_K_B(PSV)==UNDEFINED_I) THEN
            IF(NO_K) THEN
               PS_Z_B(PSV) = ZERO
            ELSE
               WRITE(ERR_MSG,1101) PSV, 'PS_Z_b and PS_K_b '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
         IF(PS_Z_T(PSV)==UNDEFINED .AND. PS_K_T(PSV)==UNDEFINED_I) THEN
            IF(NO_K) THEN
               PS_Z_T(PSV) = ZLENGTH
            ELSE
               WRITE(ERR_MSG,1101) PSV, 'PS_Z_t and PS_K_t '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
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
! Gas phase tempture and species mass fractions
      use ps, only: PS_T_g, PS_X_g

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: ZERO, ONE, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      use toleranc

      implicit none

! Dummy Arguments:
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: PSV

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: N
! Sum of solids mass fractions.
      DOUBLE PRECISION :: SUM
!......................................................................!

! Initialze the error manager.
      CALL INIT_ERR_MSG("CHECK_PS_GAS_PHASE")


! Check mass flow and velocity
      IF(PS_MASSFLOW_G(PSV) == UNDEFINED) THEN
         IF(PS_U_g(PSV) /= UNDEFINED .OR. &
            PS_V_g(PSV) /= UNDEFINED .OR. &
            PS_W_g(PSV) /= UNDEFINED) THEN

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

      ELSEIF(PS_MASSFLOW_G(PSV) == ZERO) THEN
         IF(PS_U_g(PSV) /= ZERO .OR. &
            PS_V_g(PSV) /= ZERO .OR. &
            PS_W_g(PSV) /= ZERO) THEN

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
         IF(PS_U_g(PSV) == UNDEFINED) PS_U_g(PSV) = ZERO
         IF(PS_V_g(PSV) == UNDEFINED) PS_V_g(PSV) = ZERO
         IF(PS_W_g(PSV) == UNDEFINED) PS_W_g(PSV) = ZERO

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

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum input array sizes.
      use param, only: DIM_M, DIM_N_g, DIM_N_s
! Parameter constants
      use param1, only: UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Dummy Arguments:
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: PSV

! Local Variables:
!---------------------------------------------------------------------//
! Loop counters
      INTEGER :: M, N
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_PS_OVERFLOW")


      IF(PS_MASSFLOW_G(PSV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_MASSFLOW_G',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(PS_U_g(PSV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_U_g',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(PS_V_g(PSV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_V_g',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(PS_W_g(PSV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_W_g',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1010 FORMAT('Error 1010: ',A,' specified in an undefined PS region.')

      END SUBROUTINE CHECK_PS_OVERFLOW
