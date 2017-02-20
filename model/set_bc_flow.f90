MODULE set_bc_flow_module

! Global Variables:
!---------------------------------------------------------------------//
! Total number of solids.
      use constant, only: MMAX
! Flag: BC dimensions or Type is specified
      use bc, only: BC_DEFINED
! Use specified BC type
      use bc, only: BC_TYPE
! User specifed BC solids bulk density
      use bc, only: BC_ROP_s
! Solids volume fraction at BC
      use bc, only: BC_EP_s
      use bc, only: BC_EP_g

      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: bc_u_s, bc_v_s, bc_w_s
      use bc, only: bc_plane

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: ZERO, ONE, UNDEFINED, IS_UNDEFINED, IS_DEFINED, EQUAL
! Maximum number of BCs
      use param, only: DIMENSION_BC
! Maximum number of disperse phases
      use param, only: DIM_M

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_FLOW                                             !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Check boundary condition specifications                    !
!     - convert physical locations to i, j, k's                        !
!     - compute area of boundary surfaces (GET_BC_AREA)                !
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
!     - check specification of physical quantities                     !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_FLOW

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar

      use flow_to_vel_new_module, only: flow_to_vel_new

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop counters
      INTEGER :: BCV, I
! Total number of solids phases (continuum + discrete)
      INTEGER :: MMAX_TOT
! Flag to skip checks on indexed solid phase.
      LOGICAL :: SKIP(1:DIM_M)
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_BC_FLOW")

! Total number of solids.
      MMAX_TOT = MMAX

! Loop over each defined BC and check the user data.
      DO BCV = 1, DIMENSION_BC

         IF(.NOT.BC_DEFINED(BCV)) CYCLE

! Determine which solids phases are present.
         SKIP = .FALSE.
         DO I = 1, DIM_M
            IF ((EQUAL(BC_ROP_S(BCV,I), UNDEFINED).OR.EQUAL(BC_ROP_S(BCV,I), ZERO)) &
               .AND.(EQUAL(BC_EP_S(BCV,I), UNDEFINED).OR.EQUAL(BC_EP_S(BCV,I), ZERO))) THEN
               SKIP = .TRUE.
            ENDIF
         END DO

         IF(MMAX_TOT == 1 .AND. .NOT.EQUAL(BC_EP_g(BCV), ONE)) SKIP(1) = .FALSE.

         SELECT CASE (TRIM(BC_TYPE(BCV)))

         CASE ('MASS_INFLOW')
            CALL FLOW_TO_VEL_NEW(.TRUE., MMAX_TOT, SKIP, BCV)
            CALL CHECK_BC_VEL_INFLOW(MMAX_TOT, SKIP, BCV)

         CASE ('MASS_OUTFLOW')
            CALL FLOW_TO_VEL_NEW(.TRUE., MMAX_TOT, SKIP, BCV)
            CALL CHECK_BC_VEL_OUTFLOW(MMAX_TOT, SKIP, BCV)
         END SELECT
      ENDDO

! Cleanup and exit.
      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE SET_BC_FLOW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_VEL_INFLOW                                      !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
! Comments:                                                            !
!     The velocities at the inflow face are fixed and the momentum     !
!     equations are not solved in the inflow cells. Since the flow is  !
!     into the domain all other scalars that are used need to be       !
!     specified (e.g., mass fractions, void fraction, etc.,)           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_VEL_INFLOW(M_TOT, SKIP, BCV)

      USE param, only: DIM_M

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE


      INTEGER, intent(in) :: BCV
      INTEGER, intent(in) :: M_TOT

      LOGICAL, intent(in) :: SKIP(DIM_M)

! loop/variable indices
      INTEGER :: M

      CALL INIT_ERR_MSG("CHECK_BC_VEL_INFLOW")


! Check that gas phase velocities are defined.
      IF(IS_UNDEFINED(BC_U_G(BCV))) THEN
         WRITE(ERR_MSG,1000) trim(iVar('BC_U_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF (IS_UNDEFINED(BC_V_G(BCV))) THEN
         WRITE(ERR_MSG,1000) trim(iVar('BC_V_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(IS_UNDEFINED(BC_W_G(BCV))) THEN
         WRITE(ERR_MSG,1000) trim(iVar('BC_W_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Check that solids phase velocities are defined.
      DO M = 1, M_TOT
         IF(IS_UNDEFINED(BC_U_S(BCV,M))) THEN
            IF(SKIP(M)) THEN
               BC_U_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_U_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(IS_UNDEFINED(BC_V_S(BCV,M))) THEN
            IF(SKIP(M)) THEN
               BC_V_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_V_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(IS_UNDEFINED(BC_W_S(BCV,M))) THEN
            IF(SKIP(M)) THEN
               BC_W_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_W_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ENDDO

! Check that gas phase velocities are consistent.
      SELECT CASE (bc_plane(BCV))

      CASE ('W')
         IF(BC_U_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_U_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_U_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M)), '<'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('E')
         IF(BC_U_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_U_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_U_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M)), '>'
              CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('S')
         IF(BC_V_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_V_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_V_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M)), '<'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('N')
         IF(BC_V_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_V_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_V_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M)), '>'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('B')
         IF(BC_W_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_W_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_W_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M)), '<'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('T')
         IF(BC_W_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_W_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_W_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M)), '>'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      END SELECT

 1300 FORMAT('Error 1300: Invalid flow direction. ',A,' should be ',   &
         A,' zero. ',/'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_VEL_INFLOW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_VEL_OUTFLOW                                     !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
! Comments:                                                            !
!     The velocities at the outflow face are fixed and the momentum    !
!     equations are not solved in the outflow cells. Since the flow    !
!     is out of the domain none of the other scalars should need to    !
!     be specified (e.g., mass fractions, void fraction, etc.,).       !
!     Such values will become defined according to their adjacent      !
!     fluid cell                                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_VEL_OUTFLOW(M_TOT, SKIP, BCV)

      use bc, only: dim_m
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

! loop/variable indices
      INTEGER, intent(in) :: BCV
      INTEGER, intent(in) :: M_TOT
      LOGICAL, intent(in) :: SKIP(DIM_M)

! Loop variable
      INTEGER :: M

      CALL INIT_ERR_MSG("CHECK_BC_VEL_OUTFLOW")

! Check that gas phase velocities are defined.
      IF(IS_UNDEFINED(BC_U_G(BCV))) THEN
         WRITE(ERR_MSG,1000) trim(iVar('BC_U_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF (IS_UNDEFINED(BC_V_G(BCV))) THEN
         WRITE(ERR_MSG,1000) trim(iVar('BC_V_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(IS_UNDEFINED(BC_W_G(BCV))) THEN
         WRITE(ERR_MSG,1000) trim(iVar('BC_W_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Check that solids phase velocities are defined.
      DO M = 1, M_TOT
         IF(IS_UNDEFINED(BC_U_S(BCV,M))) THEN
            IF(SKIP(M)) THEN
               BC_U_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_U_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(IS_UNDEFINED(BC_V_S(BCV,M))) THEN
            IF(SKIP(M)) THEN
               BC_V_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_V_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(IS_UNDEFINED(BC_W_S(BCV,M))) THEN
            IF(SKIP(M)) THEN
               BC_W_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_W_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ENDDO


! Check that gas phase velocities are consistent.
      SELECT CASE (bc_plane(BCV))

      CASE ('W')
         IF(BC_U_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_U_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_U_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M)), '>'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('E')
         IF(BC_U_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_U_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_U_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M)), '<'
              CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('S')
         IF(BC_V_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_V_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_V_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M)), '>'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('N')
         IF(BC_V_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_V_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_V_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M)), '<'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('B')
         IF(BC_W_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_W_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_W_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M)), '>'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('T')
         IF(BC_W_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_W_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_W_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M)), '<'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      END SELECT

 1300 FORMAT('Error 1300: Invalid flow direction. ',A,' should be ',   &
         A,' zero. ',/'Please correct the mfix.dat file.')
      CALL FINL_ERR_MSG


      RETURN


 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_VEL_OUTFLOW
END MODULE set_bc_flow_module
