MODULE FLOW_TO_VEL_NEW_MODULE

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int

      use bc, only: BC_MASSFLOW_G
      use bc, only: BC_MASSFLOW_S
      use bc, only: BC_VOLFLOW_G
      use bc, only: BC_VOLFLOW_S
      use bc, only: bc_type, bc_plane, bc_u_s, bc_v_s, bc_w_s, bc_ep_s, bc_ep_g
      use bc, only: bc_area, bc_u_g, bc_v_g, bc_w_g
      use exit_mod, only: mfix_exit
      use param, only: DIM_M
      use param1, only: UNDEFINED, ZERO, ONE, IS_DEFINED, IS_UNDEFINED
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar
      use toleranc, only: compare

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: FLOW_TO_VEL_NEW                                         !
!  Author: M. Syamlal                                 Date: 28-JUL-92  !
!                                                                      !
!  Purpose: Convert volumetric and mass flow rates to velocities       !
!     A specified mass flow rate is first converted to volumetric      !
!     flow rate. The volumetric flow rate is then converted to a       !
!     velocity.                                                        !
!                                                                      !
!    When both flow rates and velocities are specified, a consistency  !
!    check is done. The first time flow_to_vel is called in by setting !
!    the logical DO_VEL_CHECK to .TRUE.. If cut-cells are not used,    !
!    flow_to_vel is only called once.  When cut-cells are used,        !
!    flow_to_vel is called another time after the cut-cell pre-        !
!    processing stage. During, the second call, the velocity check     !
!    should not be performed, because the velocity assigned suring the !
!    first call will not match the flow rate. Therfore, when called    !
!    from cut_cell_preprocessing.f DO_VEL_CHECK is set to .FALSE..     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE FLOW_TO_VEL_NEW(DO_VEL_CHECK, M_TOT, SKIP, BCV)

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      LOGICAL, intent(in) :: DO_VEL_CHECK

! loop/variable indices
      INTEGER, intent(in) :: M_TOT
      LOGICAL, intent(in) :: SKIP(DIM_M)

! loop/variable indices
      INTEGER, intent(in) :: BCV

! Whether any volumetric flow conversion was done
      LOGICAL :: CONVERTED = .FALSE.

! Loop index
      INTEGER :: M


      CALL INIT_ERR_MSG("FLOW_TO_VEL_NEW")

! Mass flows rates are converted to volumetric flow rates.
      IF(IS_DEFINED(BC_MASSFLOW_G(BCV))) &
         CALL GAS_MASSFLOW_TO_VOLFLOW(BCV)

      DO M=1,M_TOT
         IF(IS_DEFINED(BC_MASSFLOW_S(BCV,M))) &
            CALL SOLIDS_MASSFLOW_TO_VOLFLOW(BCV,M,SKIP(M))
      ENDDO

! Volumetric flow rates are converted to velocities.
      IF(IS_DEFINED(BC_VOLFLOW_G(BCV))) THEN
         CALL GAS_VOLFLOW_TO_VELOCITY(DO_VEL_CHECK, BCV)
! Set the conversion flag.
         CONVERTED = .TRUE.
      ENDIF

      DO M=1,M_TOT
         IF(IS_DEFINED(BC_VOLFLOW_S(BCV,M))) THEN
            CALL SOLIDS_VOLFLOW_TO_VELOCITY(DO_VEL_CHECK,BCV,M,SKIP(M))
! Set the conversion flag.
            CONVERTED = .TRUE.
         ENDIF
      ENDDO

      IF(CONVERTED ) THEN
         WRITE(ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1100 FORMAT('Warning 1100: Some volumetric or mass flow rates have ', &
         'been converted',/'velocity. Ensure that the third (unused) ',&
         'dimension in 2D simulations',/'is correctly specified.')

      END SUBROUTINE FLOW_TO_VEL_NEW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_MASSFLOW_TO_VOLFLOW                                 !
!  Author: M. Syamlal                                 Date: 28-JUL-92  !
!                                                                      !
!  Purpose: Convert a gas phase BC input from a mass flow rate to      !
!  a volumetric flow rate.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_MASSFLOW_TO_VOLFLOW(BCV)

      use bc, only: BC_MASSFLOW_g
      use bc, only: BC_P_g
      use bc, only: BC_T_g
      use bc, only: BC_VOLFLOW_g
      use eos, only: EOSG
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar
      use fld_const, only: mw_avg, ro_g0
      use scales   , only: P_REF

      IMPLICIT NONE

      INTEGER, intent(in) :: BCV

! Volumetric flow rate computed from mass flow rate
      real(c_real) :: VOLFLOW
! Average molecular weight
      real(c_real) :: MW

      CALL INIT_ERR_MSG("GAS_MASSFLOW_TO_VOLFLOW")

! No need to convert if the mass flow is zero.
      IF(COMPARE(BC_MASSFLOW_G(BCV),ZERO)) THEN
         VOLFLOW = ZERO

! Incompressible gas BC.
      ELSEIF(IS_DEFINED(RO_G0)) THEN
         VOLFLOW = BC_MASSFLOW_G(BCV)/RO_G0

! Well-defined compresible gas BC.
      ELSEIF(IS_DEFINED(BC_P_G(BCV)) .AND. IS_DEFINED(BC_T_G(BCV))) THEN
         MW = mw_avg
         VOLFLOW = BC_MASSFLOW_G(BCV) / &
            EOSG(MW,(BC_P_G(BCV)-P_REF),BC_T_G(BCV))

! Fails. This shouldn't happen as previous checks should catch any
! errors leading to this routine.
      ELSE
         WRITE(ERR_MSG, 1100) BCV
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error 1100: Boundary condition ',I3,' failed sanity ',   &
      'check.',/'Please report this to the MFIX mailing list.')

      ENDIF

! Check that a specified volumetric flow matches the calculated value.
      IF(IS_DEFINED(BC_VOLFLOW_G(BCV))) THEN
         IF(.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_G(BCV))) THEN
            WRITE(ERR_MSG,1101) trim(iVar('BC_MASSFLOW_g',BCV)), BCV,  &
               VOLFLOW, BC_VOLFLOW_g(BCV)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ELSE

! Store the calculated volumetric flow rate.
         BC_VOLFLOW_G(BCV) = VOLFLOW
      ENDIF

 1101 FORMAT('Error 1101: Volumetric flow rate calculated from ',A,/   &
         'does NOT equal the specified volumetric flow rate for BC',I3,&
         /3x,'>>> Calculated: ',G14.7,/3x,'>>> Specified:  ',G14.7,/   &
         'Please correct the mfix.dat file.')

! Clean up and return.
      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE GAS_MASSFLOW_TO_VOLFLOW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOLIDS_MASSFLOW_TO_VOLFLOW                              !
!  Author: M. Syamlal                                 Date: 28-JUL-92  !
!                                                                      !
!  Purpose: Convert solids phase BC input from a mass flow rate to     !
!  a volumetric flow rate.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SOLIDS_MASSFLOW_TO_VOLFLOW(BCV,M, SKIP_M)

      USE bc, only: BC_MASSFLOW_s
      USE bc, only: BC_VOLFLOW_s
      USE constant, only: RO_s0
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

! loop/variable indices
      INTEGER, intent(in) :: BCV, M
      LOGICAL, intent(in) :: SKIP_M

! Volumetric flow rate computed from mass flow rate
      real(c_real) :: VOLFLOW

      CALL INIT_ERR_MSG("SOLIDS_MASSFLOW_TO_VOLFLOW")


      IF(SKIP_M) THEN
         WRITE(ERR_MSG,1100) M, BCV, trim(iVar("BC_MASSFLOW_S",BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Solids phase ',I2,' has a specified mass ',  &
         'flow rate',/'at BC ',I3,', ',A,'. But, both BC_ROP_s and ',&
         'BC_EP_s are zero or undefined.',/'Please correct the ',&
         'mfix.dat file.')

      VOLFLOW = BC_MASSFLOW_S(BCV,M)/RO_S0(M)

! If volumetric flow is also specified compare both
      IF(IS_DEFINED(BC_VOLFLOW_S(BCV,M))) THEN
         IF(.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_S(BCV,M))) THEN
            WRITE(ERR_MSG,1101) trim(iVar('BC_MASSFLOW_S',BCV,M)), BCV, &
               VOLFLOW, BC_VOLFLOW_S(BCV,M)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ELSE
         BC_VOLFLOW_S(BCV,M) = VOLFLOW
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1101 FORMAT('Error 1101: Volumetric flow rate calculated from ',A,/   &
         'does NOT equal the specified volumetric flow rate for BC',I3,&
         /3x,'>>> Calculated: ',G14.7,/3x,'>>> Specified:  ',G14.7,/   &
         'Please correct the mfix.dat file.')


      END SUBROUTINE SOLIDS_MASSFLOW_TO_VOLFLOW




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_VOLFLOW_TO_VELOCITY                                 !
!  Author: M. Syamlal                                 Date: 28-JUL-92  !
!                                                                      !
!  Purpose: Convert gas phase volumetric rate to a velocity.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_VOLFLOW_TO_VELOCITY(DO_VEL_CHECK, BCV)

      USE compar, only: myPE

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop/variable indices
      INTEGER, intent(in) :: BCV

! Whether any volumetric flow conversion was done
      LOGICAL, intent(in) :: DO_VEL_CHECK


      real(c_real) :: SGN, OFF

! Velocity computed from volumetric flow rate
      real(c_real) :: VEL
!-----------------------------------------------
      CALL INIT_ERR_MSG("GAS_VOLFLOW_TO_VELOCITY")

      SELECT CASE (trim(BC_TYPE(BCV)))
      CASE ('MASS_INFLOW');  SGN =  ONE; OFF = ZERO
      CASE ('MASS_OUTFLOW'); SGN = -ONE; OFF = ONE
      CASE DEFAULT
        write(*,*) 'error in GAS_VOLFLOW_TO_VELOCITY'
        call mfix_exit(myPE)
      END SELECT

      SELECT CASE (BC_PLANE(BCV))
      CASE ('W'); SGN = -SGN
      CASE ('S'); SGN = -SGN
      CASE ('B'); SGN = -SGN
      END SELECT

! Calculate the velocity based on the volumetric flow rate,
! BC area and BC volume fraction.
      VEL = SGN*BC_VOLFLOW_G(BCV)/(BC_AREA(BCV)*BC_EP_G(BCV))

! if the user also defined the boundary velocity through the plane, then
! check that the calculated value agrees with the specified value. if
! the user did not define the boundary velocity through the plane, then
! if mass_inflow set the value of the boundary velocity to the
! calculated value. otherwise do nothing.
      IF(BC_PLANE(BCV) == 'W' .OR. BC_PLANE(BCV)== 'E') THEN

         IF(IS_DEFINED(BC_U_G(BCV)) .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL,BC_U_G(BCV))) THEN
               WRITE(ERR_MSG,1100) BCV, VEL, 'BC_U_g', BC_U_G(BCV)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_U_G(BCV) = VEL
            BC_V_G(BCV) = OFF * BC_V_G(BCV)
            BC_W_G(BCV) = OFF * BC_W_G(BCV)
         ENDIF

      ELSEIF(BC_PLANE(BCV) == 'S' .OR. BC_PLANE(BCV)== 'N') THEN
         IF(IS_DEFINED(BC_V_G(BCV)) .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL,BC_V_G(BCV))) THEN
               WRITE(ERR_MSG, 1100) BCV, VEL, 'BC_V_g', BC_V_G(BCV)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_V_G(BCV) = VEL
            BC_U_G(BCV) = OFF * BC_U_G(BCV)
            BC_W_G(BCV) = OFF * BC_W_G(BCV)
         ENDIF

      ELSEIF(BC_PLANE(BCV) == 'B' .OR. BC_PLANE(BCV)== 'T') THEN
         IF(IS_DEFINED(BC_W_G(BCV)) .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL, BC_W_G(BCV))) THEN
               WRITE(ERR_MSG, 1100) BCV, VEL, 'BC_W_g', BC_W_G(BCV)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_W_G(BCV) = VEL
            BC_U_G(BCV) = OFF * BC_U_G(BCV)
            BC_V_G(BCV) = OFF * BC_V_G(BCV)
         ENDIF

      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1100 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,') = ',G14.7,/1X,70('*')/)

      END SUBROUTINE GAS_VOLFLOW_TO_VELOCITY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOLIDS_VOLFLOW_TO_VELOCITY                              !
!  Author: M. Syamlal                                 Date: 28-JUL-92  !
!                                                                      !
!  Purpose: Convert volumetric and mass flow rates to velocities       !
!     A specified mass flow rate is first converted to volumetric      !
!     flow rate. The volumetric flow rate is then converted to a       !
!     velocity.                                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SOLIDS_VOLFLOW_TO_VELOCITY(DO_VEL_CHECK, BCV, M, SKIP_M)

      USE funits, only: dmp_log, unit_log
      USE compar, only: myPE

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop/variable indices
      INTEGER, intent(in) :: BCV, M
! Logical to preform velocity check.
      LOGICAL, intent(in) :: DO_VEL_CHECK, SKIP_M

! Velocity computed from volumetric flow rate
      real(c_real) :: VEL
      real(c_real) :: SGN, OFF

!-----------------------------------------------

      CALL INIT_ERR_MSG("SOLIDS_VOLFLOW_TO_VELOCITY")

      IF(SKIP_M) THEN
         WRITE(ERR_MSG,1100) M, BCV, trim(iVar("BC_VOLFLOW_S",BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Solids phase ',I2,' has a specified ',       &
         'volumetric flow rate',/'at BC ',I3,', ',A,'. But, both ',&
         'BC_ROP_s and BC_EP_s are zero or undefined.',/'Please ',&
         'the mfix.dat file.')

      SELECT CASE (trim(BC_TYPE(BCV)))
      CASE ('MASS_INFLOW');  SGN =  ONE; OFF = ZERO
      CASE ('MASS_OUTFLOW'); SGN = -ONE; OFF = ONE
      CASE DEFAULT
        write(*,*) 'error in SOLIDS_VOLFLOW_TO_VELOCITY'
        call mfix_exit(myPE)
      END SELECT

      SELECT CASE (BC_PLANE(BCV))
      CASE ('W'); SGN = -SGN
      CASE ('S'); SGN = -SGN
      CASE ('B'); SGN = -SGN
      END SELECT

      IF(ABS(BC_EP_S(BCV,M)) > ZERO) THEN
         VEL = SGN * BC_VOLFLOW_S(BCV,M)/(BC_AREA(BCV)*BC_EP_S(BCV,M))
      ELSE
         IF(ABS(BC_VOLFLOW_S(BCV,M)) > ZERO) THEN
            VEL = ZERO
         ELSE
            IF(DMP_LOG)WRITE (UNIT_LOG, 1101) BCV, M
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

 1101 FORMAT('Error 1101: BC No:',I2,' Non-zero vol. or mass flow ',&
         'specified with BC_ROP_s', I1,' = 0.')

      IF(BC_PLANE(BCV) == 'W' .OR. BC_PLANE(BCV)== 'E') THEN
         IF(IS_DEFINED(BC_U_S(BCV,M)) .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL, BC_U_S(BCV,M))) THEN
              WRITE(ERR_MSG, 1300) BCV, (-VEL), 'BC_U_s', M, BC_U_S(BCV,M)
              CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_U_S(BCV,M) = VEL
            BC_V_S(BCV,M) = OFF * BC_V_S(BCV,M)
            BC_W_S(BCV,M) = OFF * BC_W_S(BCV,M)
         ENDIF

      ELSEIF(BC_PLANE(BCV) == 'S' .OR. BC_PLANE(BCV)== 'N') THEN
         IF(IS_DEFINED(BC_V_S(BCV,M)) .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL,BC_V_S(BCV,M))) THEN
               WRITE(ERR_MSG,1300) BCV, VEL, 'BC_V_s', M, BC_V_S(BCV,M)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_V_S(BCV,M) = VEL
            BC_U_S(BCV,M) = OFF * BC_U_S(BCV,M)
            BC_W_S(BCV,M) = OFF * BC_W_S(BCV,M)
         ENDIF

      ELSEIF(BC_PLANE(BCV) == 'B' .OR. BC_PLANE(BCV)== 'T') THEN
         IF(IS_DEFINED(BC_W_S(BCV,M)) .AND. DO_VEL_CHECK) THEN
            IF(.NOT.COMPARE(VEL,BC_W_S(BCV,M))) THEN
               WRITE(ERR_MSG, 1300) BCV, VEL, 'BC_W_s', M, BC_W_S(BCV,M)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            BC_W_S(BCV,M) = VEL
            BC_U_S(BCV,M) = OFF * BC_U_S(BCV,M)
            BC_V_S(BCV,M) = OFF * BC_V_S(BCV,M)
         ENDIF
      ENDIF


      CALL FINL_ERR_MSG

      RETURN

 1300 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,I1,') = ',G14.7,/1X,70('*')/)

      END SUBROUTINE SOLIDS_VOLFLOW_TO_VELOCITY

END MODULE FLOW_TO_VEL_NEW_MODULE
