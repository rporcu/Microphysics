MODULE CHECK_BC_INFLOW_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   use param1   , only: undefined, one, zero, is_undefined, is_defined

   CONTAINS
! !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! !                                                                      !
! ! Subroutine: CHECK_BC_INFLOW                                          !
! ! Author: J.Musser                                    Date: 01-Mar-14  !
! !                                                                      !
! ! Purpose: Provided a detailed error message on common inflow BC       !
! !                                                                      !
! !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!       SUBROUTINE CHECK_BC_INFLOW(M_TOT, SKIP, BCV)

! ! Modules
! !---------------------------------------------------------------------//
!       use param, only: dim_m
!       use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar
!       IMPLICIT NONE

! ! Dummy arguments
! !---------------------------------------------------------------------//
!       INTEGER, INTENT(in) :: BCV
!       INTEGER, INTENT(in) :: M_TOT
!       LOGICAL, INTENT(in) :: SKIP(DIM_M)

!       RETURN
!       END SUBROUTINE CHECK_BC_INFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_MASS_INFLOW                                     !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message on BC                     !
!                                                                      !
! Comments:                                                            !
!     The velocities at the inflow face are fixed and the momentum     !
!     equations are not solved in the inflow cells. Since the flow is  !
!     into the domain all other scalars that are used need to be       !
!     specified (e.g., mass fractions, void fraction, etc.,)           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_MASS_INFLOW(M_TOT, SKIP, BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_ep_g, bc_p_g
      use bc, only: bc_rop_s, bc_ep_s
      use bc, only: bc_massflow_g
      use param    , only: dim_m
      use fld_const, only: ro_g0
      use constant , only: ro_s0
      use toleranc , only: compare
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)

! Local variables
!---------------------------------------------------------------------//
! loop/variable indices
      INTEGER :: M
      real(c_real) :: SUM_EP
! Solids phase density in BC region.
      real(c_real) :: BC_ROs(DIM_M)

!---------------------------------------------------------------------//

      CALL INIT_ERR_MSG("CHECK_BC_MASS_INFLOW")

! Check gas phase volume fraction.
      IF(IS_UNDEFINED(BC_EP_G(BCV))) THEN
         WRITE(ERR_MSG, 1000) trim(iVar('BC_EP_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Verify compressible boundary condition variables.
      IF(IS_UNDEFINED(RO_G0)) THEN
         IF(IS_UNDEFINED(BC_P_G(BCV))) THEN
            IF(IS_UNDEFINED(BC_MASSFLOW_G(BCV)) .AND.                   &
               ABS(BC_MASSFLOW_G(BCV)) > ZERO) THEN
               WRITE(ERR_MSG, 1100) trim(iVar('BC_P_g',BCV))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
 1100 FORMAT('Error 1100: ',A,' must be specified for compressible ',  &
         'flows',/'when specifying BC_MASSFLOW_g to make the ',        &
         'conversion to velocity.',/'Please correct the mfix.dat file.')

         ELSEIF(BC_P_G(BCV) <= ZERO) THEN
            WRITE(ERR_MSG, 1101) BCV, trim(iVal(BC_P_G(BCV)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
 1101 FORMAT('Error 1101: Pressure must be greater than zero for ',    &
         'compressible flow',/' >>>  BC_P_g(',I3,') = ',A,/'Please ',  &
         'correct the mfix.dat file.')
      ENDIF

! Calculate the solids volume fraction from the gas phase if there is
! only one solids phase.
      IF(M_TOT == 1 .AND. IS_UNDEFINED(BC_EP_S(BCV,1))) THEN
         BC_EP_S(BCV,1) = ONE - BC_EP_g(BCV)
      ENDIF

! Bulk density or solids volume fraction must be explicitly defined
! if there are more than one solids phase.
      IF(M_TOT > 1 .AND. .NOT.COMPARE(BC_EP_g(BCV),ONE)) THEN
         DO M = 1, M_TOT
            IF(IS_UNDEFINED(BC_ROP_S(BCV,M)) .AND. &
               IS_UNDEFINED(BC_EP_S(BCV,M))) THEN
               WRITE(ERR_MSG, 1200) M, BCV, 'BC_ROP_s and BC_EP_s'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO
      ENDIF
 1200 FORMAT('Error 1200: Insufficient solids phase ',I2,' data ',     &
         'for BC',I3,'. ',/A,' not specified.',/'Please correct the ', &
         'mfix.dat file.')

! Initialize the sum of the total volume fraction.
      SUM_EP = BC_EP_G(BCV)
      DO M = 1, M_TOT

! If this phase is not present, clear out EPs and ROPs for the BC and
! cycle the solids loop. No need to continue checks.
         IF(SKIP(M)) THEN
            BC_EP_S(BCV,M)  = ZERO
            BC_ROP_S(BCV,M) = ZERO
            CYCLE
         ENDIF

! Set the solids density for the BC region.
         BC_ROs(M) = RO_s0(M)

! If both input parameters are defined. Make sure they are equivalent.
         IF(IS_DEFINED(BC_ROP_S(BCV,M)) .AND.                         &
            IS_DEFINED(BC_EP_S(BCV,M))) THEN

            IF(.NOT.COMPARE(BC_EP_S(BCV,M)*BC_ROs(M),                  &
               BC_ROP_S(BCV,M))) THEN
               WRITE(ERR_MSG,1214) BCV
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
 1214 FORMAT('Error 1214: Illegal initial condition region : ',I3,/    &
         'BC_EP_s and BC_ROP_s are inconsistent. Please correct the ',/&
         'mfix.dat file.')

! Compute BC_EP_s from BC_ROP_s
         ELSEIF(IS_UNDEFINED(BC_EP_S(BCV,M))) THEN
            BC_EP_S(BCV,M) = BC_ROP_S(BCV,M) / BC_ROs(M)

! Compute BC_ROP_s from BC_EP_s and BC_ROs
         ELSEIF(IS_UNDEFINED(BC_ROP_S(BCV,M))) THEN
            BC_ROP_S(BCV,M) = BC_EP_S(BCV,M) * BC_ROs(M)

         ENDIF
! Add this phase to the total volume fraction.
         SUM_EP = SUM_EP + BC_EP_S(BCV,M)
      ENDDO

! Verify that the volume fractions sum to one.
      IF(.NOT.COMPARE(SUM_EP,ONE)) THEN
         WRITE(ERR_MSG,1215) BCV, trim(iVal(SUM_EP))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1215 FORMAT('Error 1215: Illegal boundary condition region: ',I3,'. ',&
         'Sum of volume',/'fractions does NOT equal ONE. (SUM = ',A,   &
         ')',/'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_MASS_INFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_P_INFLOW                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided detailed error message on bc                       !
!                                                                      !
! Comments:                                                            !
!     Unlike the MI boundary, for the PI boundary the velocities at    !
!     the inflow face are calculated by solving the momentum eqns      !
!     and are not fixed. In this way, the PI is similar to the PO      !
!     except that the flow is into the domain and hence all other      !
!     scalars (e.g., mass fractions, void fraction, temperature,       !
!     etc.,) at the inflow cells need to be specified. To satisfy      !
!     the error routines at the start of the simulation, both the      !
!     tangential and normal components at the inflow also need to      !
!     be specified. The velocities values essentially serve as IC.     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_P_INFLOW(M_TOT, SKIP, BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_p_g, bc_rop_s
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: bc_u_s, bc_v_s, bc_w_s
      use param    , only: dim_m
      use fld_const, only: ro_g0
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: M

!---------------------------------------------------------------------//

      CALL INIT_ERR_MSG("CHECK_BC_P_INFLOW")

! Remove checks on bc_ep_g/bc_rop_s; using routine check_bc_outflow

      IF (IS_UNDEFINED(BC_P_G(BCV))) THEN
         WRITE(ERR_MSG,1000) 'BC_P_g', BCV
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      ELSEIF (BC_P_G(BCV)<=ZERO .AND. IS_UNDEFINED(RO_G0)) THEN
         WRITE(ERR_MSG, 1101) BCV, trim(iVal(BC_P_G(BCV)))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1101 FORMAT('Error 1101: Pressure must be greater than zero for ',    &
         'compressible flow',/' >>>  BC_P_g(',I3,') = ',A,/'Please ',  &
         'correct the mfix.dat file.')
      ENDIF

! Check that velocities are also specified. These are essentially used
! as initial conditions for the boundary region. If they are not
! specified then a default value is set here.
      IF(IS_UNDEFINED(BC_U_G(BCV))) THEN
         BC_U_G(BCV) = ZERO
         WRITE(ERR_MSG, 1300) trim(iVar('BC_U_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
      ENDIF

      IF(IS_UNDEFINED(BC_V_G(BCV))) THEN
         BC_V_G(BCV) = ZERO
         WRITE(ERR_MSG, 1300) trim(iVar('BC_V_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
      ENDIF

      IF(IS_UNDEFINED(BC_W_G(BCV))) THEN
         BC_W_G(BCV) = ZERO
         WRITE(ERR_MSG, 1300) trim(iVar('BC_W_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
      ENDIF

      DO M = 1, M_TOT
         IF (SKIP(M)) THEN
            BC_U_S(BCV,M) = ZERO
            BC_V_S(BCV,M) = ZERO
            BC_W_S(BCV,M) = ZERO
         ELSE
            IF(IS_UNDEFINED(BC_U_S(BCV,M))) THEN
               BC_U_S(BCV,M) = ZERO
               IF(ABS(BC_ROP_S(BCV,M)) > ZERO) THEN
                  WRITE(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
               ENDIF
            ENDIF

            IF(IS_UNDEFINED(BC_V_S(BCV,M))) THEN
               BC_V_S(BCV,M) = ZERO
               IF(ABS(BC_ROP_S(BCV,M)) > ZERO) THEN
                  WRITE(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
               ENDIF
            ENDIF

            IF(IS_UNDEFINED(BC_W_S(BCV,M))) THEN
               BC_W_S(BCV,M) = ZERO
               IF(ABS(BC_ROP_S(BCV,M)) > ZERO) THEN
                  WRITE(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

 1300 FORMAT('Warning 1300: ',A,' was undefined. This variable was ', &
         'set ',/ 'to zero to be used as the initial value in the BC ',&
         'region.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_BC_P_INFLOW
END MODULE CHECK_BC_INFLOW_MODULE
