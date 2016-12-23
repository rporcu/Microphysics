MODULE CHECK_INITIAL_CONDITIONS_MODULE

! Parameter constants
      use param1, only: UNDEFINED, UNDEFINED_I, IS_DEFINED, IS_UNDEFINED, ZERO, ONE

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_INITIAL_CONDITIONS                                !
!  Author: P. Nicoletti                               Date: 02-DEC-91  !
!  Author: J.Musser                                   Date: 01-MAR-14  !
!                                                                      !
!  Purpose: check the initial conditions input section                 !
!     - check geometry of any specified IC region                      !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_INITIAL_CONDITIONS

! Global Variables:
!---------------------------------------------------------------------//
! Flag: IC geometry was detected.
      use ic, only: IC_DEFINED
! Flag: DEM solids present.
      use run, only: DEM_SOLIDS
! Flag: New run or a restart.
      use run, only: RUN_TYPE

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of IC.
      use param, only: DIMENSION_IC

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival

      use calc_cell_module, only: calc_cell
      use location_check_module, only: location_check

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter for ICs
      INTEGER :: ICV
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_INITIAL_CONDITIONS")

! Determine which ICs are DEFINED
      CALL CHECK_IC_GEOMETRY

! Loop over all IC arrays.
      DO ICV=1, DIMENSION_IC

! Verify user input for defined defined IC.
         IF(IC_DEFINED(ICV)) THEN
! Gas phase checks.
            CALL CHECK_IC_GAS_PHASE(ICV)
! Generic solids phase checks.
            CALL CHECK_IC_SOLIDS_PHASES(ICV)

! Verify that no data was defined for unspecified IC. ICs are only
! defined for new runs, so these checks are restricted to new runs.
         ELSEIF(RUN_TYPE == 'NEW') THEN
            CALL CHECK_IC_OVERFLOW(ICV)
         ENDIF
      ENDDO

! Check the initial conditions for the DEM model as well
      IF(DEM_SOLIDS) CALL CHECK_IC_COMMON_DISCRETE

! Finalize the error manager.
      CALL FINL_ERR_MSG

      RETURN

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_GEOMETRY                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_IC_GEOMETRY


! Global Variables:
!---------------------------------------------------------------------//
! Flag: IC contains geometric data and/or specified type
      use ic, only: IC_DEFINED
! Flag: IC type.
      use ic, only: IC_TYPE
! User specified: IC geometry
      use ic, only: IC_X_e, IC_X_w, IC_I_e, IC_I_w
      use ic, only: IC_Y_n, IC_Y_s, IC_J_n, IC_J_s
      use ic, only: IC_Z_t, IC_Z_b, IC_K_t, IC_K_b
! User specified: System geometry
      use geometry, only: IMAX, IMIN1, IMAX1, DX
      use geometry, only: JMAX, JMIN1, JMAX1, DY
      use geometry, only: KMAX, KMIN1, KMAX1, DZ
! Flag: New run or a restart.
      use run, only: RUN_TYPE

! Global Parameters:
!---------------------------------------------------------------------//
! The max number of ICs.
      use param, only: DIMENSION_IC

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival


      implicit none


! Local Variables:
!---------------------------------------------------------------------//
! Loop/variable indices
      INTEGER :: ICV
! Local spatial indices.
      INTEGER :: I_w, I_e, J_s, J_n, K_b, K_t
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_GEOMETRY")

! Check geometry of any specified IC region
      DO ICV = 1, DIMENSION_IC

         IC_DEFINED(ICV) = .FALSE.
         IF (IS_DEFINED(IC_X_W(ICV))) IC_DEFINED(ICV) = .TRUE.
         IF (IS_DEFINED(IC_X_E(ICV))) IC_DEFINED(ICV) = .TRUE.
         IF (IS_DEFINED(IC_Y_S(ICV))) IC_DEFINED(ICV) = .TRUE.
         IF (IS_DEFINED(IC_Y_N(ICV))) IC_DEFINED(ICV) = .TRUE.
         IF (IS_DEFINED(IC_Z_B(ICV))) IC_DEFINED(ICV) = .TRUE.
         IF (IS_DEFINED(IC_Z_T(ICV))) IC_DEFINED(ICV) = .TRUE.
         IF (IS_DEFINED(IC_I_W(ICV))) IC_DEFINED(ICV) = .TRUE.
         IF (IS_DEFINED(IC_I_E(ICV))) IC_DEFINED(ICV) = .TRUE.
         IF (IS_DEFINED(IC_J_S(ICV))) IC_DEFINED(ICV) = .TRUE.
         IF (IS_DEFINED(IC_J_N(ICV))) IC_DEFINED(ICV) = .TRUE.
         IF (IS_DEFINED(IC_K_B(ICV))) IC_DEFINED(ICV) = .TRUE.
         IF (IS_DEFINED(IC_K_T(ICV))) IC_DEFINED(ICV) = .TRUE.

! An IC is defined for restart runs only if it is a 'PATCH'.
         IF(RUN_TYPE /= 'NEW' .AND. IC_TYPE(ICV) /= 'PATCH') &
            IC_DEFINED(ICV) = .FALSE.

! Ignore patched IC regions for new runs. It may be better to flag this as
! and error to avoid user confusion.
         IF(RUN_TYPE == 'NEW' .AND. IC_TYPE(ICV) == 'PATCH') &
            IC_DEFINED(ICV) = .FALSE.


         IF(.NOT.IC_DEFINED(ICV)) CYCLE

         IF (IS_UNDEFINED(IC_X_W(ICV)) .AND. IS_UNDEFINED(IC_I_W(ICV))) THEN
            WRITE(ERR_MSG, 1100) ICV, 'IC_X_w and IC_I_w'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (IS_UNDEFINED(IC_X_E(ICV)) .AND. IS_UNDEFINED(IC_I_E(ICV))) THEN
            WRITE(ERR_MSG, 1100) ICV, 'IC_X_e and IC_I_e'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (IS_UNDEFINED(IC_Y_S(ICV)) .AND. IS_UNDEFINED(IC_J_S(ICV))) THEN
            WRITE(ERR_MSG, 1100) ICV, 'IC_Y_s and IC_J_s'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (IS_UNDEFINED(IC_Y_N(ICV)) .AND. IS_UNDEFINED(IC_J_N(ICV))) THEN
            WRITE(ERR_MSG, 1100) ICV, 'IC_Y_n and IC_J_n'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (IS_UNDEFINED(IC_Z_B(ICV)) .AND. IS_UNDEFINED(IC_K_B(ICV))) THEN
            WRITE(ERR_MSG, 1100) ICV, 'IC_Z_b and IC_K_b'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (IS_UNDEFINED(IC_Z_T(ICV)) .AND. IS_UNDEFINED(IC_K_T(ICV))) THEN
            WRITE(ERR_MSG, 1100) ICV, 'IC_Z_t and IC_K_t'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

      ENDDO   ! end loop over (icv = 1,dimension_ic)

 1100 FORMAT('Error 1100: Initial condition region ',I3,' is ill-',    &
         'defined.',/' > ',A,' are not specified.',/'Please correct ', &
         'the mfix.dat file.')


      DO ICV = 1, DIMENSION_IC

! Skip this check if the IC region is not specified.
         IF(.NOT.IC_DEFINED(ICV)) CYCLE

         IF (IS_DEFINED(IC_X_W(ICV)).AND. IS_DEFINED(IC_X_E(ICV))) THEN
            CALL CALC_CELL (IC_X_W(ICV), DX, IMAX, I_W)
            I_W = I_W + 1
            CALL CALC_CELL (IC_X_E(ICV), DX, IMAX, I_E)
            IF (IC_I_W(ICV)/=UNDEFINED_I .OR. IC_I_E(ICV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK (IC_I_W(ICV), I_W, ICV, 'IC - west')
               CALL LOCATION_CHECK (IC_I_E(ICV), I_E, ICV, 'IC - east')
            ELSE
               IC_I_W(ICV) = I_W
               IC_I_E(ICV) = I_E
            ENDIF
         ENDIF

! Report problems with calculated bounds.
         IF(IC_I_W(ICV) > IC_I_E(ICV)) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_W > IC_I_E'
             write(*,*)' dump:',IC_I_W(ICV),IC_I_E(ICV)
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_W(ICV) < IMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_W < IMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_W(ICV) > IMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_W > IMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_E(ICV) < IMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_E < IMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_E(ICV) > IMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_Z_t and IC_K_t'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (IS_DEFINED(IC_Y_S(ICV)) .AND. IS_DEFINED(IC_Y_N(ICV))) THEN
            CALL CALC_CELL (IC_Y_S(ICV), DY, JMAX, J_S)
            J_S = J_S + 1
            CALL CALC_CELL (IC_Y_N(ICV), DY, JMAX, J_N)
            IF (IC_J_S(ICV)/=UNDEFINED_I .OR. IC_J_N(ICV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK (IC_J_S(ICV), J_S, ICV, 'IC - south')
               CALL LOCATION_CHECK (IC_J_N(ICV), J_N, ICV, 'IC - north')
            ELSE
               IC_J_S(ICV) = J_S
               IC_J_N(ICV) = J_N
            ENDIF
         ENDIF

         IF(IC_J_S(ICV) > IC_J_N(ICV)) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_S > IC_J_N'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_S(ICV)<JMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_S < JMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_S(ICV)>JMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_S >  JMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_N(ICV)<JMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_N < JMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_N(ICV)>JMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_N > JMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


         IF (IS_DEFINED(IC_Z_B(ICV)) .AND. IS_DEFINED(IC_Z_T(ICV))) THEN
            CALL CALC_CELL (IC_Z_B(ICV), DZ, KMAX, K_B)
            K_B = K_B + 1
            CALL CALC_CELL (IC_Z_T(ICV), DZ, KMAX, K_T)
            IF (IC_K_B(ICV)/=UNDEFINED_I .OR. IC_K_T(ICV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK (IC_K_B(ICV), K_B, ICV, 'IC - bottom')
               CALL LOCATION_CHECK (IC_K_T(ICV), K_T, ICV, 'IC - top')
            ELSE
               IC_K_B(ICV) = K_B
               IC_K_T(ICV) = K_T
            ENDIF
         ENDIF

         IF(IC_K_B(ICV) > IC_K_T(ICV)) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_B > IC_K_T'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_B(ICV) < KMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_B < KMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_B(ICV) > KMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_B > KMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_T(ICV) < KMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_T < KMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_T(ICV) > KMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_T > KMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


 1101 FORMAT('Error 1101: Initial condition region ',I2,' is ill-',    &
         'defined.',/3x,A,/'Please correct the mfix.dat file.')

      ENDDO   ! end loop over (icv=1,dimension_ic)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_IC_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_GAS_PHASE                                       !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Verify gas phase input variables in IC region.              !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_IC_GAS_PHASE(ICV)

! Global Variables:
!---------------------------------------------------------------------//
! Gas phase volume fraction, pressure, temperature
      use ic, only: IC_EP_g, IC_P_g
! Gas phase velocity components.
      use ic, only: IC_U_g, IC_V_g, IC_W_g
! IC Type: UNDEFINED or PATCH.
      use ic, only: IC_TYPE

! Specified constant gas density and viscosity.
      use fld_const, only: ro_g0

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival

      implicit none

! Dummy Arguments:
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: ICV

! Local Variables:
!---------------------------------------------------------------------//
! Flag for regular IC regions
      LOGICAL :: BASIC_IC
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_GAS_PHASE")

! Patch ICs skip various checks.
      BASIC_IC = (IC_TYPE(ICV) /= 'PATCH')

! Check that gas phase velocity components are initialized.
      IF(BASIC_IC) THEN
         IF(IS_UNDEFINED(IC_U_G(ICV))) THEN
            WRITE(ERR_MSG, 1000) trim(iVar('IC_U_g',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(IS_UNDEFINED(IC_V_G(ICV))) THEN
            WRITE(ERR_MSG, 1000) trim(iVar('IC_V_g',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(IS_UNDEFINED(IC_W_G(ICV))) THEN
            WRITE(ERR_MSG, 1000) trim(iVar('IC_W_g',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

! Check that gas phase void fraction is initialized. Patched ICs may
! have an undefined volume fration. A second check is preformed on
! the solids.
      IF(IS_UNDEFINED(IC_EP_G(ICV)) .AND. BASIC_IC) THEN
         WRITE(ERR_MSG, 1000) trim(iVar('IC_EP_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Check that if the gas phase pressure is initialized and the gas is
! compressible that the gas phase pressure is not zero or negative
      IF(IS_DEFINED(IC_P_G(ICV))) THEN
         IF(IS_UNDEFINED(RO_G0).AND. IC_P_G(ICV)<=ZERO) THEN
            WRITE(ERR_MSG, 1100) trim(iVar('IC_P_g',ICV)),             &
               iVal(IC_P_G(ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

 1100 FORMAT('Error 1100: Pressure must be greater than 0.0 for ',     &
         'compressible flow',/'Illegal value: ',A,' = ',A,/'Please ',  &
         'correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      END SUBROUTINE CHECK_IC_GAS_PHASE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_IC_SOLIDS_PHASES                                  !
!  Author: P. Nicoletti                               Date: 02-DEC-91  !
!  Author: J.Musser                                   Date: 01-MAR-14  !
!                                                                      !
!  Purpose: Verify solids phase(s) input variables in IC region.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_IC_SOLIDS_PHASES(ICV)


! Global Variables:
!---------------------------------------------------------------------//
! Solids volume fraction, bulk density
      use ic, only: IC_EP_s, IC_ROP_s
! Solids velocity components.
      use ic, only: IC_U_s, IC_V_s, IC_W_s
! Gas phase volume fraction and temperature.
      use ic, only: IC_EP_g
! IC Type: UNDEFINED or PATCH.
      use ic, only: IC_TYPE

! Number of solids phases.
      use constant, only: MMAX

      use constant, only: ro_s0

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids species.
      use param, only: DIM_M

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival

      use toleranc, only: compare

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Index of IC region.
      INTEGER, INTENT(IN) :: ICV
! Loop/variable index
      INTEGER :: M
! Various sums.
      DOUBLE PRECISION SUM_EP
! Solids phase density in IC region.
      DOUBLE PRECISION :: IC_ROs(1:DIM_M)
! Flag to skip checks on indexed solid phase.
      LOGICAL :: SKIP(1:DIM_M)
! Flag for PATCH IC regions
      LOGICAL :: BASIC_IC
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_SOLIDS_PHASES")

! Patch ICs skip various checks.
      BASIC_IC = (IC_TYPE(ICV) /= 'PATCH')

! Calculate EP_s from EP_g if there is only one solids phase.
      IF(MMAX == 1 .AND. IS_DEFINED(IC_EP_S(ICV,1))) THEN
         IF(IS_DEFINED(IC_EP_g(ICV))) IC_EP_S(ICV,1) = ONE-IC_EP_g(ICV)
      ENDIF

! Bulk density or solids volume fraction must be explicitly defined
! if there are more than one solids phase.
      IF(MMAX > 1 .AND. .NOT.COMPARE(IC_EP_g(ICV),ONE)) THEN
! IC_EP_g may be undefined for PATCH IC regions.
         IF(IS_DEFINED(IC_EP_g(ICV))) THEN
            DO M = 1, MMAX
               IF(IS_DEFINED(IC_ROP_S(ICV,M)) .AND. &
                  IS_DEFINED(IC_EP_S(ICV,M))) THEN
                  WRITE(ERR_MSG, 1400) M, ICV, 'IC_ROP_s and IC_EP_s'
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDDO

! If IC_EP_G is undefined, then ROP_s and EP_s should be too.
         ELSE
            DO M = 1, MMAX
               IF(IS_DEFINED(IC_ROP_S(ICV,M)) .AND. &
                  IS_DEFINED(IC_EP_S(ICV,M))) THEN
                  WRITE(ERR_MSG, 1400) M, ICV, 'IC_ROP_s and IC_EP_s'
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

 1400 FORMAT('Error 1400: Insufficient solids phase ',I2,' ',          &
         'information for IC',/'region ',I3,'. ',A,' not specified.',/ &
         'Please correct the mfix.dat file.')

! Determine which solids phases are present.
      DO M = 1, MMAX
         SKIP(M)=(IS_UNDEFINED(IC_ROP_S(ICV,M)).OR.ABS(IC_ROP_S(ICV,M))<EPSILON(ZERO)) &
            .AND.(IS_UNDEFINED(IC_EP_S(ICV,M)) .OR.ABS(IC_EP_S(ICV,M))<EPSILON(ZERO))
      ENDDO

      IF(MMAX == 1 .AND. ABS(IC_EP_g(ICV)-ONE)>ZERO) SKIP(1) = .FALSE.

      DO M=1, MMAX

! check that solids phase m velocity components are initialized
         IF(BASIC_IC) THEN
            IF(IS_UNDEFINED(IC_U_S(ICV,M))) THEN
               IF (SKIP(M)) THEN
                  IC_U_S(ICV,M) = ZERO
               ELSE
                  WRITE(ERR_MSG, 1000)trim(iVar('IC_U_s',ICV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDIF

            IF(IS_UNDEFINED(IC_V_S(ICV,M))) THEN
               IF(SKIP(M)) THEN
                  IC_V_S(ICV,M) = ZERO
               ELSE
                  WRITE(ERR_MSG, 1000)trim(iVar('IC_V_s',ICV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDIF

            IF(IS_UNDEFINED(IC_W_S(ICV,M))) THEN
               IF(SKIP(M)) THEN
                  IC_W_S(ICV,M) = ZERO
               ELSE
                  WRITE(ERR_MSG, 1000)trim(iVar('IC_W_s',ICV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDIF

         ENDIF


         IC_ROs(M) = RO_s0(M)

      ENDDO   ! end loop over (m=1)


! Initialize the sum of the total volume fraction.
      SUM_EP = IC_EP_G(ICV)

      DO M=1, MMAX

! Clear out both varaibles if this phase is skipped.
         IF(BASIC_IC .AND. SKIP(M)) THEN
            IC_EP_S(ICV,M)  = ZERO
            IC_ROP_S(ICV,M) = ZERO

! Leave everything undefined for PATCH ICs that are not specifed.
         ELSEIF(.NOT.BASIC_IC .AND. (IS_UNDEFINED(IC_ROP_S(ICV,M))      &
            .AND. IS_UNDEFINED(IC_EP_S(ICV,M)))) THEN

! If both input parameters are defined. Make sure they are equivalent.
         ELSEIF(IS_DEFINED(IC_ROP_S(ICV,M)) .AND.                     &
            IS_DEFINED(IC_EP_S(ICV,M))) THEN

            IF(.NOT.COMPARE(IC_EP_S(ICV,M)*IC_ROs(M),                  &
               IC_ROP_S(ICV,M))) THEN

! BASIC_IC regions require that the IC_ROP_s and IC_EP_s specifications
! match although it is unlikely that anyone would specify both.
               IF(BASIC_IC) THEN

                  WRITE(ERR_MSG,1406) M, ICV
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

1406 FORMAT('Error 1406: IC_EP_s and IC_ROP_s are inconsistent for ',&
         'phase ',I2,/,'in IC region ', I3,'. Please correct the ',&
         'mfix.dat file.')


! PachedeIC regions defer to IC_EP_s if the values do not match. This
! prevents a dead lock or the need to define both. This case is rather
! common as a defined IC_EP_s is converted to IC_ROP_s. Therefore, if
! a patch region is used more than once, these values may not match.
               ELSE

                  WRITE(ERR_MSG,1407) trim(iVar('IC_ROP_s',ICV,M)), &
                     trim(iVAL(IC_ROP_S(ICV,M))), trim(iVar('IC_EP_s',&
                     ICV,M)), trim(iVAL(IC_EP_S(ICV,M)))
                  CALL FLUSH_ERR_MSG()

1407 FORMAT('Warning 1407: IC_EP_s and IC_ROP_s are inconsistent:',    &
         2(/3x,A,' = ',A),/'Deferring to IC_EP_s to overcome conflict.')

                  IC_ROP_S(ICV,M) = IC_EP_S(ICV,M)*IC_ROs(M)

               ENDIF
            ENDIF


! Compute IC_EP_s from IC_ROP_s
         ELSEIF(IS_UNDEFINED(IC_EP_S(ICV,M)))THEN
            IC_EP_S(ICV,M) = IC_ROP_S(ICV,M) / IC_ROs(M)

! Compute IC_ROP_s from IC_EP_s and IC_ROs
         ELSEIF(IS_UNDEFINED(IC_ROP_S(ICV,M))) THEN
            IC_ROP_S(ICV,M) = IC_EP_S(ICV,M) * IC_ROs(M)
! This is a sanity check.
         ELSE

         ENDIF
! Add this phase to the total volume fraction.
         SUM_EP = SUM_EP + IC_EP_S(ICV,M)
      ENDDO

! Verify that the volume fractions sum to one.
      IF(BASIC_IC .AND. .NOT.COMPARE(SUM_EP,ONE)) THEN
         WRITE(ERR_MSG,1410) ICV
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1410 FORMAT('Error 1410: Illegal initial condition region : ',I3,/    &
         'Sum of volume fractions does NOT equal ONE. Please correct',/&
         'the mfix.dat file.')


      CALL FINL_ERR_MSG

      RETURN


 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      END SUBROUTINE CHECK_IC_SOLIDS_PHASES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_OVERFLOW                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Verify that no data was defined for unspecified IC.         !
!                                                                      !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!

      SUBROUTINE CHECK_IC_OVERFLOW(ICV)

! Global Variables:
!---------------------------------------------------------------------//
! Gas phase volume fraction, pressure, temperature, species.
      use ic, only: IC_EP_g
! Gas phase velocity components.
      use ic, only: IC_U_g, IC_V_g, IC_W_g
! Solids volume fraction, bulk density
      use ic, only: IC_ROP_s
! Solids velocity components.
      use ic, only: IC_U_s, IC_V_s, IC_W_s
! IC Type: UNDEFINED or PATCH.
      use ic, only: IC_TYPE

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases
      use param, only: DIM_M

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival

      implicit none

! Dummy Arguments:
!---------------------------------------------------------------------//
! IC region index.
      INTEGER, INTENT(IN) :: ICV

! Local Variables:
!---------------------------------------------------------------------//
! Loop counters
      INTEGER :: M
!......................................................................!

      IF (IC_TYPE(ICV) == 'PATCH') RETURN

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_OVERFLOW")

! GAS PHASE quantities
! -------------------------------------------->>>
      IF(IS_DEFINED(IC_U_G(ICV))) THEN
          WRITE(ERR_MSG, 1010) trim(iVar('IC_U_g',ICV))
          CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IS_DEFINED(IC_V_G(ICV))) THEN
         WRITE(ERR_MSG, 1010) trim(iVar('IC_V_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IS_DEFINED(IC_W_G(ICV))) THEN
         WRITE(ERR_MSG, 1010) trim(iVar('IC_W_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IS_DEFINED(IC_EP_G(ICV))) THEN
         WRITE(ERR_MSG, 1010) trim(iVar('IC_EP_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
! --------------------------------------------<<<


! SOLIDS PHASE quantities
! -------------------------------------------->>>
      DO M=1, DIM_M
         IF(IS_DEFINED(IC_ROP_S(ICV,M))) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_ROP_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IS_DEFINED(IC_U_S(ICV,M))) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_U_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IS_DEFINED(IC_V_S(ICV,M))) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_V_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IS_DEFINED(IC_W_S(ICV,M))) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_W_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO
! --------------------------------------------<<<

      CALL FINL_ERR_MSG
      RETURN

 1010 FORMAT('Error 1010: ',A,' specified in an undefined IC region')

      END SUBROUTINE CHECK_IC_OVERFLOW

END SUBROUTINE CHECK_INITIAL_CONDITIONS
END MODULE CHECK_INITIAL_CONDITIONS_MODULE
