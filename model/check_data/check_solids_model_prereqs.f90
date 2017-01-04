MODULE CHECK_SOLIDS_MODEL_PREREQS_MODULE
CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_MODEL_PREREQS                              !
!  Purpose: Check the distributed parallel namelist variables.         !
!                                                                      !
!  Author: P. Nicoletti                               Date: 14-DEC-99  !
!  Reviewer: J.Musser                                 Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_MODEL_PREREQS


! Global Variables:
!---------------------------------------------------------------------//
! Number of ranks.
      use run, only: SOLIDS_MODEL

! Flag: Use DES E-E solids model
      use run, only: TFM_SOLIDS, TFM_COUNT
      use run, only: DEM_SOLIDS, DEM_COUNT
      use run, only: PIC_SOLIDS, PIC_COUNT

! Flag: Use DES E-L model

! Flag: gas/solids E-L simulation, otherwise granular flow.
      use discretelement, only: DES_CONTINUUM_COUPLED
! Flag: Fluid affects particles, but particles do not impact fluid.
      use discretelement, only: DES_ONEWAY_COUPLED
! Number of phases specified by the user.
      use constant, only: MMAX
! User specified constant gas density
      use fld_const, only: ro_g0
! Print E-L data.
      use discretelement, only: PRINT_DES_DATA

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases.
      use param, only: DIM_M

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, &
         init_err_msg, ivar, ival, err_msg

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: M ! Phase index
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_MODEL_PREREQS")

! Loop over the phases to see what was specified.
      DO M=1, DIM_M
         SOLIDS_MODEL(M) = trim(adjustl(SOLIDS_MODEL(M)))
         SELECT CASE(SOLIDS_MODEL(M))
         CASE ('TFM'); TFM_COUNT = TFM_COUNT + 1
         CASE ('DEM'); DEM_COUNT = DEM_COUNT + 1
         CASE ('PIC'); PIC_COUNT = PIC_COUNT + 1
         CASE ('---')
         CASE DEFAULT
            WRITE(ERR_MSG,1001) iVar('SOLIDS_MODEL',M), SOLIDS_MODEL(M)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1001 FORMAT('Error 1001: Unknown solids model: ',A,' = ',A)

         END SELECT
      ENDDO

! Clear out the unused phases.
      MMAX =  TFM_COUNT + DEM_COUNT + PIC_COUNT

! Set the runtime flags:
      TFM_SOLIDS = (TFM_COUNT > 0)
      DEM_SOLIDS = (DEM_COUNT > 0)
      PIC_SOLIDS = (PIC_COUNT > 0)

      IF(TFM_SOLIDS)THEN
         WRITE(ERR_MSG, 1002)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(PIC_SOLIDS)THEN
         WRITE(ERR_MSG, 1003)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1002 FORMAT('Error 1002: TFM solids are not supported in this&
         & version of MFIX.',/'Please correct the mfix.dat file.')

 1003 FORMAT('Error 1003: PIC solids are not supported in this&
         & version of MFIX.',/'Please correct the mfix.dat file.')

! Set flag for coupled simulations
      DES_CONTINUUM_COUPLED = dem_solids .and. (ABS(RO_g0) > 0.0d0)

! Overwrite user settings if no Lagrangian solids
      IF(.NOT.DEM_SOLIDS) THEN
         DES_CONTINUUM_COUPLED = .FALSE.
         PRINT_DES_DATA = .FALSE.
         DES_ONEWAY_COUPLED = .FALSE.
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_MODEL_PREREQS
END MODULE CHECK_SOLIDS_MODEL_PREREQS_MODULE
