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
      use discretelement, only: DISCRETE_ELEMENT
! Flag: Use MPPIC E-L model
      use mfix_pic, only: MPPIC
! Flag: Use TFM and DEM solids models.
      use discretelement, only: DES_CONTINUUM_HYBRID
! Flag: gas/solids E-L simulation, otherwise granular flow.
      use discretelement, only: DES_CONTINUUM_COUPLED
! Flag: Fluid affects particles, but particles do not impact fluid.
      use discretelement, only: DES_ONEWAY_COUPLED
! Number of discrete solids phases.
      use discretelement, only: DES_MMAX
! Number of phases specified by the user.
      use physprop, only: MMAX, SMAX
! User specified constant gas density
      use fldvar, only: RO_g0
! Print E-L data.
      use discretelement, only: PRINT_DES_DATA

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases.
      use param, only: DIM_M

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: M ! Phase index
! Error indicator for mmax
      LOGICAL :: ERR_MMAX
! Dummy integer
      INTEGER :: N
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_MODEL_PREREQS")

! Moved this here but will need to discuss if this is best route
! going forward
      ERR_MMAX = .FALSE.
      IF (MMAX < 0) ERR_MMAX = .TRUE.

! Check MMAX
      IF (ERR_MMAX) THEN
         WRITE(ERR_MSG, 1000) iVal(DIM_M)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1000 FORMAT('Error 1000: MMAX out of range. Min: 0, Max: ',A)

! Loop over the phases to see what was specified.
      DO M=1, MMAX
         SOLIDS_MODEL(M) = trim(adjustl(SOLIDS_MODEL(M)))
         SELECT CASE(SOLIDS_MODEL(M))
         CASE ('TFM'); TFM_COUNT = TFM_COUNT + 1
         CASE ('DEM'); DEM_COUNT = DEM_COUNT + 1
         CASE ('PIC'); PIC_COUNT = PIC_COUNT + 1

         CASE DEFAULT
            WRITE(ERR_MSG,1001) iVar('SOLIDS_MODEL',M), SOLIDS_MODEL(M)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1001 FORMAT('Error 1001: Unknown solids model: ',A,' = ',A)

         END SELECT
      ENDDO

! Set the number of discrete phases.
      DES_MMAX = DEM_COUNT + PIC_COUNT

! Set the number of TFM phases.  (should be equivalent to tfm_count)
      MMAX = MMAX - DES_MMAX
      SMAX = MMAX   ! USE smax in the code for the number of TFM phases

! Clear out the unused phases.
      SOLIDS_MODEL((MMAX+DES_MMAX)+1:DIM_M) = '---'

! Set the runtime flags:
      TFM_SOLIDS = (TFM_COUNT > 0)
      DEM_SOLIDS = (DEM_COUNT > 0)
      PIC_SOLIDS = (PIC_COUNT > 0)

! MPPIC and continuum solids don't mix.
      IF(PIC_SOLIDS .AND. TFM_SOLIDS)THEN
         WRITE(ERR_MSG, 1002)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1002 FORMAT('Error 1002: MPPIC solids and continuum solids cannot ',&
         'be combined.',/'Please correct the mfix.dat file.')


! MPPIC and DEM solids don't mix.
      IF(PIC_SOLIDS .AND. DEM_SOLIDS)THEN
         WRITE(ERR_MSG, 1003)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1003 FORMAT('Error 1003: MPPIC solids and DES solids cannot be ',     &
         'combined.',/'Please correct the mfix.dat file.')

! temporary move for now since these rely on definition of mmax/smax
      M = max(SMAX, DES_MMAX)

 1004 FORMAT('Error 1004: Illegal or unknown input: ',A,' = ',A,/  &
         'Please correct the mfix.dat file.')

! Set the DEM runtime flag.
      DISCRETE_ELEMENT = DEM_SOLIDS .OR. PIC_SOLIDS
! Set the MMPIC runtime flag.
      MPPIC = PIC_SOLIDS
! Set the Hybird flag.
      DES_CONTINUUM_HYBRID = (DEM_SOLIDS .AND. TFM_SOLIDS)
! Set flag for coupled simulations
      DES_CONTINUUM_COUPLED = .NOT.(RO_g0 == 0.0d0)

      IF(DES_CONTINUUM_HYBRID) CALL HYBRID_HACK

! Overwrite user settings if no Lagrangian solids
      IF(.NOT.DISCRETE_ELEMENT) THEN
         DES_CONTINUUM_COUPLED = .FALSE.   ! This keyword might get removed.
         PRINT_DES_DATA = .FALSE.
         DES_ONEWAY_COUPLED = .FALSE.
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_MODEL_PREREQS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: Hybrid_HACK                                             !
!  Purpose: Check the distributed parallel namelist variables.         !
!                                                                      !
!  Author: P. Nicoletti                               Date: 14-DEC-99  !
!  Reviewer: J.Musser                                 Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE HYBRID_HACK

! Number of ranks.
      use run, only: SOLIDS_MODEL

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases.
      use param, only: DIM_M

      use error_manager

      integer :: TFM_MAX
      integer :: DEM_MIN
      integer :: M
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("HYBRID_HACK")

! Initialize the loop variables.
      TFM_MAX = -DIM_M
      DEM_MIN =  DIM_M

! Loop over the phases to see what was specified.
      DO M=1, DIM_M
         SELECT CASE(SOLIDS_MODEL(M))
         CASE ('TFM'); TFM_MAX = max(M, TFM_MAX)
         CASE ('DEM'); DEM_MIN = min(M, DEM_MIN)
         END SELECT
      ENDDO

      write(*,*)'  MAX TFM index: ', TFM_MAX
      write(*,*)'  MIN DEM index: ', DEM_MIN

      IF(DEM_MIN < TFM_MAX) THEN
         WRITE(ERR_MSG, 2000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 2000 FORMAT('Error 2000: Illegal phase specification for hybrid ',    &
         'model. All TFM',/'solids must be defined before DEM solids.',&
         /'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE HYBRID_HACK
