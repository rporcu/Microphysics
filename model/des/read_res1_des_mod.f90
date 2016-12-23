      MODULE READ_RES1_DES

      use compar, only: PE_IO
      use compar, only: myPE
      use desmpi, only: dprocbuf, iprocbuf, drootbuf, irootbuf, dpar_pos, idispls, iscr_recvcnt, iscattercnts
      use discretelement, only: entering_ghost, particle_state, exiting_ghost, normal_particle, normal_ghost
      use error_manager, only: err_msg, init_err_msg, flush_err_msg, finl_err_msg, ival
      use exit_mod, only: mfix_exit
      use in_binary_512, only: in_bin_512
      use in_binary_512i, only: in_bin_512i
      use run, only: bDist_IO

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: INIT_READ_RES_DES
      PUBLIC :: FINL_READ_RES_DES

      PUBLIC :: READ_PAR_POS
      PUBLIC :: READ_RES_DES
      PUBLIC :: READ_RES_pARRAY
      PUBLIC :: READ_RES_cARRAY

      INTERFACE READ_RES_DES
         MODULE PROCEDURE READ_RES_DES_0I
         MODULE PROCEDURE READ_RES_DES_1I
         MODULE PROCEDURE READ_RES_DES_0D
         MODULE PROCEDURE READ_RES_DES_1D
         MODULE PROCEDURE READ_RES_DES_0L
         MODULE PROCEDURE READ_RES_DES_1L
      END INTERFACE

      INTERFACE READ_RES_pARRAY
         MODULE PROCEDURE READ_RES_pARRAY_1I
         MODULE PROCEDURE READ_RES_pARRAY_1D
         MODULE PROCEDURE READ_RES_pARRAY_1L
      END INTERFACE

      INTERFACE READ_RES_cARRAY
         MODULE PROCEDURE READ_RES_cARRAY_1I
         MODULE PROCEDURE READ_RES_cARRAY_1D
         MODULE PROCEDURE READ_RES_cARRAY_1L
      END INTERFACE


      INTEGER, PARAMETER :: RDES_UNIT = 901

      INTEGER :: pIN_COUNT
      INTEGER :: cIN_COUNT

! Send/Recv parameters for Particle arrays:
      INTEGER :: pROOTCNT, pPROCCNT
      INTEGER :: pRECV
      INTEGER, allocatable :: pSCATTER(:)
      INTEGER, allocatable :: pDISPLS(:)

! Variables used for reading restart file
      INTEGER, ALLOCATABLE :: pRestartMap(:)
      INTEGER, ALLOCATABLE :: cRestartMap(:)

! Send/Recv parameters for Particle arrays:
      INTEGER :: cROOTCNT, cPROCCNT
      INTEGER :: cRECV
      INTEGER, allocatable :: cSCATTER(:)
      INTEGER, allocatable :: cDISPLS(:)

      INTEGER, ALLOCATABLE :: iPAR_COL(:,:)

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_READ_RES_DES                                        !
!                                                                      !
! Purpose: Construct the file name and open the DES RES file.          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE INIT_READ_RES_DES(BASE, lVERSION, lNEXT_REC)

      use discretelement, only: MAX_PIP, PIP
      use discretelement, only: iGHOST_CNT
      use discretelement, only: NEIGH_MAX,NEIGH_NUM

      use compar, only: numPEs
      use machine, only: OPEN_N1


      implicit none

      CHARACTER(len=*), INTENT(IN)  :: BASE
      DOUBLE PRECISION, INTENT(OUT) :: lVERSION
      INTEGER, INTENT(OUT) :: lNEXT_REC

      CHARACTER(len=32) :: lFNAME

! Integer Error Flag
      INTEGER :: IER


      allocate(pSCATTER(0:numPEs-1))
      allocate(pDISPLS(0:numPEs-1))

      allocate(cSCATTER(0:numPEs-1))
      allocate(cDISPLS(0:numPEs-1))


      WRITE(lFNAME,'(A,A)') BASE//'_DES.RES'
      OPEN(UNIT=RDES_UNIT, FILE=lFNAME,     &
         FORM='UNFORMATTED', STATUS='UNKNOWN', ACCESS='DIRECT',  &
         RECL=OPEN_N1)

      READ(RDES_UNIT, REC=1) pIN_COUNT

      READ(RDES_UNIT, REC=1) lVERSION
      READ(RDES_UNIT, REC=2) pIN_COUNT
!     READ(RDES_UNIT, REC=3) -NOTHING-
      READ(RDES_UNIT, REC=4) cIN_COUNT


      IER = 0

! Allocate the particle restart map. This is used in determining were
! particle data is sent. Only process zero needs this array.
      allocate( pRestartMap(pIN_COUNT), STAT=IER)
      IF(IER/=0) THEN
         WRITE(ERR_MSG, 1200) 'pRestartMap', trim(iVAL(pIN_COUNT))
         CALL FLUSH_ERR_MSG
      ENDIF

! Allocate the collision restart map array. All ranks allocatet this
! array so that mapping the collision data can be done in parallel.
      allocate( cRestartMap(cIN_COUNT), STAT=IER)
      IF(IER/=0) THEN
         WRITE(ERR_MSG, 1200) 'cRestartMap', trim(iVAL(cIN_COUNT))
         CALL FLUSH_ERR_MSG
      ENDIF

 1200 FORMAT('Error 1200: Unable to allocate sufficient memory to ',&
         'read in DES',/'restart file. size(',A,') = ',A)

! CALL GLOBAL_ALL_SUM(IER)
      IF(IER/=0) CALL MFIX_EXIT(myPE)


      lNEXT_REC = 5

      RETURN
      END SUBROUTINE INIT_READ_RES_DES


!``````````````````````````````````````````````````````````````````````!
! Subroutine: CLOSE_RES_DES                                            !
!                                                                      !
! Purpose: Close the DES RES file.                                     !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE FINL_READ_RES_DES


      close(RDES_UNIT)

      IF(allocated(dPROCBUF)) deallocate(dPROCBUF)
      IF(allocated(dROOTBUF)) deallocate(dROOTBUF)
      IF(allocated(iPROCBUF)) deallocate(iPROCBUF)
      IF(allocated(iROOTBUF)) deallocate(iROOTBUF)

      IF(allocated(pRestartMap)) deallocate(pRestartMap)
      IF(allocated(cRestartMap)) deallocate(cRestartMap)

      IF(allocated(pSCATTER)) deallocate(pSCATTER)
      IF(allocated(pDISPLS)) deallocate(pDISPLS)

      IF(allocated(cSCATTER)) deallocate(cSCATTER)
      IF(allocated(cDISPLS)) deallocate(cDISPLS)




      RETURN
      END SUBROUTINE FINL_READ_RES_DES



!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_PAR_POS                                             !
!                                                                      !
! Purpose: Generates the mapping used by the scatter routines to send  !
! read data to the correct rank.                                       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_PAR_POS(lNEXT_REC)

      use discretelement, only: PIP
!      use discretelement, only: DES_POS_NEW
      use compar, only: numPEs

      implicit none

      INTEGER, INTENT(INOUT) :: lNEXT_REC

      INTEGER :: lDIMN
      INTEGER :: LC1, lPROC
      INTEGER :: lScatterCNTS(0:NUMPEs-1)
! The number of particles on each process.
      INTEGER :: PAR_CNT(0:NUMPEs-1)

!-----------------------------------------------

! All process read positions for distributed IO restarts.
      DO LC1 = 1, lDIMN
!         CALL READ_RES_DES(lNEXT_REC, DES_POS_NEW(:,LC1))
      ENDDO
      RETURN
      END SUBROUTINE READ_PAR_POS




!``````````````````````````````````````````````````````````````````````!
! Subroutine: GLOBAL_TO_LOC_COL                                        !
!                                                                      !
! Purpose: Generates the mapping used by the scatter routines to send  !
! read data to the correct rank.                                       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE GLOBAL_TO_LOC_COL(iglobal_id)

      use discretelement, only: PIP


      use funits, only: DMP_LOG

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      implicit none

      INTEGER, DIMENSION(:), INTENT(OUT) :: iglobal_id

! Loop counters.
      INTEGER :: LC1, LC2, LC3, IER
      INTEGER :: UNMATCHED
      INTEGER, ALLOCATABLE :: iLOCAL_ID(:)

! Max global id.
      INTEGER :: MAX_ID, lSTAT
! Debug flags.
      LOGICAL :: dFlag
      LOGICAL, parameter :: setDBG = .FALSE.

      CALL INIT_ERR_MSG("GLOBAL_TO_LOC_COL")

! Initialize the error flag.
      IER = 0

! Set the local debug flag.
      dFlag = (DMP_LOG .AND. setDBG)

      MAX_ID = maxval(IGLOBAL_ID(1:PIP))
      ! CALL GLOBAL_ALL_MAX(MAX_ID)

      allocate(iLOCAL_ID(MAX_ID), STAT=lSTAT)
      ! CALL GLOBAL_ALL_SUM(lSTAT)

! All ranks successfully allocated the array. This permits a crude
! but much faster collision owner detection.
      IF(lSTAT /= 0) THEN
         WRITE(ERR_MSG,1000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1000 FORMAT('Error 1000: Unable to allocate sufficient memory to ',&
         'generate the',/'map from global to local particle IDs.')

      iLOCAL_ID = 0
      DO LC1=1, PIP
         iLOCAL_ID(iGLOBAL_ID(LC1)) = LC1
      ENDDO

! Store the particle data.
      LC3 = 1
      LC2 = 0
      UNMATCHED = 0

! FIXME Fix Restart
!      LP1: DO LC1 = 1, cIN_COUNT
!          IF(cRestartMap(LC1) == myPE) THEN
!             LC2 = LC2 + 1
!             NEIGHBORS(LC2) = iLOCAL_ID(iPAR_COL(1,LC1))
!             NEIGHBORS(2,LC2) = iLOCAL_ID(iPAR_COL(2,LC1))
! ! Verify that the local indices are valid. If they do not match it is
! ! likely because one of the neighbor was removed via an outlet at the time
! ! the RES file was written but the ghost data wasn't updated.
!             IF(NEIGHBORS(1,LC2) == 0 .OR. NEIGHBORS(2,LC2) == 0) THEN
!                UNMATCHED = UNMATCHED + 1
!                IF(dFLAG) THEN
!                   WRITE(ERR_MSG,1100) iPAR_COL(1,LC1), NEIGHBORS(1,LC2),   &
!                      iPAR_COL(2,LC1), NEIGHBORS(2,LC2)
!                   CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
!                ENDIF
!                DO WHILE(PEA(LC3,1))
!                   LC3 = LC3 + 1
!                ENDDO
!                NEIGHBORS(2,LC2) = LC3
!             ENDIF
!          ENDIF
!       ENDDO LP1

! 1100 FORMAT('Error 1100: Particle neighbor local indices are invalid.',/  &
!         5x,'Global-ID    Local-ID',/' 1:  ',2(3x,I9),/' 2:  ',2(3x,I9))

      ! CALL GLOBAL_ALL_SUM(UNMATCHED)
      IF(UNMATCHED /= 0) THEN
         WRITE(ERR_MSG,1101) trim(iVal(UNMATCHED))
         CALL FLUSH_ERR_MSG
      ENDIF

 1101 FORMAT(' Warning: 1101: ',A,' particle neighbor datasets were ',&
         'not matched',/' during restart.')

      IF(allocated(iLOCAL_ID)) deallocate(iLOCAL_ID)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE GLOBAL_TO_LOC_COL



!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_0I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_0I(lNEXT_REC, INPUT_I)


      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(OUT) :: INPUT_I

      IF(bDIST_IO) THEN
         READ(RDES_UNIT, REC=lNEXT_REC) INPUT_I
      ELSE
         IF(myPE == PE_IO) READ(RDES_UNIT, REC=lNEXT_REC) INPUT_I
         ! !CALL BCAST(INPUT_I, PE_IO)
      ENDIF

      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE READ_RES_DES_0I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_1I                                              !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1I(lNEXT_REC, INPUT_I)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(OUT) :: INPUT_I(:)

      INTEGER :: lSIZE

      lSIZE = size(INPUT_I)

      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) &
            CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
         ! !CALL BCAST(INPUT_I, PE_IO)
      ENDIF


      RETURN
      END SUBROUTINE READ_RES_DES_1I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_0D                                          !
!                                                                      !
! Purpose: Write scalar double percision values to RES file.           !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_0D(lNEXT_REC, INPUT_D)


      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: INPUT_D

      IF(bDIST_IO) THEN
         READ(RDES_UNIT, REC=lNEXT_REC) INPUT_D
      ELSE
         IF(myPE == PE_IO) READ(RDES_UNIT, REC=lNEXT_REC) INPUT_D
         ! !CALL BCAST(INPUT_D, PE_IO)
      ENDIF
      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE READ_RES_DES_0D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1D                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1D(lNEXT_REC, INPUT_D)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: INPUT_D(:)

      INTEGER :: lSIZE

      lSIZE = size(INPUT_D)

      IF(bDIST_IO) THEN
         CALL IN_BIN_512(RDES_UNIT, INPUT_D, lSIZE, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) &
            CALL IN_BIN_512(RDES_UNIT, INPUT_D, lSIZE, lNEXT_REC)
         ! !CALL BCAST(INPUT_D, PE_IO)
      ENDIF


      RETURN
      END SUBROUTINE READ_RES_DES_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_0L                                          !
!                                                                      !
! Purpose: Write scalar logical values to RES file.                    !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_0L(lNEXT_REC, OUTPUT_L)


      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: OUTPUT_L

      INTEGER :: OUTPUT_I

      OUTPUT_L = .TRUE.

      IF(bDIST_IO)THEN
         READ(RDES_UNIT, REC=lNEXT_REC) OUTPUT_I
      ELSE
         IF(myPE == PE_IO) READ(RDES_UNIT, REC=lNEXT_REC) OUTPUT_I
         ! !CALL BCAST(OUTPUT_I, PE_IO)
      ENDIF

      IF(OUTPUT_I == 1) OUTPUT_L = .TRUE.
      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE READ_RES_DES_0L


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1L                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1L(lNEXT_REC, INPUT_L)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: INPUT_L(:)

      INTEGER, ALLOCATABLE :: INPUT_I(:)

      INTEGER :: lSIZE, LC1

      lSIZE = size(INPUT_I)
      ALLOCATE( INPUT_I(lSIZE))

      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) &
            CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
         ! !CALL BCAST(INPUT_I, PE_IO)
      ENDIF

      DO LC1=1, LSIZE
         IF(INPUT_I(LC1) == 1) THEN
            INPUT_L(LC1) = .TRUE.
         ELSE
            INPUT_L(LC1) = .FALSE.
         ENDIF
      ENDDO

      IF(allocated(INPUT_I)) deallocate(INPUT_I)

      RETURN
      END SUBROUTINE READ_RES_DES_1L


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_pARRAY_1I(lNEXT_REC, OUTPUT_I)

      use discretelement, only: PIP

      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf

      use compar, only: numPEs

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(OUT) :: OUTPUT_I(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(iPROCBUF(pPROCCNT))
      allocate(iROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iScr_RecvCNT = pRECV
      iScatterCNTS = pSCATTER

      CALL IN_BIN_512i(RDES_UNIT, OUTPUT_I, pIN_COUNT, lNEXT_REC)

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_pARRAY_1I



!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_pARRAY_1D                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_pARRAY_1D(lNEXT_REC, OUTPUT_D)

      use discretelement, only: PIP
      use desmpi, only: dRootBuf
      use desmpi, only: dProcBuf
      use compar, only: numPEs

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: OUTPUT_D(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      DOUBLE PRECISION, ALLOCATABLE :: lBUF_D(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(dPROCBUF(pPROCCNT))
      allocate(dROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iScr_RecvCNT = pRECV
      iScatterCNTS = pSCATTER

      CALL IN_BIN_512(RDES_UNIT, OUTPUT_D, pIN_COUNT, lNEXT_REC)

      deallocate(dPROCBUF)
      deallocate(dROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_pARRAY_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_pARRAY_1L                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_pARRAY_1L(lNEXT_REC, OUTPUT_L)

      use discretelement, only: PIP
      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf
      use compar, only: numPEs

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: OUTPUT_L(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)

      allocate(iPROCBUF(pPROCCNT))
      allocate(iROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iScr_RecvCNT = pRECV
      iScatterCNTS = pSCATTER

      allocate(lBUF_I(pIN_COUNT))
      CALL IN_BIN_512i(RDES_UNIT, lBUF_I, pIN_COUNT, lNEXT_REC)
      DO LC1=1,pIN_COUNT
         IF(lBUF_I(LC1) == 1) THEN
            OUTPUT_L(LC1) = .TRUE.
         ELSE
            OUTPUT_L(LC1) = .FALSE.
         ENDIF
      ENDDO
      deallocate(lBUF_I)

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_pARRAY_1L


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_cARRAY_1I(lNEXT_REC, OUTPUT_I)

      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf
      use compar, only: numPEs
      use discretelement, only: NEIGH_NUM

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(OUT) :: OUTPUT_I(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(iPROCBUF(cPROCCNT))
      allocate(iROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iScr_RecvCNT = cRECV
      iScatterCNTS = cSCATTER

      CALL IN_BIN_512i(RDES_UNIT, OUTPUT_I, cIN_COUNT, lNEXT_REC)

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_cARRAY_1I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_cARRAY_1D                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_cARRAY_1D(lNEXT_REC, OUTPUT_D)

      use compar, only: numPEs
      use discretelement, only: NEIGH_NUM
      use desmpi, only: dRootBuf
      use desmpi, only: dProcBuf

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: OUTPUT_D(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      DOUBLE PRECISION, ALLOCATABLE :: lBUF_D(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(dPROCBUF(cPROCCNT))
      allocate(dROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iScr_RecvCNT = cRECV
      iScatterCNTS = cSCATTER


      CALL IN_BIN_512(RDES_UNIT, OUTPUT_D, cIN_COUNT, lNEXT_REC)

      deallocate(dPROCBUF)
      deallocate(dROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_cARRAY_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_pARRAY_1L                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_cARRAY_1L(lNEXT_REC, OUTPUT_L)

      use compar, only: numPEs
      use discretelement, only: NEIGH_NUM
      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: OUTPUT_L(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)

      allocate(iPROCBUF(cPROCCNT))
      allocate(iROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iScr_RecvCNT = cRECV
      iScatterCNTS = cSCATTER

      allocate(lBUF_I(cIN_COUNT))
      CALL IN_BIN_512i(RDES_UNIT, lBUF_I, cIN_COUNT, lNEXT_REC)
      DO LC1=1,cIN_COUNT
         IF(lBUF_I(LC1) == 1) THEN
            OUTPUT_L(LC1) = .TRUE.
         ELSE
            OUTPUT_L(LC1) = .FALSE.
         ENDIF
      ENDDO
      deallocate(lBUF_I)

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_cARRAY_1L

      END MODULE READ_RES1_DES
