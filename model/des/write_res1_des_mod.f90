      MODULE WRITE_RES1_DES

      use compar, only: PE_IO
      use compar, only: myPE
      use discretelement, only: nonexistent, particle_state
      use out_bin_512_mod, only: out_bin_512
      use out_bin_512i_mod, only: out_bin_512i
      use run, only: bDist_IO

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: INIT_WRITE_RES_DES
      PUBLIC :: FINL_WRITE_RES_DES

      PUBLIC :: WRITE_RES_DES
      PUBLIC :: WRITE_RES_pARRAY

! Write scalar and data WITHOUT MPI data collection.
      INTERFACE WRITE_RES_DES
         MODULE PROCEDURE WRITE_RES_DES_0I, WRITE_RES_DES_1I
         MODULE PROCEDURE WRITE_RES_DES_0D, WRITE_RES_DES_1D
         MODULE PROCEDURE WRITE_RES_DES_0L, WRITE_RES_DES_1L
      END INTERFACE

! Write particle array data.
      INTERFACE WRITE_RES_pARRAY
         MODULE PROCEDURE WRITE_RES_pARRAY_1I
         MODULE PROCEDURE WRITE_RES_pARRAY_1D
         MODULE PROCEDURE WRITE_RES_pARRAY_1L
      END INTERFACE

      INTEGER, PARAMETER :: RDES_UNIT = 901

! Send/Recv parameters for Particle arrays:
      INTEGER :: pROOTCNT, pPROCCNT
      INTEGER :: pSEND

! Send/Recv parameters for Collision/Neighbor arrays:
      INTEGER :: cROOTCNT, cPROCCNT
      INTEGER :: cSEND

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Subroutine: OPEN_RES_DES                                             !
!                                                                      !
! Purpose: Construct the file name and open the DES RES file.          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE OPEN_RES_DES(BASE)

      use machine, only: OPEN_N1

      CHARACTER(len=*), INTENT(IN)  :: BASE
      CHARACTER(len=32) :: lFNAME

      IF(bDIST_IO) THEN
         WRITE(lFNAME,'(A,I4.4,A)') BASE//'_DES_',myPE,'.RES'
         OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',         &
            STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)

      ELSEIF(myPE == PE_IO) THEN
         WRITE(lFNAME,'(A,A)') BASE//'_DES.RES'
         OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',         &
            STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)
      ENDIF

      END SUBROUTINE OPEN_RES_DES

!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_WRITE_RES_DES                                       !
!                                                                      !
! Purpose: Construct the file name and open the DES RES file.          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE INIT_WRITE_RES_DES(BASE, lVERSION, lNEXT_REC)

      use compar, only: numPEs

      use discretelement, only: PIP, iGHOST_CNT
      use discretelement, only: NEIGHBORS, NEIGHBOR_INDEX, NEIGH_NUM

      CHARACTER(len=*), INTENT(IN)  :: BASE
      DOUBLE PRECISION, INTENT(IN) :: lVERSION
      INTEGER, INTENT(OUT) :: lNEXT_REC

! Number of real particles on local rank
      INTEGER :: lPROC
! Total number of real particles.
      INTEGER :: lGHOST_CNT
! Local gather counts for send/recv
      INTEGER :: lGatherCnts(0:NUMPEs-1)
! Loop counters
      INTEGER :: LC1,part

      CALL OPEN_RES_DES(BASE)


      pROOTCNT = PIP
      pPROCCNT = pROOTCNT

      lGHOST_CNT = iGHOST_CNT

      cROOTCNT = NEIGH_NUM
      cPROCCNT = cROOTCNT

! Write out the initial data.
      lNEXT_REC = 1
      CALL WRITE_RES_DES(lNEXT_REC, lVERSION)    ! RES file version
      CALL WRITE_RES_DES(lNEXT_REC, pROOTCNT)    ! Number of Particles
      CALL WRITE_RES_DES(lNEXT_REC, lGHOST_CNT)  ! Number of Ghosts
      CALL WRITE_RES_DES(lNEXT_REC, cROOTCNT)    ! Number of neighbors

      RETURN
      END SUBROUTINE INIT_WRITE_RES_DES

!``````````````````````````````````````````````````````````````````````!
! Subroutine: CLOSE_RES_DES                                            !
!                                                                      !
! Purpose: Close the DES RES file.                                     !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE FINL_WRITE_RES_DES

      IF(bDIST_IO .OR. myPE == PE_IO) close(RDES_UNIT)


      RETURN
      END SUBROUTINE FINL_WRITE_RES_DES

!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_0I                                         !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_0I(lNEXT_REC, INPUT_I)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(IN) :: INPUT_I

      WRITE(RDES_UNIT, REC=lNEXT_REC) INPUT_I
      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE WRITE_RES_DES_0I

!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_1I                                         !
!                                                                      !
! Purpose: Write an array of integers to RES file. Note that data is   !
! not collected and hsould be on rank 0.                               !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_1I(lNEXT_REC, INPUT_I)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(IN) :: INPUT_I(:)

      INTEGER :: lSIZE

      lSIZE = size(INPUT_I)
      CALL OUT_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)

      RETURN
      END SUBROUTINE WRITE_RES_DES_1I

!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_0D                                         !
!                                                                      !
! Purpose: Write scalar double percision values to RES file.           !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_0D(lNEXT_REC, INPUT_D)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(IN) :: INPUT_D

      WRITE(RDES_UNIT, REC=lNEXT_REC) INPUT_D
      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE WRITE_RES_DES_0D

!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_1D                                         !
!                                                                      !
! Purpose: Write an array of double percision values to RES file. Note !
! that data is not collected and should be on rank 0.                  !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_1D(lNEXT_REC, INPUT_D)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(IN) :: INPUT_D(:)

      INTEGER :: lSIZE

      lSIZE = size(INPUT_D)
      CALL OUT_BIN_512(RDES_UNIT, INPUT_D, lSIZE, lNEXT_REC)

      RETURN
      END SUBROUTINE WRITE_RES_DES_1D

!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_0L                                         !
!                                                                      !
! Purpose: Write scalar logical values to RES file.                    !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_0L(lNEXT_REC, INPUT_L)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(IN) :: INPUT_L

      INTEGER :: INPUT_I

      INPUT_I = merge(1,0,INPUT_L)

      WRITE(RDES_UNIT, REC=lNEXT_REC) INPUT_I
      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE WRITE_RES_DES_0L

!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_1L                                         !
!                                                                      !
! Purpose: Write an array of integers to RES file. Note that data is   !
! not collected and hsould be on rank 0.                               !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_1L(lNEXT_REC, INPUT_L)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(IN) :: INPUT_L(:)

      INTEGER, ALLOCATABLE :: INPUT_I(:)

      INTEGER :: lSIZE, LC1

      lSIZE = size(INPUT_L)
      ALLOCATE(INPUT_I(lSIZE))

      DO LC1=1, lSIZE
         INPUT_I(LC1) = merge(1,0,INPUT_L(LC1))
      ENDDO

      CALL OUT_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)

      IF(allocated(INPUT_I)) deallocate(INPUT_I)

      RETURN
      END SUBROUTINE WRITE_RES_DES_1L


!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_PARRAY_1I                                      !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_PARRAY_1I(lNEXT_REC, INPUT_I, pLOC2GLB)

      use desmpi, only: iProcBuf
      use discretelement, only: MAX_PIP, PIP
      use discretelement, only: iGLOBAL_ID

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(IN) :: INPUT_I(:)
      LOGICAL, INTENT(IN), OPTIONAL :: pLOC2GLB

      LOGICAL :: lLOC2GLB
! Loop counters
      INTEGER :: LC1, LC2

      lLOC2GLB = .FALSE.
      IF(present(pLOC2GLB)) lLOC2GLB = pLOC2GLB

      CALL OUT_BIN_512i(RDES_UNIT, INPUT_I(:PIP), PIP, lNEXT_REC)



      RETURN
      END SUBROUTINE WRITE_RES_PARRAY_1I

!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_PARRAY_1D                                      !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_PARRAY_1D(lNEXT_REC, INPUT_D)

      use discretelement, only: MAX_PIP, PIP

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(IN) :: INPUT_D(:)

! Loop counters
      INTEGER :: LC1, LC2

      CALL OUT_BIN_512(RDES_UNIT, INPUT_D(:pip), pip, lNEXT_REC)

      RETURN
      END SUBROUTINE WRITE_RES_PARRAY_1D

!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_PARRAY_1D                                      !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_PARRAY_1L(lNEXT_REC, INPUT_L)

      use desmpi, only: iProcBuf
      use discretelement, only: MAX_PIP, PIP

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(IN) :: INPUT_L(:)

! Loop counters
      INTEGER :: LC1, LC2

      write(6,*)'death in write res parray 1l'
      stop 8832

      RETURN
      END SUBROUTINE WRITE_RES_PARRAY_1L


      END MODULE WRITE_RES1_DES
