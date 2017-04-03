      MODULE vtp

      use compar, only: numpes, mype, pe_io
      use discretelement, only: s_time, vtp_findex
      use error_manager, only: err_msg, ival, flush_err_msg, init_err_msg, finl_err_msg

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      IMPLICIT NONE

      integer, PRIVATE :: GLOBAL_CNT
      integer, PRIVATE :: LOCAL_CNT

      integer :: DES_UNIT = 2000

! file unit for ParaView *.pvd data
      integer, PARAMETER :: PVD_UNIT = 2050

! formatted file name
      CHARACTER(LEN=64) :: FNAME_VTP

      INTERFACE VTP_WRITE_DATA
         MODULE PROCEDURE VTP_WRITE_DP1
         MODULE PROCEDURE VTP_WRITE_DP2
         MODULE PROCEDURE VTP_WRITE_I1
      END INTERFACE

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_DP1                                            !
!                                                                      !
! Purpose: Collect and write 1D double percision arrays to the VTP     !
! file. This routine is designed to collect the data for parallel and  !
! serial runs. This routine also manages the distribted IO case.       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_DP1(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      real(c_real), INTENT(in) :: DATA(:)

      integer :: LC

      WRITE(DES_UNIT,1000) NAME
      DO LC=1, GLOBAL_CNT
         WRITE(DES_UNIT,1001,ADVANCE="NO") real(data(LC))
      ENDDO
      WRITE(DES_UNIT,1002)

 1000 FORMAT('<DataArray type="Float32" Name="',A,'" format="ascii">')
 1001 FORMAT(ES14.6,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_DP1

!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_DP2                                            !
!                                                                      !
! Purpose: Collect and write 2D double percision arrays to the VTP     !
! file. This routine is designed to collect the data for parallel and  !
! serial runs. This routine also manages the distribted IO case.       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_DP2(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      real(c_real), INTENT(in) :: DATA(:,:)

      CHARACTER(len=16) :: NOC
      integer :: LB, UB
      integer :: LC1, LC2

      LB = LBOUND(DATA,2)
      UB = UBOUND(DATA,2)
      NOC=''; WRITE(NOC,*) (UB-LB)+1

      WRITE(DES_UNIT,1000) NAME, trim(adjustl(NOC))
      DO LC1=1, GLOBAL_CNT
         DO LC2=LB, UB
            WRITE(DES_UNIT,1001,ADVANCE="NO") &
               real(data(LC1,LC2))
         ENDDO
      ENDDO
      WRITE(DES_UNIT,1002)


 1000 FORMAT('<DataArray type="Float32" Name="',A,'" NumberOf',        &
         'Components="',A,'" format="ascii">')
 1001 FORMAT(ES14.6,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_DP2



!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_I1                                             !
!                                                                      !
! Purpose: Collect and write 1D integer arrays to the VTP file. This   !
! routine is designed to collect the data for parallel and serial      !
! runs. This routine also manages the distribted IO case.              !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_I1(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      integer, INTENT(in) :: DATA(:)

      integer :: LC


      IF(myPE == PE_IO) THEN
         WRITE(DES_UNIT,1000) NAME
         DO LC=1, GLOBAL_CNT
            WRITE(DES_UNIT,1001,ADVANCE="NO") data(LC)
         ENDDO
         WRITE(DES_UNIT,1002)
      ENDIF


 1000 FORMAT('<DataArray type="Float32" Name="',A,'" format="ascii">')
 1001 FORMAT(I10,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_I1


!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_ELEMENT                                        !
!                                                                      !
! Purpose: Write a string to the VTP file. It masks the need to check  !
! the logical before flushing.                                         !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_ELEMENT(ELEMENT)

      CHARACTER(len=*), INTENT(in) :: ELEMENT

      IF(myPE == PE_IO) &
         WRITE(DES_UNIT,"(A)") ELEMENT

      RETURN
      END SUBROUTINE VTP_WRITE_ELEMENT



!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_OPEN_FILE                                            !
!                                                                      !
! Purpose: This routine opens the VTP file and calcualtes the offsets  !
! for dmp data collection.                                             !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_OPEN_FILE(max_pip, particle_state, NoPc)

! Modules
!-----------------------------------------------
      use run, only: run_type, run_name
      use discretelement, only: normal_particle

      IMPLICIT NONE

      integer         , intent(in   ) :: max_pip
      integer         , intent(in   ) :: particle_state(max_pip)
      character(len=*), intent(  out) :: nopc

      integer :: NumberOfPoints

! check whether an error occurs in opening a file
      integer :: IOS
! Integer error flag.
      integer :: IER

! logical used for testing is the data file already exists
      logical :: EXISTS_VTP
! status of the vtp file to be written
      CHARACTER(LEN=8) :: STATUS_VTP

      integer :: lc


! Calculate the number of 'real' particles on the local process.
      local_cnt = 0
      do lc=1,max_pip
         if(particle_state(lc) == normal_particle) &
            local_cnt = local_cnt + 1
      enddo

! Calculate the total number of particles system-wide.
      global_cnt = local_cnt
      NumberOfPoints = GLOBAL_CNT
      WRITE(NoPc,"(I10.10)") NumberOfPoints

! set the file name and unit number and open file
      WRITE(fname_vtp,'(A,"_DES_",I5.5,".vtp")') &
         trim(run_name), vtp_findex

      IER = 0
      IF(myPE == PE_IO) THEN

! The file should be new but could exist due to restarting.
         STATUS_VTP = 'NEW'
! Check to see if the file already exists.
         INQUIRE(FILE=FNAME_VTP,EXIST=EXISTS_VTP)
! The given file should not exist if the run type is NEW.
         IF(EXISTS_VTP)THEN
! The VTP should never exist for a NEW run.
            IF(RUN_TYPE == 'NEW')THEN
               IER = 1
! The file may exist during a RESTART.
            ELSE
               STATUS_VTP = 'REPLACE'
            ENDIF
         ENDIF

! Open the file and record any errors.
         IF(IER == 0) THEN
            OPEN(UNIT=DES_UNIT, FILE=FNAME_VTP,                        &
               STATUS=STATUS_VTP, IOSTAT=IOS)
            IF(IOS /= 0) IER = 2
         ENDIF
      ENDIF

      IF(IER /= 0) THEN
         CALL INIT_ERR_MSG("VTP_MOD --> OPEN_VTP")
         WRITE(ERR_MSG, 1100) IER
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Unable to open VTP file. This could be ',    &
         'caused by a VTP',/'file with the same file name already ',   &
         'existing. or an error code',/' returned by the OPEN ',       &
         'function.'/'Error code: ',I2,4x,'Aborting.')

      END SUBROUTINE VTP_OPEN_FILE


!......................................................................!
! SUBROUTINE: VTP_CLOSE_FILE                                           !
!                                                                      !
! Purpose: This routine closes the vtp file.                           !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_CLOSE_FILE


      VTP_FINDEX=VTP_FINDEX+1

      IF(myPE .eq.pe_IO) CLOSE(des_unit)


      END SUBROUTINE VTP_CLOSE_FILE


!......................................................................!
! SUBROUTINE: ADD_VTP_TO_PVD                                           !
!                                                                      !
! Purpose: This routine opens the pvd file.                            !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE ADD_VTP_TO_PVD

      use run, only: RUN_TYPE, RUN_NAME

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Index position of desired character
      integer IDX_f, IDX_b
! logical used for testing is the data file already exists
      logical :: EXISTS_PVD
! Generic input limited to 256 characters
      CHARACTER(LEN=256) INPUT

! formatted file name
      CHARACTER(LEN=64) :: FNAME_PVD = ''

      logical, SAVE :: FIRST_PASS = .TRUE.

! IO Status flag
      integer :: IOS

! Variables related to gather
      integer :: IER

!-----------------------------------------------

      CALL INIT_ERR_MSG('VTP_MOD --> ADD_VTP_TO_PVD')

! Initialize the error flag.
      IER = 0

! Obtain the file name and open the pvd file
      FNAME_PVD = TRIM(RUN_NAME)//'_DES.pvd'

! The PVD file is only written by PE_IO with serial IO.
      IF(myPE == PE_IO) THEN

! Check to see if the file already exists.
         INQUIRE(FILE=FNAME_PVD,EXIST=EXISTS_PVD)

         IF(FIRST_PASS) THEN

! Open the "NEW" file and write the necessary header information.
            IF(RUN_TYPE /= 'RESTART_1')THEN

! The file exists but first_pass is also true so most likely an existing
! file from an earlier/other run is present in the directory. Exit to
! prevent accidently overwriting the existing file.
               IF(EXISTS_PVD) THEN
                  IER = 1
               ELSE
                  OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,STATUS='NEW')
                  WRITE(PVD_UNIT,"(A)")'<?xml version="1.0"?>'
                  WRITE(PVD_UNIT,"(A)")'<VTKFile type="Collection" &
                     &version="0.1" byte_order="LittleEndian">'
                  WRITE(PVD_UNIT,"(3X,'<Collection>')")
               ENDIF

! This is the first pass of a restart run. Extra care is needed to make
! sure that the pvd file is ready to accept new data.
            ELSE ! a restart run
               IF(EXISTS_PVD) THEN
! Open the file at the beginning.
                  OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,&
                     POSITION="REWIND",STATUS='OLD',IOSTAT=IOS)
                  IF(IOS /= 0) IER = 2
               ELSE ! a pvd file does not exist
                  IER = 3
               ENDIF

               IF(IER == 0) THEN
! Loop over the entries in the PVD file, looking for a match to the
! file that is being written. If no match is found, the data will be
! appended to the end of the pvd file, otherwise, the old data will
! be over-written.
                  DO
! Read in the entires of the PVD file.
                     READ(PVD_UNIT,"(A)",IOSTAT=IOS)INPUT
                     IF(IOS > 0) THEN
                        IER = 4
                        EXIT
                     ELSEIF(IOS<0)THEN
! The end of the pvd file has been reached without finding an entry
! matching the current record. Exit the loop.
                        BACKSPACE(PVD_UNIT)
                        EXIT
                     ENDIF
! Find the first instances of file=" and "/> in the read data.
                     IDX_f = INDEX(INPUT,'file="')
                     IDX_b = INDEX(INPUT,'"/>')
! Skip rows that do not contain file data
                     IF(IDX_f == 0 .AND. IDX_b == 0) CYCLE
! Truncate the file name from the read data
                     WRITE (INPUT,"(A)") INPUT(IDX_f+6:IDX_b-1)
! If the file name matches the current VTP record, break the loop to
! over-write this record.
                     IF(TRIM(FNAME_VTP) == TRIM(INPUT)) THEN
                        BACKSPACE(PVD_UNIT)
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF ! No errors
            ENDIF ! run_type new or restart

         ELSE ! not FIRST_PASS
            OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,&
               POSITION="APPEND",STATUS='OLD',IOSTAT=IOS)
            IF (IOS /= 0) IER = 2
         ENDIF

      ENDIF ! if myPE == PE_IO and not distributed IO


      ! CAlL GLOBAL_ALL_SUM(IER)
      IF(IER /= 0) THEN
         SELECT CASE(IER)
         CASE(1); WRITE(ERR_MSG,1101) trim(FNAME_PVD)
         CASE(2); WRITE(ERR_MSG,1102) trim(FNAME_PVD)
         CASE(3); WRITE(ERR_MSG,1103) trim(FNAME_PVD)
         CASE(4); WRITE(ERR_MSG,1104) trim(FNAME_PVD)
         CASE DEFAULT; WRITE(ERR_MSG,1105) trim(FNAME_PVD)
         END SELECT
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1101 FORMAT('Error 1101: A PVD file was detected in the run ',        &
         'directory which should',/'not exist for a NEW run.',/        &
         'File: ',A)

 1102 FORMAT('Error 1102: Fatal error status returned while OPENING ', &
         'PVD file.',/'File: ', A)

 1103 FORMAT('Error 1103: PVD file MISSING from run directory.',/      &
         'File: ',A)

 1104 FORMAT('Error 1104: Fatal error status returned while READING ', &
         'PVD file.',/'File: ', A)

 1105 FORMAT('Error 1105:: Fatal unclassified error when processing ', &
         'PVD file.',/'File: ', A)


! If there were no errors, updated the file.
      IF(myPE == PE_IO) THEN

! Remove the last two lines written so that additional data can be added
         IF(.NOT.FIRST_PASS) THEN
            BACKSPACE(PVD_UNIT)
            BACKSPACE(PVD_UNIT)
         ENDIF

! Write the data to the file
         WRITE(PVD_UNIT,"(6X,A,A,A,A,A,A,A)")&
         '<DataSet timestep="',TRIM(iVal(S_TIME)),'" ',&
         'group="" part="0" ',& ! necessary file data
         'file="',TRIM(FNAME_VTP),'"/>' ! file name of vtp

! Write the closing tags
         WRITE(PVD_UNIT,"(3X,A)")'</Collection>'
         WRITE(PVD_UNIT,"(A)")'</VTKFile>'

         CLOSE(PVD_UNIT)
      ENDIF

! Identify that the files has been created and opened for next pass
      FIRST_PASS = .FALSE.

      CALL FINL_ERR_MSG

! Return to the calling routine
      RETURN

      END SUBROUTINE ADD_VTP_TO_PVD
      END MODULE VTP
