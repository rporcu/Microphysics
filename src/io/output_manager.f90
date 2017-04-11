module output_manager_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   use param1, only: IS_DEFINED, IS_UNDEFINED

  contains
!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: output_manager                                          !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Purpose: Relocate calls to write output files (RES,VTP). This was   !
!  done to simplify the time_march code.                               !
!                                                                      !
!----------------------------------------------------------------------!
     subroutine output_manager(max_pip, time, dt, &
        xlength, ylength, zlength, nstep, &
        particle_state, &
        des_radius, des_pos_new, des_vel_new, &
        omega_new, finished) &
        bind(C, name="mfix_output_manager")

! Global Variables:
!---------------------------------------------------------------------//

      use compar, only: myPE, PE_IO

      use machine, only: wall_time
      use output, only: USR_TIME, USR_DT
      use output, only: VTP_TIME, VTP_DT
      use param, only: DIMENSION_USR
      use run, only: tstop
      use run, only: dem_solids
      use time_cpu, only: CPU_IO
      use write_des_data_module, only: write_des_data


      IMPLICIT NONE

      integer(c_int), intent(in   ) :: max_pip
      real(c_real)  , intent(in   ) :: time, dt, xlength, ylength, zlength
      integer(c_int), intent(in   ) :: nstep

      real(c_real), intent(in) :: des_radius(max_pip)

      real(c_real), intent(in) :: des_pos_new(max_pip,3)
      real(c_real), intent(in) :: des_vel_new(max_pip,3)
      real(c_real), intent(in) :: omega_new(max_pip,3)


      integer(c_int), intent(inout) :: particle_state(max_pip)

! Dummy Arguments:
!---------------------------------------------------------------------//
! Flag that a steady state case is completed.
      integer(c_int), intent(in) :: finished

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter and counter
      integer :: LC, IDX
! Flag that the header (time) has not be written.
      logical :: HDR_MSG
! Wall time at the start of IO operations.
      real(c_real) :: WALL_START
! SPX file extensions.
      CHARACTER(LEN=35) ::  EXT_END
!......................................................................!

! Initialize the file extension array.
      EXT_END = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

! Initial the header flag.
      HDR_MSG = .TRUE.

! Get the current time before any IO operations begin
      WALL_START = WALL_TIME()

! Write special output, if needed
      IDX = 0
      DO LC = 1, DIMENSION_USR
         IF(CHECK_TIME(USR_TIME(LC))) THEN
            USR_TIME(LC) = NEXT_TIME(USR_DT(LC))
            call write_usr1 (LC, time, dt, max_pip, des_pos_new,&
                             des_vel_new, omega_new, xlength, ylength, zlength)
            call notify_user('.USR:',EXT_END(LC:LC))
            IDX = IDX + 1
         ENDIF
      ENDDO
      IF(IDX /=0) CALL FLUSH_LIST

      IF(DEM_SOLIDS) THEN
         IF(CHECK_TIME(VTP_TIME)) THEN
            VTP_TIME = NEXT_TIME(VTP_DT)
            CALL WRITE_DES_DATA(max_pip, particle_state, des_radius, &
               des_pos_new, des_vel_new)
            CALL NOTIFY_useR('DES.vtp;')
         ENDIF
      ENDIF

      CALL FLUSH_NOTIFY_useR

      ! Add the amount of time needed for all IO operations to total.
      CPU_IO = CPU_IO + (WALL_TIME() - WALL_START)

      contains

!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      logical FUNCTION CHECK_TIME(lTIME)

      real(c_real), INTENT(IN) :: lTIME

      IF(IS_UNDEFINED(DT)) THEN
         CHECK_TIME = (FINISHED == 1)
      ELSE
         CHECK_TIME = (TIME+0.1d0*DT>=lTIME).OR.(TIME+0.1d0*DT>=TSTOP)
      ENDIF

      RETURN
      END FUNCTION CHECK_TIME


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      DOUBLE PRECISION FUNCTION NEXT_TIME(lWRITE_DT)

      real(c_real), INTENT(IN) :: lWRITE_DT

      IF (IS_DEFINED(DT)) THEN
         NEXT_TIME = (INT((TIME + 0.1d0*DT)/lWRITE_DT)+1)*lWRITE_DT
      ELSE
         NEXT_TIME = lWRITE_DT
      ENDIF

      RETURN
      END FUNCTION NEXT_TIME


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      subroutine NOTIFY_useR(MSG, EXT)

      use output, only: FULL_LOG

      CHARACTER(len=*), INTENT(IN) :: MSG
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: EXT


      logical :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(HDR_MSG) THEN
         IF(SCR_LOG) WRITE(*, 1000, ADVANCE='NO') TIME
         HDR_MSG = .FALSE.
      ENDIF

 1000 FORMAT(' ',/' t=',F12.6,' Wrote')

      IF(.NOT.present(EXT)) THEN
         IF(SCR_LOG) WRITE(*, 1100, ADVANCE='NO') MSG
      ELSE
         IF(IDX == 0) THEN
            IF(SCR_LOG) WRITE(*, 1110, ADVANCE='NO') MSG, EXT
         ELSE
            IF(SCR_LOG) WRITE(*, 1120, ADVANCE='NO') EXT
         ENDIF
      ENDIF

 1100 FORMAT(1X,A)
 1110 FORMAT(1X,A,1x,A)
 1120 FORMAT(',',A)

      RETURN
      END subroutine NOTIFY_useR

!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      subroutine FLUSH_LIST

      use output, only: FULL_LOG

      logical :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(SCR_LOG) WRITE(*,1000, ADVANCE='NO')

 1000 FORMAT(';')

      RETURN
      END subroutine FLUSH_LIST


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      subroutine flush_notify_user

      use discretelement, only: des_continuum_coupled
      use discretelement, only: DTSOLID
      use error_manager, only: err_msg, flush_err_msg, ival
      use tunit_module, only: get_tunit
      use machine, only: wall_time
      use output, only: FULL_LOG
      use output, only: NLOG
      use time_cpu, only: TIME_START
      use time_cpu, only: WALL_START

      real(c_real) :: WALL_ELAP, WALL_LEFT, WALL_NOW
      CHARACTER(LEN=9) :: CHAR_ELAP, CHAR_LEFT
      CHARACTER(LEN=4) :: UNIT_ELAP, UNIT_LEFT

      integer :: TNITS
      logical :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(.NOT.HDR_MSG) THEN
         IF(SCR_LOG) WRITE(*,1000)
      ENDIF

 1000 FORMAT(' ',/' ')

! Write the elapsed time and estimated remaining time
      IF(MOD(NSTEP,NLOG) == 0) THEN

         IF(DEM_SOLIDS .AND. .NOT.des_continuum_coupled) THEN
            TNITs = CEILING(real((TSTOP-TIME)/DTSOLID))
            WRITE(ERR_MSG, 1100) TIME, DTSOLID, trim(iVal(TNITs))
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)
         ENDIF
 1100 FORMAT(/'Time: ',g12.5,3x,'DT: ',g12.5,3x,'DEM NITs: ',A)

         WALL_NOW = WALL_TIME()
! Calculate the elapsed wall time.
         WALL_ELAP = WALL_NOW - WALL_START
         CALL GET_TUNIT(WALL_ELAP, UNIT_ELAP)
         CHAR_ELAP=''; WRITE(CHAR_ELAP,"(F9.2)") WALL_ELAP
         CHAR_ELAP = trim(adjustl(CHAR_ELAP))
! Estimate the remaining wall time.
         WALL_LEFT = (WALL_NOW-WALL_START)*(TSTOP-TIME)/               &
            max(TIME-TIME_START,1.0d-6)
         CALL GET_TUNIT(WALL_LEFT, UNIT_LEFT)

         IF (IS_DEFINED(DT)) THEN
            CHAR_LEFT=''; WRITE(CHAR_LEFT,"(F9.2)") WALL_LEFT
            CHAR_LEFT = trim(adjustl(CHAR_LEFT))
         ELSE
            CHAR_LEFT = '0.0'
            UNIT_LEFT = 's'
         ENDIF

! Notify the user of usage/remaining wall times.
         WRITE(ERR_MSG,2000)                                           &
            'Elapsed:', trim(CHAR_ELAP), trim(UNIT_ELAP),              &
            'Est. Remaining:',trim(CHAR_LEFT), trim(UNIT_LEFT)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 2000 FORMAT('Wall Time - ',2(A,1X,A,A,4X))

      end subroutine flush_notify_user

      end subroutine output_manager

!----------------------------------------------------------------------!
! Subroutine: BACKUP_RES                                               !
! Purpose: Shift existing RES file backup files by one index, then     !
! create a copy of the current RES file.                               !
!----------------------------------------------------------------------!
      subroutine BACKUP_RES

      use compar, only: myPE, PE_IO
      use output, only: RES_BACKUPS
      use run, only: DEM_SOLIDS

      IMPLICIT NONE

      INTERFACE
         FUNCTION rename(input,output) BIND(C,name="rename") RESULT(r)
            use iso_c_binding, only: c_char, c_int
            CHARACTER(kind=c_char),INTENT(in) :: input(*)
            CHARACTER(kind=c_char),INTENT(in) :: output(*)
            integer(c_int)        :: r
         END FUNCTION rename
      END INTERFACE

      CHARACTER(len=256) :: FNAME0, FNAME1
      CHARACTER :: CHAR

      integer :: LC, I, IERR, IREC
      integer, PARAMETER :: ISRC = 1234567
      integer, PARAMETER :: IDST = 12345678

      IF(myPE /= PE_IO) RETURN

! Shift all the existing backups by one.
      DO LC=RES_BACKUPS,2,-1
         CALL SET_FNAME(FNAME0,'.RES', LC-1)
         CALL SET_FNAME(FNAME1,'.RES', LC)
            i = rename(FNAME0, FNAME1)

         IF(DEM_SOLIDS) THEN
            CALL SET_FNAME(FNAME0,'_DES.RES', LC-1)
            CALL SET_FNAME(FNAME1,'_DES.RES', LC)
            i = rename(FNAME0, FNAME1)
         ENDIF
      ENDDO

! Copy RES to RES1
      CALL SET_FNAME(FNAME0, '.RES')
      CALL SET_FNAME(FNAME1, '.RES' ,1)

      OPEN(UNIT=ISRC, FILE=FNAME0, ACCESS='DIRECT', STATUS='OLD', ACTION='READ', IOSTAT=IERR, RECL=1)
      OPEN(UNIT=IDST, FILE=FNAME1, ACCESS='DIRECT', STATUS='REPLACE', ACTION='WRITE', IOSTAT=IERR, RECL=1)
      IREC = 1
      DO
         READ(UNIT=ISRC, REC=IREC, IOSTAT=IERR) CHAR
         IF (IERR.NE.0) EXIT
         WRITE(UNIT=IDST, REC=I) CHAR
         IREC = IREC + 1
      END DO

      IF(DEM_SOLIDS) THEN
         CALL SET_FNAME(FNAME0, '_DES.RES')
         CALL SET_FNAME(FNAME1, '_DES.RES' ,1)

         OPEN(UNIT=ISRC, FILE=FNAME0, ACCESS='DIRECT', STATUS='OLD', ACTION='READ', IOSTAT=IERR, RECL=1)
         OPEN(UNIT=IDST, FILE=FNAME1, ACCESS='DIRECT', STATUS='REPLACE', ACTION='WRITE', IOSTAT=IERR, RECL=1)
         IREC = 1
         DO
            READ(UNIT=ISRC, REC=IREC, IOSTAT=IERR) CHAR
            IF (IERR.NE.0) EXIT
            WRITE(UNIT=IDST, REC=I) CHAR
            IREC = IREC + 1
         END DO

      ENDIF

      RETURN

      contains

!----------------------------------------------------------------------!
! Subroutine: SET_FNAME                                                !
! Purpose: Set the backup RES file name based on pINDX.                !
!----------------------------------------------------------------------!
      subroutine SET_FNAME(pFNAME, pEXT, pINDX)

      use run, only: RUN_NAME

      implicit none

      CHARACTER(LEN=*), INTENT(OUT) :: pFNAME
      CHARACTER(LEN=*), INTENT(IN) ::  pEXT
      integer, INTENT(IN), OPTIONAL :: pINDX

! Set the file format for backup copies
      pFNAME=''
      IF(.NOT.PRESENT(pINDX)) THEN
         WRITE(pFNAME,1000) trim(RUN_NAME),pEXT
      ELSE
         IF(RES_BACKUPS < 10) THEN
            WRITE(pFNAME,1001) trim(RUN_NAME), pEXT, pINDX
         ELSEIF(RES_BACKUPS < 100) THEN
            WRITE(pFNAME,1002) trim(RUN_NAME), pEXT, pINDX
         ELSEIF(RES_BACKUPS < 1000) THEN
            WRITE(pFNAME,1003) trim(RUN_NAME), pEXT, pINDX
         ELSEIF(RES_BACKUPS < 10000) THEN
            WRITE(pFNAME,1004) trim(RUN_NAME), pEXT, pINDX
         ELSEIF(RES_BACKUPS < 10000) THEN
            WRITE(pFNAME,1005) trim(RUN_NAME), pEXT, pINDX
         ELSE
            WRITE(pFNAME,1006) trim(RUN_NAME), pEXT, pINDX
         ENDIF
      ENDIF

 1000 FORMAT(2A)
 1001 FORMAT('BACKUP_RES/',2A,I1.1)
 1002 FORMAT('BACKUP_RES/',2A,I2.2)
 1003 FORMAT('BACKUP_RES/',2A,I3.3)
 1004 FORMAT('BACKUP_RES/',2A,I4.4)
 1005 FORMAT('BACKUP_RES/',2A,I5.5)
 1006 FORMAT('BACKUP_RES/',2A,I6.6)

      end subroutine SET_FNAME

      end subroutine BACKUP_RES

end module output_manager_module
