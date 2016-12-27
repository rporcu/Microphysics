module output_manager_module

   use iso_c_binding, only: c_double, c_int
   use param1, only: UNDEFINED, UNDEFINED_I, IS_DEFINED, IS_UNDEFINED
   use write_out1_module, only: write_out1

  contains
!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: OUTPUT_MANAGER                                          !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Purpose: Relocate calls to write output files (RES,VTP). This was   !
!  done to simplify the time_march code.                               !
!                                                                      !
!----------------------------------------------------------------------!
     SUBROUTINE OUTPUT_MANAGER(ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g, &
         particle_state, des_radius, ro_sol, des_pos_new, &
         des_vel_new, des_usr_var, omega_new, exit_signal, finished)&
         bind(C, name="mfix_output_manager")

! Global Variables:
!---------------------------------------------------------------------//

      use compar, only: myPE, PE_IO
      use compar  , only:  istart3, iend3, jstart3, jend3, kstart3, kend3

      use machine, only: wall_time
      use output, only: OUT_TIME, OUT_DT
      use output, only: RES_BACKUP_TIME, RES_BACKUP_DT
      use output, only: RES_TIME, RES_DT
      use output, only: USR_TIME, USR_DT
      use output, only: VTP_TIME, VTP_DT
      use param, only: DIMENSION_USR
      use run, only: TIME, DT, TSTOP
      use run, only: dem_solids
      use time_cpu, only: CPU_IO
      use write_des_data_module, only: write_des_data
      use write_res0_des_module, only: write_res0_des
      use write_res1_mod, only: write_res1
      use discretelement, only: max_pip

      IMPLICIT NONE

      real(c_double), INTENT(IN   ) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), INTENT(IN   ) :: p_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), INTENT(IN   ) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), INTENT(IN   ) :: u_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), INTENT(IN   ) :: v_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      real(c_double), INTENT(IN   ) :: w_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      real(c_double), intent(in) :: des_radius(max_pip)
      real(c_double), intent(in) :: ro_sol(max_pip)

      real(c_double), intent(in) :: des_pos_new(max_pip,3)
      real(c_double), intent(in) :: des_vel_new(max_pip,3)
      real(c_double), intent(in) :: des_usr_var(max_pip,1)
      real(c_double), intent(in) :: omega_new(max_pip,3)


      integer(c_int), intent(inout) :: particle_state(max_pip)

! Dummy Arguments:
!---------------------------------------------------------------------//
! Flag that the the user specified batch time (plus buffer) is met.
      integer(c_int), intent(in) :: exit_signal
! Flag that a steady state case is completed.
      integer(c_int), intent(in) :: finished

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter and counter
      INTEGER :: LC, IDX
! Flag that the header (time) has not be written.
      LOGICAL :: HDR_MSG
! Wall time at the start of IO operations.
      DOUBLE PRECISION :: WALL_START
! SPX file extensions.
      CHARACTER(LEN=35) ::  EXT_END
!......................................................................!

! Initialize the file extension array.
      EXT_END = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

! Initial the header flag.
      HDR_MSG = .TRUE.

! Get the current time before any IO operations begin
      WALL_START = WALL_TIME()

! Create a backup copy of the RES file.
      IF(TIME+0.1d0*DT>=RES_BACKUP_TIME) THEN
         RES_BACKUP_TIME = NEXT_TIME(RES_BACKUP_DT)
         CALL BACKUP_RES
      ENDIF

! Write restart file, if needed
      IF(CHECK_TIME(RES_TIME) .OR. EXIT_SIGNAL == 1) THEN

         RES_TIME = NEXT_TIME(RES_DT)
         CALL WRITE_RES1(ep_g, p_g, ro_g, rop_g, u_g, v_g, w_g)
         CALL NOTIFY_USER('.RES;')

         IF(DEM_SOLIDS) THEN
            CALL WRITE_RES0_DES( particle_state, &
               des_radius, ro_sol, des_usr_var, &
               des_pos_new, des_vel_new, omega_new)
            CALL NOTIFY_USER('DES.RES;')
         ENDIF

      ENDIF

! Write standard output, if needed
      IF(CHECK_TIME(OUT_TIME)) THEN
         OUT_TIME = NEXT_TIME(OUT_DT)
         CALL WRITE_OUT1(ep_g,p_g,ro_g)

         CALL NOTIFY_USER('.OUT;')
      ENDIF

! Write special output, if needed
      IDX = 0
      DO LC = 1, DIMENSION_USR
         IF(CHECK_TIME(USR_TIME(LC))) THEN
            USR_TIME(LC) = NEXT_TIME(USR_DT(LC))
            CALL WRITE_USR1 (LC, des_pos_new, des_vel_new, omega_new)
            CALL NOTIFY_USER('.USR:',EXT_END(LC:LC))
            IDX = IDX + 1
         ENDIF
      ENDDO
      IF(IDX /=0) CALL FLUSH_LIST

      IF(DEM_SOLIDS) THEN
         IF(CHECK_TIME(VTP_TIME)) THEN
            VTP_TIME = NEXT_TIME(VTP_DT)
            CALL WRITE_DES_DATA(des_radius, des_pos_new, des_vel_new, des_usr_var)
            CALL NOTIFY_USER('DES.vtp;')
         ENDIF
      ENDIF

      CALL FLUSH_NOTIFY_USER


! Add the amount of time needed for all IO operations to total.
      CPU_IO = CPU_IO + (WALL_TIME() - WALL_START)

      RETURN

      contains


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION CHECK_TIME(lTIME)

      DOUBLE PRECISION, INTENT(IN) :: lTIME

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

      DOUBLE PRECISION, INTENT(IN) :: lWRITE_DT

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
      SUBROUTINE NOTIFY_USER(MSG, EXT)

      use output, only: FULL_LOG
      use funits, only: DMP_LOG
      use funits, only: UNIT_LOG

      CHARACTER(len=*), INTENT(IN) :: MSG
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: EXT


      LOGICAL :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(HDR_MSG) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1000, ADVANCE='NO') TIME
         IF(SCR_LOG) WRITE(*, 1000, ADVANCE='NO') TIME
         HDR_MSG = .FALSE.
      ENDIF

 1000 FORMAT(' ',/' t=',F12.6,' Wrote')

      IF(.NOT.present(EXT)) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1100, ADVANCE='NO') MSG
         IF(SCR_LOG) WRITE(*, 1100, ADVANCE='NO') MSG
      ELSE
         IF(IDX == 0) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1110, ADVANCE='NO') MSG, EXT
            IF(SCR_LOG) WRITE(*, 1110, ADVANCE='NO') MSG, EXT
         ELSE
            IF(DMP_LOG) WRITE(UNIT_LOG, 1120, ADVANCE='NO') EXT
            IF(SCR_LOG) WRITE(*, 1120, ADVANCE='NO') EXT
         ENDIF
      ENDIF

 1100 FORMAT(1X,A)
 1110 FORMAT(1X,A,1x,A)
 1120 FORMAT(',',A)

      RETURN
      END SUBROUTINE NOTIFY_USER

!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE FLUSH_LIST

      use output, only: FULL_LOG
      use funits, only: DMP_LOG
      use funits, only: UNIT_LOG

      LOGICAL :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(DMP_LOG) WRITE(UNIT_LOG,1000, ADVANCE='NO')
      IF(SCR_LOG) WRITE(*,1000, ADVANCE='NO')

 1000 FORMAT(';')

      RETURN
      END SUBROUTINE FLUSH_LIST


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE FLUSH_NOTIFY_USER

      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: DTSOLID
      use error_manager, only: err_msg, flush_err_msg, ival
      use funits, only: DMP_LOG
      use funits, only: UNIT_LOG
      use tunit_module, only: get_tunit
      use machine, only: wall_time
      use output, only: FULL_LOG
      use output, only: NLOG
      use run, only: TIME, NSTEP
      use time_cpu, only: TIME_START
      use time_cpu, only: WALL_START

      DOUBLE PRECISION :: WALL_ELAP, WALL_LEFT, WALL_NOW
      CHARACTER(LEN=9) :: CHAR_ELAP, CHAR_LEFT
      CHARACTER(LEN=4) :: UNIT_ELAP, UNIT_LEFT

      INTEGER :: TNITS
      LOGICAL :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(.NOT.HDR_MSG) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG,1000)
         IF(SCR_LOG) WRITE(*,1000)
      ENDIF

 1000 FORMAT(' ',/' ')

! Write the elapsed time and estimated remaining time
      IF(MOD(NSTEP,NLOG) == 0) THEN

         IF(DEM_SOLIDS .AND. .NOT.DES_CONTINUUM_COUPLED) THEN
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



      RETURN
      END SUBROUTINE FLUSH_NOTIFY_USER

      END SUBROUTINE OUTPUT_MANAGER


!----------------------------------------------------------------------!
! Subroutine: INIT_OUTPUT_VARS                                         !
! Purpose: Initialize variables used for controling ouputs of the      !
! various files.                                                       !
!----------------------------------------------------------------------!
      SUBROUTINE INIT_OUTPUT_VARS

      use machine, only: wall_time
      use output, only: OUT_TIME, OUT_DT
      use output, only: RES_TIME, RES_DT
      use output, only: USR_TIME, USR_DT
      use output, only: VTP_TIME, VTP_DT
      use output, only: RES_BACKUP_TIME, RES_BACKUP_DT
      use output, only: RES_BACKUPS
      use param, only: DIMENSION_USR
      use run, only: RUN_TYPE
      use run, only: TIME, DT
      use time_cpu, only: CPU_IO
      use time_cpu, only: TIME_START
      use time_cpu, only: WALL_START
      use discretelement, only: PRINT_DES_DATA

      use funits, only: CREATE_DIR

      IMPLICIT NONE

! Loop counter
      INTEGER :: LC

! Initialize times for writing outputs
      OUT_TIME = merge(TIME, UNDEFINED, IS_DEFINED(OUT_DT))

! Initialize the amount of time spent on IO
      CPU_IO = 0.0d0

! Initizle RES
      IF (RUN_TYPE == 'NEW') THEN
         RES_TIME = TIME
      ELSE
         IF (IS_DEFINED(DT)) THEN
            RES_TIME = RES_DT *                                        &
               (INT((TIME + 0.1d0*DT)/RES_DT) + 1)
         ENDIF
      ENDIF

! Initizle RES_BACKUP_TIME
      RES_BACKUP_TIME = UNDEFINED
      IF(IS_DEFINED(RES_BACKUP_DT)) RES_BACKUP_TIME =                 &
         RES_BACKUP_DT * (INT((TIME+0.1d0*DT)/RES_BACKUP_DT)+1)

! Initialize USR_TIME
      DO LC = 1, DIMENSION_USR
         USR_TIME(LC) = UNDEFINED
         IF (IS_DEFINED(USR_DT(LC))) THEN
            IF (RUN_TYPE == 'NEW') THEN
               USR_TIME(LC) = TIME
            ELSE
               USR_TIME(LC) = USR_DT(LC) *                             &
                  (INT((TIME+0.1d0*DT)/USR_DT(LC))+1)
            ENDIF
         ENDIF
      ENDDO

      VTP_TIME = UNDEFINED
      IF(IS_DEFINED(VTP_DT)) THEN
         PRINT_DES_DATA = .TRUE.
         IF (RUN_TYPE == 'NEW'.OR.RUN_TYPE=='RESTART_2') THEN
            VTP_TIME = TIME
         ELSE
            VTP_TIME = VTP_DT*(INT((TIME + 0.1d0*DT)/VTP_DT)+1)
         ENDIF
      ENDIF

! Create a subdir for RES backup files.
      IF(RES_BACKUPS /= UNDEFINED_I) CALL CREATE_DIR('BACKUP_RES')


      WALL_START = WALL_TIME()
      TIME_START = TIME

      RETURN
      END SUBROUTINE INIT_OUTPUT_VARS



!----------------------------------------------------------------------!
! Subroutine: BACKUP_RES                                               !
! Purpose: Shift existing RES file backup files by one index, then     !
! create a copy of the current RES file.                               !
!----------------------------------------------------------------------!
      SUBROUTINE BACKUP_RES

      use compar, only: myPE, PE_IO
      use output, only: RES_BACKUPS
      use run, only: DEM_SOLIDS

      IMPLICIT NONE

      INTERFACE
         FUNCTION rename(input,output) BIND(C,name="rename") RESULT(r)
            USE iso_c_binding, only: c_char, c_int
            CHARACTER(kind=c_char),INTENT(in) :: input(*)
            CHARACTER(kind=c_char),INTENT(in) :: output(*)
            INTEGER(c_int)        :: r
         END FUNCTION rename
      END INTERFACE

      CHARACTER(len=256) :: FNAME0, FNAME1
      CHARACTER :: CHAR

      INTEGER :: LC, I, IERR, IREC
      INTEGER, PARAMETER :: ISRC = 1234567
      INTEGER, PARAMETER :: IDST = 12345678

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
      SUBROUTINE SET_FNAME(pFNAME, pEXT, pINDX)

      use run, only: RUN_NAME

      implicit none

      CHARACTER(LEN=*), INTENT(OUT) :: pFNAME
      CHARACTER(LEN=*), INTENT(IN) ::  pEXT
      INTEGER, INTENT(IN), OPTIONAL :: pINDX

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

      RETURN
      END SUBROUTINE SET_FNAME

      END SUBROUTINE BACKUP_RES
end module output_manager_module
