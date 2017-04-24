module output_manager_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   use param, only: IS_DEFINED, IS_UNDEFINED

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

      use output, only: USR_TIME, USR_DT
      use param, only: DIM_USR
      use run, only: tstop
      use run, only: dem_solids

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
! SPX file extensions.
      CHARACTER(LEN=35) ::  EXT_END
!......................................................................!

! Initialize the file extension array.
      EXT_END = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

! Initial the header flag.
      HDR_MSG = .TRUE.

! Write special output, if needed
      IDX = 0
      DO LC = 1, DIM_USR
         IF(CHECK_TIME(USR_TIME(LC))) THEN
            USR_TIME(LC) = NEXT_TIME(USR_DT(LC))
            call write_usr1 (LC, time, dt, max_pip, des_pos_new,&
                             des_vel_new, omega_new, xlength, ylength, zlength)
            call notify_user('.USR:',EXT_END(LC:LC))
            IDX = IDX + 1
         ENDIF
      ENDDO
      IF(IDX /=0) CALL FLUSH_LIST

      CALL FLUSH_NOTIFY_useR

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

      use run, only: FULL_LOG

      CHARACTER(len=*), INTENT(IN) :: MSG
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: EXT


      logical :: SCR_LOG

      SCR_LOG = (FULL_LOG)

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

      use run, only: FULL_LOG

      logical :: SCR_LOG

      SCR_LOG = (FULL_LOG)

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
      use run, only: FULL_LOG, nlog

      real(c_real) :: WALL_ELAP, WALL_LEFT, WALL_NOW
      CHARACTER(LEN=9) :: CHAR_ELAP, CHAR_LEFT
      CHARACTER(LEN=4) :: UNIT_ELAP, UNIT_LEFT

      integer :: TNITS
      logical :: SCR_LOG

      SCR_LOG = (FULL_LOG)

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

      ENDIF

 2000 FORMAT('Wall Time - ',2(A,1X,A,A,4X))

      end subroutine flush_notify_user

      end subroutine output_manager

end module output_manager_module
