      MODULE FUNITS

! Whether this processor should write the log file in DMP mode.
! Usually this flag is true only for PE_IO.  All the PEs may be forced
! to write a log file by setting ENABLE_DMP_LOG to .true. in output_mod.f.
      LOGICAL :: DMP_LOG

! Flag indicated that the log was opened globally.
      LOGICAL :: LOG_WAS_CLOSED = .FALSE.

! RRATES debug file unit number
      integer, PARAMETER :: UNIT_RRATES = 43

! mfix.dat file unit number
      integer, PARAMETER :: UNIT_DAT = 51

! RUN_NAME.OUT file unit number
      integer, PARAMETER :: UNIT_OUT = 52

! RUN_NAME.LOG file unit number. (DEFAULT/Serial 53)
      integer, PARAMETER :: UNIT_LOG = 53

! Temporary (scratch) file unit number
      integer, PARAMETER :: UNIT_TMP = 54

! RUN_NAME.RES file unit number
      integer, PARAMETER :: UNIT_RES = 55

      CONTAINS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: NEWUNIT                                                !
!  Author: A. Choudhary                               Date: 01/21/2015 !
!                                                                      !
!  Purpose: Finds an open i/o unit number; Usage:                      !
!   integer myunit                                                     !
!   open(convert='big_endian',unit=newunit(myunit),file='filename')    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      integer FUNCTION newunit(unit)
      IMPLICIT NONE

! optional variable to hold unit number
      integer, INTENT(OUT), OPTIONAL  :: unit
! lower and upper limits to search for available units
      integer, PARAMETER  :: lun_min = 100, lun_max= 999
! check to see if the unit is open
      LOGICAL :: is_open
! looping variable
      integer :: lun

      newunit = -1

      DO lun = lun_min, lun_max
        INQUIRE(UNIT=lun, OPENED=is_open)
        IF(.NOT.is_open) THEN
          newunit = lun
          EXIT
        END IF
      END DO

      IF(present(unit)) unit=newunit

      RETURN
      END FUNCTION NEWUNIT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CREATE_DIR                                             !
!  Author: J.Musser                                      Date: 06/2015 !
!                                                                      !
!  Purpose: Create the directory pDIR in the run directory.            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CREATE_DIR(PDIR)

      use compar, only: myPE, PE_IO
      use iso_c_binding, only: c_int16_t

      IMPLICIT NONE

      interface
         function mkdir(path,mode) bind(c,name="mkdir")
            use iso_c_binding, only: c_int16_t, c_int, c_char
            integer(c_int) :: mkdir
            character(kind=c_char,len=1) :: path(*)
            integer(c_int16_t), value :: mode
         end function mkdir

         function rmdir(path,mode) bind(c,name="rmdir")
            use iso_c_binding, only: c_int16_t, c_int, c_char
            integer(c_int) :: rmdir
            character(kind=c_char,len=1) :: path(*)
            integer(c_int16_t), value :: mode
         end function rmdir
      end interface

      CHARACTER(LEN=*), INTENT(IN) :: pDIR

      CHARACTER(LEN=256) :: CMD
      integer :: IOS, I
      integer, PARAMETER :: tUNIT = 9638

      IF(myPE /= PE_IO) RETURN

      OPEN(FILE=trim(pDIR)//'/tmp',UNIT=tUNIT,STATUS='NEW',IOSTAT=IOS)

      IF(IOS == 0 )THEN
         write(*,*) 'The directory already exists.'
         close(tUNIT)
         WRITE(CMD,"(A,'/tmp')")adjustl(trim(pDIR))
         i = rmdir(CMD, int(o'777',c_int16_t))
      ELSE
         write(*,*) 'Creating the directory.'
         WRITE(CMD,"('mkdir ',A)")pDIR
         i = mkdir(pDIR, int(o'777',c_int16_t))
      ENDIF

      RETURN
      END SUBROUTINE CREATE_DIR


      END MODULE FUNITS
