MODULE read_namelist_module

   use parse_line_module, only: parse_line

   integer, private :: argc = 0
   character(len=80), private :: argv(32)

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: READ_NAMELIST(POST)                                 !
!     Author: P. Nicoletti                            Date: 25-NOV-91  !
!                                                                      !
!     Purpose: Read in the NAMELIST variables                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_NAMELIST(dt_inout)

      use bc
      use drag, only: drag_c1, drag_d1
      use constant, only: gravity
      use deprecated_or_unknown_module, only: deprecated_or_unknown
      use discretelement, only: des_coll_model, des_en_input, des_en_wall_input, des_et_input, des_et_wall_input, particle_types
      use discretelement, only: des_etat_fac, des_etat_w_fac, v_poisson, vw_poisson
      use discretelement, only: des_explicitly_coupled, des_oneway_coupled, e_young, ew_young
      use discretelement, only: kn, kn_w, kt_fac, kt_w_fac, mew, mew_w, des_etat_w_fac
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ival

      use fld_const, only: mu_g0, mw_avg
      use fld_const, only: ro_g0
      use ic, only: ic_ep_g, ic_ep_s, ic_pack_type, ic_p_g, ic_x_w
      use ic, only: ic_u_g, ic_u_s, ic_v_g, ic_v_s, ic_w_g, ic_w_s
      use ic, only: ic_x_e, ic_y_n, ic_y_s, ic_z_b, ic_z_t
      use ic, only: ic_dp_dist, ic_dp_mean, ic_dp_std, ic_dp_min, ic_dp_max
      use ic, only: ic_ro_s_dist, ic_ro_s_mean, ic_ro_s_std, ic_ro_s_min, ic_ro_s_max
      use run, only: full_log, nlog
      use output, only: usr_x_w, usr_x_e, usr_y_n, usr_y_s, usr_z_b, usr_z_t
      use output, only: usr_dt
      use ps, only: ps_massflow_g
      use ps, only: ps_t_g, ps_u_g, ps_v_g, ps_w_g
      use ps, only: ps_x_e, ps_x_g, ps_y_n, ps_y_s, ps_z_b, ps_z_t, ps_x_w
      use run, only: call_usr, description, tstop
      use run, only: dt_fac, dt_max, dt_min, run_name
      use drag, only: drag_type
      use scales, only: p_ref, p_scale
      use usr
      use utilities, only: blank_line, line_too_big, seek_comment
      use utilities, only: make_upper_case, replace_tab
      use param, only: undefined

      use remove_comment_module, only: remove_comment
      use remove_comment_module, only: remove_par_blanks

      implicit none

! Dummy Arguments:
!------------------------------------------------------------------------//
      real(rt), intent(  out) :: dt_inout

! Local Variables:
!------------------------------------------------------------------------//
      integer, parameter :: unit_dat = 51
! LINE_STRING(1:MAXCOL) has valid input data
      integer, PARAMETER :: MAXCOL = 80
! Holds one line in the input file
      CHARACTER(LEN=512) :: LINE_STRING
! Length of noncomment string
      integer :: LINE_LEN
! Line number
      integer :: LINE_NO
! Indicate whether to do a namelist read on the line
      logical :: READ_FLAG
! Logical to check if file exits.
      logical :: lEXISTS
! Error flag
      logical :: ERROR

      CHARACTER(len=256) :: STRING
      integer :: IOS, II

      real(rt) :: dt

! Flags restricting what data from the mfix.dat to process

      READ_FLAG = .TRUE.
      LINE_NO = 0

      dt = undefined

      ! Open the mfix.dat file. Report errors if the file is not located or
      ! there is difficulties opening it.
      inquire(file='mfix.dat',exist=lEXISTS)
      IF(.NOT.lEXISTS) THEN
         WRITE(*,1000)
         stop 20010

 1000 FORMAT(2/,1X,70('*')/' From: READ_NAMELIST',/' Error 1000: ',    &
         'The input data file, mfix.dat, is missing. Aborting.',/1x,   &
         70('*'),2/)

      ELSE
         OPEN(UNIT=UNIT_DAT, FILE='mfix.dat', STATUS='OLD', IOSTAT=IOS)
         IF(IOS /= 0) THEN
            WRITE (*,1100)
            stop 20011
         ENDIF

      ENDIF
! Loop through the mfix.dat file and process the input data.
      READ_LP: DO
         READ (UNIT_DAT,"(A)",IOSTAT=IOS) LINE_STRING
         IF(IOS < 0) EXIT READ_LP

         LINE_NO = LINE_NO + 1

         LINE_LEN = SEEK_COMMENT(LINE_STRING,LEN(LINE_STRING)) - 1
         call remove_comment(LINE_STRING, LINE_LEN+1, LEN(LINE_STRING))

         IF(LINE_LEN <= 0) CYCLE READ_LP           ! comment line
         IF(BLANK_LINE(LINE_STRING)) CYCLE READ_LP ! blank line

         IF(LINE_TOO_BIG(LINE_STRING,LINE_LEN,MAXCOL) > 0) THEN
            write (*, 1100)  trim(iVAL(LINE_NO)), &
                 &  trim(trim(ival(MAXCOL))), LINE_STRING(1:MAXCOL)
            stop 20012
         ENDIF

 1100 FORMAT(//1X,70('*')/1x,'From: READ_NAMELIST',/1x,'Error 1100: ', &
         'Line ',A,' in file mfix.dat is too long. Input lines should', &
         /1x,'not pass column ',A,'.',2/3x,A,2/1x,&
         'Please correct input deck.',/1X,70('*'),2/)

! All subsequent lines are thermochemical data
         IF(LINE_STRING(1:11) == 'THERMO DATA') EXIT READ_LP

         CALL SET_KEYWORD(ERROR)
         IF (ERROR) THEN
! At this point, the keyword was not identified therefore it is
! either deprecated or unknown.
            CALL DEPRECATED_OR_UNKNOWN(LINE_NO, LINE_STRING(1:LINE_LEN))
         ENDIF

      ENDDO READ_LP

      DO II=1, argc
         line_string = argv(ii)
         line_len = len(line_string)
         CALL SET_KEYWORD(ERROR)
         IF (ERROR) THEN
            CALL DEPRECATED_OR_UNKNOWN(LINE_NO, LINE_STRING(1:LINE_LEN))
         ENDIF
      ENDDO

      CLOSE(UNIT=UNIT_DAT)

      dt_inout = dt

      RETURN

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: SET_KEYWORD(ERROR)                                       !
! Author: P. Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
! Purpose: Process LINE_STRING for MFIX keyword data.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_KEYWORD(ERROR)

      IMPLICIT NONE

      logical, intent(OUT) ::ERROR

! External namelist files:
!---------------------------------------------------------------------//
      include 'run_control.inc'
      include 'physical_params.inc'
      include 'geometry.inc'
      include 'gas_phase.inc'
      include 'initial_conditions.inc'
      include 'boundary_conditions.inc'
      include 'point_sources.inc'
      include 'usr_hooks.inc'
      include 'desnamelist.inc'
      include 'usrnlst.inc'

      ERROR = .FALSE.

      CALL MAKE_UPPER_CASE (LINE_STRING, LINE_LEN)
      CALL REPLACE_TAB (LINE_STRING, LINE_LEN)
      call remove_par_blanks(LINE_STRING)

! Complete arithmetic operations and expand line
      CALL PARSE_LINE (LINE_STRING, LINE_LEN, READ_FLAG)

! Write the current line to a scratch file
! and read the scratch file in NAMELIST format
      IF(.NOT.READ_FLAG) RETURN

! Run control keywords
      STRING=''; STRING = '&RUN_CONTROL '//&
           trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=RUN_CONTROL,  IOSTAT=IOS)
      IF(IOS == 0)  RETURN

! Physical parameter keywords
      STRING=''; STRING = '&PHYSICAL_PARAM '//&
           trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=PHYSICAL_PARAM, IOSTAT=IOS)
      IF(IOS == 0)  RETURN

! Geometry and discretization keywords
      STRING=''; STRING = '&GEOMETRY '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=GEOMETRY, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Gas phase keywords
      STRING=''; STRING = '&GAS_PHASE '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=GAS_PHASE, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Initial condtion keywords
      STRING=''; STRING = '&INITIAL_CONDITIONS '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=INITIAL_CONDITIONS, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Boundary condition keywords
      STRING=''; STRING = '&BOUNDARY_CONDITIONS '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=BOUNDARY_CONDITIONS, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Point source keywords
      STRING=''; STRING = '&POINT_SOURCES '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=POINT_SOURCES, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! User hook keywords
      STRING=''; STRING = '&USER_HOOKS '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=USER_HOOKS, IOSTAT=IOS)
      IF(IOS == 0)  RETURN

! Discrete Element model input parameters.
       STRING=''; STRING = '&DES_INPUT_DATA '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
       READ(STRING, NML=DES_INPUT_DATA, IOSTAT=IOS)
       IF(IOS == 0)  RETURN

! User defined input parameters.
       STRING=''; STRING = '&USR_INPUT_DATA '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
       READ(STRING, NML=USR_INPUT_DATA, IOSTAT=IOS)
       IF(IOS == 0)  RETURN

       ERROR = .TRUE.

      RETURN
      END SUBROUTINE SET_KEYWORD

END SUBROUTINE READ_NAMELIST


subroutine add_argument(fname, nlen) &
   bind(c,name='mfix_add_argument')

   use iso_c_binding, only: c_int, c_float, c_char
   implicit none
   integer(c_int), intent(in) :: nlen
   character(kind=c_char), intent(in) :: fname(nlen)
   integer :: lc, bnd

   argc = argc + 1
   if(argc > 32) then
      write(*,*) 'from add_argument:'
      write(*,*) 'too many arguments!!!!'
      stop 986532
   endif

! Copy over the string, character-by-character.
   bnd = min(nlen,80)
   do lc=1, bnd
      argv(argc)(lc:lc) = fname(lc)
   enddo
! Clear out the remaining string.
   do lc=bnd+1,80
      argv(argc)(lc:lc) = ' '
   enddo

end subroutine add_argument


END MODULE read_namelist_module
