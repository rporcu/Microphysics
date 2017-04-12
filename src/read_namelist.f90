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
      SUBROUTINE READ_NAMELIST(time, dt)

      use bc
      use compar, only: mype, pe_io
      use drag, only: drag_c1, drag_d1
      use constant, only: d_p0, gravity, ro_s0
      use deprecated_or_unknown_module, only: deprecated_or_unknown
      use discretelement, only: des_coll_model, des_en_input, des_en_wall_input, des_et_input, des_et_wall_input
      use discretelement, only: des_etat_fac, des_etat_w_fac, v_poisson, vw_poisson
      use discretelement, only: des_explicitly_coupled, des_intg_method, des_oneway_coupled, e_young, ew_young
      use discretelement, only: kn, kn_w, kt_fac, kt_w_fac, mew, mew_w, particles, des_etat_w_fac, print_des_data
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar
      use exit_mod, only: mfix_exit
      use fld_const, only: mu_g0, mw_avg
      use fld_const, only: ro_g0
      use bc, only: cyclic_x, cyclic_y, cyclic_z
      use bc, only: cyclic_x_pd, cyclic_y_pd, cyclic_z_pd
      use ic, only: ic_ep_g, ic_ep_s, ic_p_g, ic_x_w, ic_type
      use ic, only: ic_u_g, ic_u_s, ic_v_g, ic_v_s, ic_w_g, ic_w_s
      use ic, only: ic_x_e, ic_y_n, ic_y_s, ic_z_b, ic_z_t
      use leqsol, only: do_transpose, leq_it, leq_method
      use leqsol, only: leq_pc, leq_sweep, leq_tol, max_nit, ival
      use output, only: full_log, nlog
      use output, only: usr_dt
      use ps, only: ps_massflow_g
      use ps, only: ps_t_g, ps_u_g, ps_v_g, ps_w_g
      use ps, only: ps_x_e, ps_x_g, ps_y_n, ps_y_s, ps_z_b, ps_z_t, ps_x_w
      use remove_comment_module, only: remove_comment
      use remove_comment_module, only: remove_par_blanks
      use residual, only: group_resid, resid_string
      use run, only: call_usr, description, detect_stall, discretize, tstop
      use run, only: dt_fac, dt_max, dt_min, run_name, run_type, solids_model
      use drag, only: drag_type
      use scales, only: p_ref, p_scale
      use residual, only: norm_g, tol_diverge, tol_resid
      use ur_facs, only: ur_fac
      use usr
      use utilities, only: blank_line, line_too_big, seek_comment
      use utilities, only: make_upper_case, replace_tab
      use param1, only: undefined

      IMPLICIT NONE

! Dummy Arguments:
!------------------------------------------------------------------------//
      real(c_real), intent(  out) :: time, dt

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
! Coefficient of restitution (old symbol)
      real(c_real) :: e
! Indicate whether to do a namelist read on the line
      logical :: READ_FLAG
! Logical to check if file exits.
      logical :: lEXISTS
! Error flag
      logical :: ERROR

      CHARACTER(len=256) :: STRING
      integer :: IOS, II

! Flags restricting what data from the mfix.dat to process
      logical :: READ_LOCKED, READ_FULL

      E = UNDEFINED
      READ_FLAG = .TRUE.
      LINE_NO = 0

      READ_LOCKED = .TRUE.
      READ_FULL = .TRUE.

      time = undefined
      dt = undefined

      ! Open the mfix.dat file. Report errors if the file is not located or
      ! there is difficulties opening it.
      inquire(file='mfix.dat',exist=lEXISTS)
      IF(.NOT.lEXISTS) THEN
         IF(myPE == PE_IO) WRITE(*,1000)
         call mfix_exit(myPE)

 1000 FORMAT(2/,1X,70('*')/' From: READ_NAMELIST',/' Error 1000: ',    &
         'The input data file, mfix.dat, is missing. Aborting.',/1x,   &
         70('*'),2/)

      ELSE
         OPEN(UNIT=UNIT_DAT, FILE='mfix.dat', STATUS='OLD', IOSTAT=IOS)
         IF(IOS /= 0) THEN
            IF(myPE == PE_IO) WRITE (*,1100)
            call mfix_exit(myPE)
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
            call mfix_exit(myPE)
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
         write(*,*) 'string:',line_string, line_len
         CALL SET_KEYWORD(ERROR)
         IF (ERROR) THEN
            CALL DEPRECATED_OR_UNKNOWN(LINE_NO, LINE_STRING(1:LINE_LEN))
         ENDIF
      ENDDO

      CLOSE(UNIT=UNIT_DAT)

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
      include 'numerical_params.inc'
      include 'geometry.inc'
      include 'gas_phase.inc'
      include 'solids_phase.inc'
      include 'initial_conditions.inc'
      include 'boundary_conditions.inc'
      include 'point_sources.inc'
      include 'output_control.inc'
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
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&RUN_CONTROL_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=RUN_CONTROL_LOCKED,  IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&RUN_CONTROL_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=RUN_CONTROL_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Physical parameter keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&PHYSICAL_PARAM_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=PHYSICAL_PARAM_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&PHYSICAL_PARAM_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=PHYSICAL_PARAM_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Numerical parameter keywords
      STRING=''; STRING = '&NUMERICAL_PARAM_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=NUMERICAL_PARAM_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Geometry and discretization keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&GEOMETRY_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=GEOMETRY_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&GEOMETRY_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=GEOMETRY_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Gas phase keywords
      STRING=''; STRING = '&GAS_PHASE_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=GAS_PHASE_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Solidss phase keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&SOLIDS_PHASE_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=SOLIDS_PHASE_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&SOLIDS_PHASE_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=SOLIDS_PHASE_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Initial condtion keywords
      STRING=''; STRING = '&INITIAL_CONDITIONS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=INITIAL_CONDITIONS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Boundary condition keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&BOUNDARY_CONDITIONS_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=BOUNDARY_CONDITIONS_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&BOUNDARY_CONDITIONS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=BOUNDARY_CONDITIONS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Point source keywords
      STRING=''; STRING = '&POINT_SOURCES_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=POINT_SOURCES_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Output control keywords
      STRING=''; STRING = '&OUTPUT_CONTROL_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=OUTPUT_CONTROL_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! User hook keywords
      STRING=''; STRING = '&useR_HOOKS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=useR_HOOKS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN

! Stop processing keyword inputs if runing POST_MFIX
       IF(.NOT.READ_FULL) RETURN

      IF(READ_LOCKED) THEN

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


      ENDIF

      IF(READ_LOCKED) ERROR = .TRUE.

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
