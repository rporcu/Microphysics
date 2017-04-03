
      !----------------------------------------------------------------------!
      ! Subroutine: init_output_vars                                         !
      ! Purpose: Initialize variables used for controlling ouputs of the      !
      ! various files.                                                       !
      !----------------------------------------------------------------------!

      subroutine init_output_vars(time, dt) &
        bind(C, name="init_output_vars")

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use param , only: dimension_usr
      use param1, only: undefined, undefined_i, is_defined

      use machine, only: wall_time
      use output, only: OUT_TIME, OUT_DT
      use output, only: RES_TIME, RES_DT
      use output, only: USR_TIME, USR_DT
      use output, only: VTP_TIME, VTP_DT
      use output, only: RES_BACKUP_TIME, RES_BACKUP_DT
      use output, only: RES_BACKUPS
      use run, only: RUN_TYPE
      use time_cpu, only: CPU_IO
      use time_cpu, only: TIME_START
      use time_cpu, only: WALL_START
      use discretelement, only: PRINT_DES_DATA

      use funits, only: create_dir

      implicit none

      real(c_real), intent(in) :: time, dt

      ! Loop counter
      integer :: LC

      ! Initialize times for writing outputs
      OUT_TIME = merge(TIME, undefined, is_defined(OUT_DT))

      ! Initialize the amount of time spent on IO
      CPU_IO = 0.0d0

      ! Initialize RES
      IF (RUN_TYPE == 'NEW') THEN
         RES_TIME = TIME
      ELSE
         IF (is_defined(DT)) THEN
            RES_TIME = RES_DT *                                        &
               (INT((TIME + 0.1d0*DT)/RES_DT) + 1)
         ENDIF
      ENDIF

! Initizle RES_BACKUP_TIME
      RES_BACKUP_TIME = undefined
      IF(is_defined(RES_BACKUP_DT)) RES_BACKUP_TIME =                 &
         RES_BACKUP_DT * (INT((TIME+0.1d0*DT)/RES_BACKUP_DT)+1)

! Initialize USR_TIME
      DO LC = 1, dimension_usr
         USR_TIME(LC) = undefined
         IF (is_defined(USR_DT(LC))) THEN
            IF (RUN_TYPE == 'NEW') THEN
               USR_TIME(LC) = TIME
            ELSE
               USR_TIME(LC) = USR_DT(LC) *                             &
                  (INT((TIME+0.1d0*DT)/USR_DT(LC))+1)
            ENDIF
         ENDIF
      ENDDO

      VTP_TIME = undefined
      IF(is_defined(VTP_DT)) THEN
         PRINT_DES_DATA = .TRUE.
         IF (RUN_TYPE == 'NEW'.OR.RUN_TYPE=='RESTART_2') THEN
            VTP_TIME = TIME
         ELSE
            VTP_TIME = VTP_DT*(INT((TIME + 0.1d0*DT)/VTP_DT)+1)
         ENDIF
      ENDIF

! Create a subdir for RES backup files.
      IF(RES_BACKUPS /= undefined_i) CALL create_dir('BACKUP_RES')


      WALL_START = WALL_TIME()
      TIME_START = TIME

      end subroutine init_output_vars

