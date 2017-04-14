      !----------------------------------------------------------------------!
      ! Subroutine: init_output_vars                                         !
      ! Purpose: Initialize variables used for controlling ouputs of the      !
      ! various files.                                                       !
      !----------------------------------------------------------------------!

      subroutine init_output_vars(time, dt) &
        bind(C, name="init_output_vars")

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      use param , only: dim_usr
      use param, only: undefined, undefined_i, is_defined

      use machine, only: wall_time
      use output, only: USR_TIME, USR_DT
      use run, only: RUN_TYPE
      use time_cpu, only: CPU_IO
      use time_cpu, only: TIME_START
      use time_cpu, only: WALL_START

      implicit none

      real(c_real), intent(in) :: time, dt

      ! Loop counter
      integer :: LC

      ! Initialize the amount of time spent on IO
      CPU_IO = 0.0d0

! Initialize USR_TIME
      DO LC = 1, dim_usr
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

      WALL_START = WALL_TIME()
      TIME_START = TIME

      end subroutine init_output_vars
