MODULE WRITE_HEADER_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_HEADER                                           C
!  Purpose: read and verify input data, open files                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-APR-97  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE WRITE_HEADER

      USE compar, only: mype, pe_io
      USE funits, only: unit_log, dmp_log
      USE machine, only: id_hour, id_minute, id_day, id_month, id_year, id_node
      USE output, only: full_log, onemeg
      USE param, only: dimension_3, dimension_n_g, dimension_n_s
      USE run, only: id_version, run_name
      use bl_fort_module, only : c_real

      IMPLICIT NONE

!     Memory required for the run
      real(c_real) :: MEMORY
!-----------------------------------------------

      IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '
      IF(DMP_LOG)WRITE (UNIT_LOG, 1005)
      IF(DMP_LOG)WRITE (UNIT_LOG,1010)RUN_NAME,ID_HOUR,ID_MINUTE,ID_MONTH,ID_DAY,ID_YEAR

      IF (FULL_LOG .and. myPE.eq.PE_IO) THEN    !//d
         WRITE (*, *) ' '
         WRITE (*, 1005)
         WRITE(*,1010)RUN_NAME,ID_HOUR,ID_MINUTE,ID_MONTH,ID_DAY,ID_YEAR
      ENDIF


      IF(DMP_LOG)WRITE (UNIT_LOG, 1015)

      RETURN
 1005 FORMAT(1X,'MFIX-Exa simulation: ')
 1010 FORMAT(1X,'Run name: ',A20,2X,'Time: ',I2,':',I2.0,20X,'Date: ',I2,'-',I2&
         ,'-',I4)
 1015 FORMAT(72('_'))
      END SUBROUTINE WRITE_HEADER
END MODULE WRITE_HEADER_MODULE
