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

      use compar, only: mype, pe_io
      use machine, only: id_hour, id_minute, id_day, id_month, id_year
      use output, only: full_log
      use run, only: run_name

      IMPLICIT NONE

!-----------------------------------------------


      RETURN
 1005 FORMAT(1X,'MFIX-Exa simulation: ')
 1010 FORMAT(1X,'Run name: ',A20,2X,'Time: ',I2,':',I2.0,20X,'Date: ',I2,'-',I2&
         ,'-',I4)
 1015 FORMAT(72('_'))
      END SUBROUTINE WRITE_HEADER
END MODULE WRITE_HEADER_MODULE
