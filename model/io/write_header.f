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
      USE param, only: dimension_3, dimension_m, dimension_n_g, dimension_n_s
      USE run, only: id_version, run_name

      IMPLICIT NONE

!     Memory required for the run
      DOUBLE PRECISION :: MEMORY
!-----------------------------------------------

      IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '
      IF(DMP_LOG)WRITE (UNIT_LOG, 1005) ID_VERSION, ID_NODE
      IF(DMP_LOG)WRITE (UNIT_LOG,1010)RUN_NAME,ID_HOUR,ID_MINUTE,ID_MONTH,ID_DAY,ID_YEAR

      IF (FULL_LOG .and. myPE.eq.PE_IO) THEN    !//d
         WRITE (*, *) ' '
         WRITE (*, 1005) ID_VERSION, ID_NODE
         WRITE(*,1010)RUN_NAME,ID_HOUR,ID_MINUTE,ID_MONTH,ID_DAY,ID_YEAR
      ENDIF

!   Calculate the memory requirement for the present run

      MEMORY = 9. + (8.*DIMENSION_3/ONEMEG)*(95. + 32.*DIMENSION_M + 4.*&
         DIMENSION_N_G + 4.*DIMENSION_M*DIMENSION_N_S)
      IF(DMP_LOG)WRITE (UNIT_LOG, '(1X,A,F7.2,A)') 'Memory required: ', MEMORY, ' Mb'
      IF (FULL_LOG .and. myPE.eq.PE_IO) THEN     !//d
         WRITE (*, '(1X,A,F7.2,A)') 'Memory required: ', MEMORY, ' Mb'
         WRITE (*, 1015)
      ENDIF

      IF(DMP_LOG)WRITE (UNIT_LOG, 1015)

      RETURN
 1005 FORMAT(1X,'MFIX (',A10,') simulation on computer: ',A20)
 1010 FORMAT(1X,'Run name: ',A20,2X,'Time: ',I2,':',I2.0,20X,'Date: ',I2,'-',I2&
         ,'-',I4)
 1015 FORMAT(72('_'))
      END SUBROUTINE WRITE_HEADER
