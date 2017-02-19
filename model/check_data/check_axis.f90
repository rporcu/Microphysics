MODULE CHECK_AXIS_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_AXIS                                              !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!  Purpose: check geometry data for one axis                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_AXIS(NA, DIMEN, ALENGTH, DA, AXIS, &
         AXIS_INDEX)

      USE param1, only: is_undefined
      USE error_manager, only: finl_err_msg, err_msg, ival, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! number of axis cells (IMAX, JMAX, KMAX)
      INTEGER, INTENT(INOUT) :: NA
! maximum number of cells along axis based on domain decomposition
      INTEGER, INTENT(IN) :: DIMEN
! axis length (XLENGTH, YLENGTH, ZLENGTH)
      real(c_real), INTENT(INOUT) :: ALENGTH
! axis checked ('X','Y','Z')
      CHARACTER, INTENT(IN) :: AXIS
! index associated with AXIS ('I','J','K')
      CHARACTER, INTENT(IN) :: AXIS_INDEX
! cell sizes (DX,DY,DZ);
! use explicit dimension for DA
! DA should be dimensioned DA(DIMEN) rather than DA(0:DIMEN+1) to be
! able to use the logic from previous versions that assumed DA(1)
! as the first element.  An error check has been added to ensure that
! DX, DY and DZ definitions in mfix.dat starts with the zeroth
! element; i.e. DA(1).
      real(c_real), INTENT(OUT) :: DA
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! percent error allowed in axis length checks
      real(c_real), PARAMETER :: PERCENT_ERROR = 1.0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop counter
      INTEGER :: LC
! temporary storage
      real(c_real) :: lSUM, lERR
!-----------------------------------------------


      CALL INIT_ERR_MSG("CHECK_AXIS")

! 1) MAKE SURE AT LEAST TWO OF NA, ALENGTH, DA ARE SPECIFIED
      IF (IS_UNDEFINED(NA) .or. IS_UNDEFINED(ALENGTH)) THEN
         WRITE(ERR_MSG, 1101) AXIS, AXIS, AXIS_INDEX
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1101 FORMAT('Error 1101: Insufficient grid information for ',A1,'-',   &
         'axis. You must',/'specify the following: ',  &
         A1,'LENGTH and ',A1,'MAX',/'Please correct the ',    &
         'mfix.dat file.')


      IF(NA>=0 .AND. NA<=DIMEN) THEN

! 4) CELL SIZE NOT SPECIFIED - calculate NON_VARIABLE DA based on
!    input that was specified
            DA = ALENGTH/DBLE(NA)
      ENDIF

! 6) CHECK CONSISTENCY OF AXIS INPUT

! This must be a legacy check because the code shouldn't get here
! without exiting and DIMEN is calculated, not a hard-coded param.
      IF (NA<0 .OR. NA>DIMEN-2) THEN
         WRITE(ERR_MSG, 1001) AXIS_INDEX//'MAX', trim(iVal(NA))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(ALENGTH <= 0.0) THEN
         WRITE(ERR_MSG, 1001) AXIS//'LENGTH'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      lSUM = 0.0
      DO LC = 1, NA
         IF (DA<=0.0 .OR. IS_UNDEFINED(DA)) THEN
            WRITE(ERR_MSG, 1201) trim(iVar(AXIS,LC))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         lSUM = lSUM + DA
      ENDDO

 1201 FORMAT('Error 1201: D',A,' is not specified or negative. ',      &
         'Please correct',/'the mfix.dat file.')


      lERR = 100.0*ABS(lSUM - ALENGTH)/ALENGTH
      IF(lERR > PERCENT_ERROR) THEN
         WRITE(ERR_MSG,1202) AXIS, AXIS, AXIS, ALENGTH, AXIS, lSUM, &
            lERR, PERCENT_ERROR
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

         DO LC = 1, NA
            WRITE(ERR_MSG,"(4x,A,' = ',A)") trim(iVar('D'//AXIS,LC)),   &
               trim(iVal(DA))
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDDO
         WRITE(ERR_MSG,"('Please correct the mfix.dat file')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)
      ENDIF

 1202 FORMAT('Error 1202: ',A1,'LENGTH and sum(D',A1,') are not ',     &
         'consistent.',/3x,A1,'LENGTH = ',g12.5,3x,'sum(D',A1,') = ',  &
         g12.5,/3x,'ERROR   = ',g12.5,3x,'ERR TOL = ',g12.5,/'  ')


      CALL FINL_ERR_MSG

      RETURN


 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_AXIS
END MODULE CHECK_AXIS_MODULE
