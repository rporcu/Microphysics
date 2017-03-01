module check_axis_module

  use bl_fort_module, only: c_real
  use iso_c_binding , only: c_int
  use run,            only: IFILE_NAME
  use error_manager,  only: finl_err_msg, err_msg, ival, flush_err_msg, &
                          & init_err_msg, ivar

  implicit none
  private

  public check_axis

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: CHECK_AXIS                                              !
  !  Author: P. Nicoletti                               Date: 27-NOV-91  !
  !                                                                      !
  !  Purpose: check geometry data for one axis                           !
  !                                                                      !
  !  Comments:                                                           !
  !    cell sizes (DX,DY,DZ);                                            !
  !    use explicit dimension for DA                                     !
  !    DA should be dimensioned DA(DIMEN) rather than DA(0:DIMEN+1) to be!
  !    able to use the logic from previous versions that assumed DA(1)   !
  !    as the first element.  An error check has been added to ensure    !
  !    that DX, DY and DZ definitions in the input file starts with the  !
  !    zeroth element; i.e. DA(1).                                       !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_axis(NA, DIMEN, ALENGTH, DA, AXIS,AXIS_INDEX)

    use param1, only: is_undefined

    integer,      intent(inout) :: NA
    integer,      intent(in)    :: DIMEN
    real(c_real), intent(inout) :: ALENGTH
    character,    intent(in)    :: AXIS
    character,    intent(in)    :: AXIS_INDEX
    real(c_real), intent(out)   :: DA
    real(c_real), parameter     :: PERCENT_ERROR = 1.0 !percent error allowed in axis length checks
    integer                     :: LC
    real(c_real)                :: lSUM, lERR

    return


    call init_err_msg("CHECK_AXIS")

    ! 1) MAKE SURE AT LEAST TWO OF NA, ALENGTH, DA ARE SPECIFIED
    if (IS_UNDEFINED(NA) .or. IS_UNDEFINED(ALENGTH)) then
       write(ERR_MSG, 1101) AXIS, AXIS, AXIS_INDEX, trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif
    
1101 format('Error 1101: Insufficient grid information for ',A1,'-',   &
          'axis. You must',/'specify the following: ',  &
          A1,'LENGTH and ',A1,'MAX',/'Please correct the ',    &
          A, '  file.')    
    
    if(NA>=0 .and. NA<=DIMEN) then       
       ! 4) CELL SIZE NOT SPECIFIED - calculate NON_VARIABLE DA based on
       !    input that was specified
       DA = ALENGTH/dble(NA)
    endif
    
    ! This must be a legacy check because the code shouldn't get here
    ! without exiting and DIMEN is calculated, not a hard-coded param.
    if (NA<0 .or. NA>DIMEN-2) then
       write(ERR_MSG, 1001) AXIS_INDEX//'MAX', trim(iVal(NA)), trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif
    
    if(ALENGTH <= 0.0) then
       write(ERR_MSG, 1001) AXIS//'LENGTH', trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    endif
    
    lSUM = 0.0
    do LC = 1, NA
       if (DA<=0.0 .or. IS_UNDEFINED(DA)) then
          write(ERR_MSG, 1201) trim(iVar(AXIS,LC)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
       lSUM = lSUM + DA
    enddo
    
1201 format('Error 1201: D',A,' is not specified or negative. ',      &
          'Please correct',/'the ',A,' file.')    
    
    lERR = 100.0*abs(lSUM - ALENGTH)/ALENGTH
    if(lERR > PERCENT_ERROR) then
       write(ERR_MSG,1202) AXIS, AXIS, AXIS, ALENGTH, AXIS, lSUM, &
            lERR, PERCENT_ERROR
       call flush_err_msg(FOOTER=.false.)
       
       do LC = 1, NA
          write(ERR_MSG,"(4x,A,' = ',A)") trim(iVar('D'//AXIS,LC)),   &
               trim(iVal(DA))
          call flush_err_msg(HEADER=.false., FOOTER=.false.)
       enddo
       write(ERR_MSG,"('Please correct ',A,' file')") trim(IFILE_NAME)
       call flush_err_msg(HEADER=.false., ABORT=.true.)
    endif
    
1202 format('Error 1202: ',A1,'LENGTH and sum(D',A1,') are not ',     &
          'consistent.',/3x,A1,'LENGTH = ',g12.5,3x,'sum(D',A1,') = ',  &
          g12.5,/3x,'ERROR   = ',g12.5,3x,'ERR TOL = ',g12.5,/'  ')
    
    
    call finl_err_msg
    
    
1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
          'Please correct the ',A,' file.')
    
  end subroutine check_axis


end module check_axis_module
