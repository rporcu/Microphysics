module check_numerics_module

  use bl_fort_module, only : c_real
  use iso_c_binding , only: c_int

  implicit none
  private

  public check_numerics

contains

  !  Subroutine: CHECK_NUMERICS                                          !
  !  Purpose: Check the numerics control namelist section                !
  !                                                                      !
  !  Author: P. Nicoletti                               Date: 27-NOV-91  !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_numerics

    use run,           only: discretize, IFILE_NAME
    use param,         only: dim_eqs
    use error_manager, only: finl_err_msg, err_msg, flush_err_msg, &
                           & init_err_msg, ivar, ival

    integer  :: l

    do l = 1,dim_eqs
       if(discretize(l) > 9 .or. discretize(l) < 0) then
          write(err_msg,2002) trim(ivar('DISCRETIZE',l)),&
               trim(ival(discretize(l))), trim(IFILE_NAME)
          call flush_err_msg(abort=.true.)
       endif
    enddo

2002 format('Error 2002: Invalid option ', A,' = ', A, '.',/  &
         'Please correct the ',A,' file.')

  end subroutine check_numerics

end module check_numerics_module
