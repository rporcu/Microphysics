!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR0                                             C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Purpose: This routine is called before the time loop starts and is  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine write_usr0() bind(C, name="write_usr0")

  call write_dat_header('POST_GRAN_TEMP.dat')

contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
  subroutine write_dat_header(fname)

    use run, only: description

    implicit none

    character(len=*) :: fname

! logical used for testing is the data file already exists
    logical :: exists
! file unit for heat transfer data
    integer, parameter :: funit = 2030

    inquire(file=fname,exist=exists)
    if (.not.exists) then
       open(unit=funit,file=fname,status='new')
       write(funit, 1000) trim(description)
    else
       open(unit=funit,file=fname,position="append",status='old')
    endif

    write(funit, 1100)

1000 format(2/,12x,A)
1100 format(2/,2x,'t*sqrt(T0)/dp',5x,'T(t)/T0',10x,'MFIX')

    close(funit)
    return
  end subroutine write_dat_header
end subroutine write_usr0
