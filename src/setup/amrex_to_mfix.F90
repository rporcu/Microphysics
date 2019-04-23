module amrex_to_mfix_module
! _________________________________________________________________

  use amrex_fort_module, only: rt => amrex_real
  use iso_c_binding,     only: c_int, c_char

  implicit none

contains

!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_get_data( &
       fluid, dem, steady_state, dt, dt_minC, dt_maxC, tstopC, call_udf, &
       namelen, mfix_datC) &
     bind(C, name="mfix_get_data")

    use fld_const, only: ro_g0
    use get_data_module, only: get_data
    use param, only: is_undefined
    use run, only: dem_solids, call_usr
    use run, only: dt_min, dt_max, tstop

    use iso_c_binding, only: C_CHAR, c_null_char

    implicit none

    integer(c_int), intent(out) :: fluid
    integer(c_int), intent(out) :: dem, call_udf
    integer(c_int), intent(out) :: steady_state
    real(rt),   intent(out) :: dt_minC, dt_maxC
    real(rt),   intent(in ) :: tstopC
    real(rt),   intent(out) :: dt

    integer(c_int), intent(in   ) :: namelen
    character(kind=c_char, len=1), dimension (namelen), intent (in) :: mfix_datC

    character(len=namelen) :: mfix_dat
    integer :: lc

    mfix_dat = " "
    do lc=1,namelen
       if( mfix_datC(lc) == c_null_char) then
          exit
       else
          mfix_dat(lc:lc) = mfix_datC(lc)
       endif
    enddo
    ! write(*,*) "Full file name: >"//trim(mfix_dat)//"<"

    call get_data(mfix_dat, dt)

! Flags for fluid setup
    fluid =  merge(1,0,abs(ro_g0) > tiny(0.0d0))

    dem      =  merge(1,0,dem_solids)
    call_udf =  merge(1,0,call_usr)

    steady_state = merge(1,0,is_undefined(dt))

    dt_minC  = dt_min
    dt_maxC  = dt_max

    ! We now set tstop in the Fortran from the C++
    tstop    = tstopC

  end subroutine mfix_get_data

!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_usr0() &
       bind(C, name="mfix_usr0")

    call usr0

  end subroutine mfix_usr0

!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_usr1(time) &
       bind(C, name="mfix_usr1")
    real(rt), intent(in   ) :: time
     
    call usr1(time)

  end subroutine mfix_usr1

!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_usr2() &
       bind(C, name="mfix_usr2")
     
    call usr2

  end subroutine mfix_usr2

!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_finl_err_msg() &
       bind(C, name="mfix_finl_err_msg")

    use error_manager, only: finl_err_msg

    call finl_err_msg

  end subroutine mfix_finl_err_msg


!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_usr3(vel_g, ulo, uhi, &
     p_g, slo, shi, dx, dy, dz) bind(C, name="mfix_usr3")

    integer(c_int), intent(in   ) :: ulo(3),uhi(3)
    integer(c_int), intent(in   ) :: slo(3),shi(3)

    real(rt), intent(inout) :: vel_g&
        (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)
    real(rt), intent(inout) :: p_g&
        (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(rt), intent(in   ) :: dx, dy, dz

    call usr3(vel_g(ulo(1):,ulo(2):,ulo(3):,1), ulo, uhi, &
              vel_g(ulo(1):,ulo(2):,ulo(3):,2), ulo, uhi, &
              vel_g(ulo(1):,ulo(2):,ulo(3):,3), ulo, uhi, &
              p_g, slo, shi, dx, dy, dz)

  end subroutine mfix_usr3

end module amrex_to_mfix_module
