module amrex_to_mfix_module
! _________________________________________________________________

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int, c_char

  implicit none

contains

!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_get_data(fluid, &
     dem, steady_state, dt, dt_minC, dt_maxC, tstopC, &
     max_nitC, normg, set_normg, call_udf) &
     bind(C, name="mfix_get_data")

    use fld_const, only: ro_g0
    use get_data_module, only: get_data
    use leqsol, only: max_nit
    use param, only: is_undefined
    use run, only: dem_solids, call_usr
    use run, only: dt_min, dt_max, tstop
    use residual, only: norm_g

    implicit none

    real(rt),   intent(in ) :: tstopC
    integer(c_int), intent(out) :: fluid
    integer(c_int), intent(out) :: dem, call_udf
    integer(c_int), intent(out) :: steady_state
    real(rt),   intent(out) :: dt_minC, dt_maxC
    real(rt),   intent(out) :: dt
    integer(c_int), intent(out) :: max_nitC
    real(rt),   intent(out) :: normg
    integer(c_int), intent(out) :: set_normg

    call get_data(dt)

! Flags for fluid setup
    fluid =  merge(1,0,abs(ro_g0) > tiny(0.0d0))

    dem      =  merge(1,0,dem_solids)
    call_udf =  merge(1,0,call_usr)

    steady_state = merge(1,0,is_undefined(dt))

    dt_minC  = dt_min
    dt_maxC  = dt_max
    max_nitC = max_nit

    ! We now set tstop in the Fortran from the C++
    tstop    = tstopC

    normg = norm_g
    set_normg = merge(1,0, abs(norm_g - 1.0d0) > tiny(0.0d0))

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
  subroutine mfix_usr1(ltime) bind(C, name="mfix_usr1")

    real(rt),   intent(in ) :: ltime

    call usr1(ltime)

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
  subroutine mfix_usr3(u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, &
     p_g, slo, shi, dx, dy, dz) bind(C, name="mfix_usr3")

    integer(c_int), intent(in   ) :: ulo(3),uhi(3)
    integer(c_int), intent(in   ) :: vlo(3),vhi(3)
    integer(c_int), intent(in   ) :: wlo(3),whi(3)
    integer(c_int), intent(in   ) :: slo(3),shi(3)

    real(rt), intent(inout) :: u_g&
        (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
    real(rt), intent(inout) :: v_g&
        (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
    real(rt), intent(inout) :: w_g&
        (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
    real(rt), intent(inout) :: p_g&
        (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(rt), intent(in   ) :: dx, dy, dz

    call usr3(u_g, ulo, uhi, v_g, vlo, vhi, w_g, wlo, whi, p_g, slo, shi,&
              dx, dy, dz)

  end subroutine mfix_usr3

!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!

end module amrex_to_mfix_module
