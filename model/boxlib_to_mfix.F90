module boxlib_to_mfix_module
! _________________________________________________________________

  use bl_fort_module, only : c_real
  use iso_c_binding , only: c_int

  implicit none

contains


!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_get_data(imax_to_c, jmax_to_c, kmax_to_c, fluid, &
     dem, steady_state, dt, dt_minC, dt_maxC, tstopC, &
     time, max_nitC, normg, set_normg, call_udf, &
     cyclic_xC, cyclic_yC, cyclic_zC, cyclic_mf, &
     xlength_C,ylength_C, zlength_C, coord_C) &
     bind(C, name="mfix_get_data")

    use fld_const, only: ro_g0
    use geometry, only: coordinates
    use geometry, only: cyclic_x,    cyclic_y,    cyclic_z
    use geometry, only: cyclic_x_mf, cyclic_y_mf, cyclic_z_mf
    use geometry, only: cyclic_x_pd, cyclic_y_pd, cyclic_z_pd
    use geometry, only: imax, jmax, kmax
    use geometry, only: xlength, ylength, zlength
    use get_data_module, only: get_data
    use leqsol, only: max_nit
    use param1, only: is_undefined
    use run, only: dem_solids, call_usr
    use run, only: dt_min, dt_max, tstop
    use toleranc, only: norm_g

    implicit none

    integer(c_int), intent(out) :: imax_to_c, jmax_to_c, kmax_to_c
    integer(c_int), intent(out) :: fluid
    integer(c_int), intent(out) :: dem, call_udf
    integer(c_int), intent(out) :: steady_state
    real(c_real), intent(out) :: dt_minC, dt_maxC, tstopC
    real(c_real), intent(out) :: dt, time
    real(c_real), intent(out) :: xlength_C, ylength_C, zlength_C
    integer(c_int)         , intent(out) :: max_nitC, coord_C
    real(c_real), intent(out) :: normg
    integer(c_int), intent(out) :: set_normg
    integer(c_int), intent(out) :: cyclic_xC, cyclic_yC, cyclic_zC, cyclic_mf

    call get_data(time, dt)

    imax_to_c = imax
    jmax_to_c = jmax
    kmax_to_c = kmax

! Flags for fluid setup
    fluid =  merge(1,0,ro_g0 /= 0.0d0)

    dem      =  merge(1,0,dem_solids)
    call_udf =  merge(1,0,call_usr)

    steady_state = merge(1,0,is_undefined(dt))

    dt_minC  = dt_min
    dt_maxC  = dt_max
    tstopC   = tstop
    max_nitC = max_nit

    normg = norm_g
    set_normg = merge(1,0,norm_g /= 1.0d0)

    cyclic_xC = merge(1,0,cyclic_x .or. cyclic_x_pd .or. cyclic_x_mf)
    cyclic_yC = merge(1,0,cyclic_y .or. cyclic_y_pd .or. cyclic_y_mf)
    cyclic_zC = merge(1,0,cyclic_z .or. cyclic_z_pd .or. cyclic_z_mf)

    cyclic_mf = merge(1,0,cyclic_x_mf .or. cyclic_y_mf .or. cyclic_z_mf)

    xlength_C = xlength
    ylength_C = ylength
    zlength_C = zlength

    if (coordinates .eq. 'CARTESIAN') then
       coord_C = 0
    else
       print *,'UNKNOWN COORDINATES'
       stop
    end if

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
  subroutine mfix_usr1() &
       bind(C, name="mfix_usr1")

    call usr1

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
  subroutine mfix_usr3(slo, shi, u_g, v_g, w_g, p_g, dx, dy, dz) &
       bind(C, name="mfix_usr3")

    integer     , intent(in   ) :: slo(3),shi(3)

    real(c_real), intent(inout) :: u_g&
        (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(c_real), intent(inout) :: v_g&
        (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(c_real), intent(inout) :: w_g&
        (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(c_real), intent(inout) :: p_g&
        (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(c_real), intent(in   ) :: dx, dy, dz

    call usr3(slo, shi, u_g, v_g, w_g, p_g, dx, dy, dz)

  end subroutine mfix_usr3

end module boxlib_to_mfix_module
