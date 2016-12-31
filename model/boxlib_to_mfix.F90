module boxlib_to_mfix_module
! _________________________________________________________________

  use compar
  use iso_c_binding
  use bl_fort_module, only : c_real

  implicit none

contains


!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_get_data(imax_to_c, jmax_to_c, kmax_to_c, fluid, &
     dem, steady_state, dt, dt_minC, dt_maxC, tstopC, &
     time, max_nitC, normg, set_normg, &
     cyclic_xC, cyclic_yC, cyclic_zC, cyclic_mf) &
     bind(C, name="mfix_get_data")

    use get_data_module, only: get_data
    use geometry, only: imax, jmax, kmax
    use run, only: dem_solids
    use run, only: dt_min, dt_max, tstop
    use param1, only: is_undefined
    use fld_const, only: ro_g0
    use toleranc, only: norm_g
    use leqsol, only: max_nit
    use geometry, only: cyclic_x,    cyclic_y,    cyclic_z
    use geometry, only: cyclic_x_pd, cyclic_y_pd, cyclic_z_pd
    use geometry, only: cyclic_x_mf, cyclic_y_mf, cyclic_z_mf

    implicit none

    integer, intent(out) :: imax_to_c, jmax_to_c, kmax_to_c
    integer, intent(out) :: fluid
    integer, intent(out) :: dem
    integer, intent(out) :: steady_state
    double precision, intent(out) :: dt_minC, dt_maxC, tstopC
    double precision, intent(out) :: dt, time
    integer, intent(out) :: max_nitC
    double precision, intent(out) :: normg
    integer, intent(out) :: set_normg
    integer, intent(out) :: cyclic_xC, cyclic_yC, cyclic_zC, cyclic_mf

    call get_data(time, dt)

    imax_to_c = imax
    jmax_to_c = jmax
    kmax_to_c = kmax

! Flags for fluid setup
    fluid =  merge(1,0,ro_g0 /= 0.0d0)

    dem   =  merge(1,0,dem_solids)

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

  end subroutine mfix_get_data


!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_set_domain(flag) &
       bind(C, name="mfix_set_domain")

    use set_domain_module, only: set_domain

    use compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

    integer(c_int), intent(inout) :: flag(istart3:iend3,jstart3:jend3,kstart3:kend3,4)

    call set_domain(flag)

  end subroutine mfix_set_domain

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
  subroutine mfix_usr3(u_g, v_g, w_g, p_g) &
       bind(C, name="mfix_usr3")

    use compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

    real(c_double), intent(inout) :: u_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3)
    real(c_double), intent(inout) :: v_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3)
    real(c_double), intent(inout) :: w_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3)
    real(c_double), intent(inout) :: p_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3)

    call usr3(u_g, v_g, w_g, p_g)

  end subroutine mfix_usr3

end module boxlib_to_mfix_module
