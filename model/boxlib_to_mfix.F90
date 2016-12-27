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
   subroutine mfix_MAIN( &
      u_g, v_g, w_g, u_go, v_go, w_go, &
      p_g, p_go, pp_g, ep_g, ep_go, &
      ro_g, ro_go, rop_g, rop_go, &
      rop_ge, rop_gn, rop_gt, &
      d_e, d_n, d_t, &
      tau_u_g ,tau_v_g, tau_w_g,&
      flux_ge, flux_gn, flux_gt, &
      trD_g, lambda_g, mu_g, &
      f_gds, A_m, b_m, drag_bm, flag, &
      particle_state, particle_phase, des_radius, ro_sol, pvol, pmass, &
      omoi, des_pos_new, des_vel_new, des_usr_var, omega_new, des_acc_old,&
      rot_acc_old, drag_fc, fc, tow)  &
      bind(C, name="mfix_MAIN")

      use discretelement, only: max_pip
      use compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

    integer(c_int)         , intent(inout) :: flag&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,4)
    real(c_double), intent(inout) :: A_m&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,-3:3)
    real(c_double), intent(inout) :: b_m&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: ep_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: ep_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: p_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: p_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: ro_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: ro_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: rop_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: rop_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: u_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: u_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: v_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: v_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: w_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: w_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: pp_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: d_e&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: d_t&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: d_n&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: mu_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: lambda_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: trD_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: tau_u_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: tau_v_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: tau_w_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: flux_gE&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: flux_gN&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: flux_gT&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: rop_gE&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: rop_gN&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: rop_gT&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: f_gds&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    real(c_double), intent(inout) :: drag_bm&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1:3)



    integer(c_int), intent(inout) :: particle_state(max_pip)
    integer(c_int), intent(inout) :: particle_phase(max_pip)

    real(c_double), intent(inout) :: des_radius(max_pip)
    real(c_double), intent(inout) :: ro_sol(max_pip)
    real(c_double), intent(inout) :: pvol(max_pip)
    real(c_double), intent(inout) :: pmass(max_pip)
    real(c_double), intent(inout) :: omoi(max_pip)

    real(c_double), intent(inout) :: des_pos_new(max_pip,3)
    real(c_double), intent(inout) :: des_vel_new(max_pip,3)
    real(c_double), intent(inout) :: des_usr_var(max_pip,1)
    real(c_double), intent(inout) :: omega_new(max_pip,3)

    real(c_double), intent(inout) :: des_acc_old(max_pip,3)
    real(c_double), intent(inout) :: rot_acc_old(max_pip,3)
    real(c_double), intent(inout) :: drag_fc(max_pip,3)
    real(c_double), intent(inout) :: fc(max_pip,3)
    real(c_double), intent(inout) :: tow(max_pip,3)

    call mfix(u_g, v_g, w_g, u_go, v_go, w_go, &
              p_g, p_go, pp_g, ep_g, ep_go, &
              ro_g, ro_go, rop_g, rop_go, &
              rop_ge, rop_gn, rop_gt, d_e, d_n, d_t, &
              tau_u_g, tau_v_g, tau_w_g,&
              flux_ge, flux_gn, flux_gt, &
              trd_g, lambda_g, mu_g,  &
              f_gds, A_m, b_m, drag_bm,  &
              flag, &
       particle_state, particle_phase, des_radius,  ro_sol, pvol, pmass, &
       omoi, des_pos_new, des_vel_new, des_usr_var, omega_new, des_acc_old,&
       rot_acc_old, drag_fc, fc, tow)

  end subroutine mfix_MAIN

!**************************************************************************!


!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_get_data(imax_to_c,jmax_to_c,kmax_to_c,dem_solids_to_c) &
       bind(C, name="mfix_get_data")

    use get_data_module, only: get_data
    use geometry, only: imax, jmax, kmax
    use run     , only: dem_solids

    integer(c_int) imax_to_c, jmax_to_c, kmax_to_c, dem_solids_to_c

    call get_data()

    imax_to_c = imax
    jmax_to_c = jmax
    kmax_to_c = kmax

    IF (dem_solids) THEN
       dem_solids_to_c = 1
    ELSE
       dem_solids_to_c = 0
    ENDIF

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
