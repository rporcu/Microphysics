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
   subroutine mfix_MAIN(flag, vol_surr, A_m, b_m, ep_g, ep_go, p_g, p_go, &
      ro_g, ro_go, rop_g, rop_go, u_g, u_go, v_g,v_go, w_g, w_go, &
      pp_g, d_e, d_n, d_t, mu_g, lambda_g, trD_g, tau_u_g ,tau_v_g, tau_w_g, &
      flux_ge, flux_gn, flux_gt, rop_ge, rop_gn, rop_gt, &
      f_gds, drag_am, drag_bm) &
      bind(C, name="mfix_MAIN")

    use compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

    integer         , intent(inout) :: flag&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,4)
    double precision, intent(inout) :: vol_surr&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: A_m&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,-3:3)
    double precision, intent(inout) :: b_m&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: ep_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: ep_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: p_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: p_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: ro_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: ro_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: rop_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: rop_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: u_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: u_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: v_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: v_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: w_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: w_go&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: pp_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: d_e&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: d_t&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: d_n&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: mu_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: lambda_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: trD_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: tau_u_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: tau_v_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: tau_w_g&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: flux_gE&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: flux_gN&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: flux_gT&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: rop_gE&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: rop_gN&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: rop_gT&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: f_gds&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: drag_am&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1)
    double precision, intent(inout) :: drag_bm&
       (istart3:iend3,jstart3:jend3,kstart3:kend3,1:3)

    call mfix(flag, vol_surr, A_m, b_m, ep_g, ep_go, p_g, p_go, &
       ro_g, ro_go, rop_g, rop_go, u_g, u_go, v_g,v_go, w_g, w_go, &
       pp_g, d_e, d_n, d_t, mu_g, lambda_g, trD_g, tau_u_g ,tau_v_g, tau_w_g,&
       flux_ge, flux_gn, flux_gt, rop_ge, rop_gn, rop_gt, &
       f_gds, drag_am, drag_bm)

  end subroutine mfix_MAIN


!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_get_data(imax_to_c,jmax_to_c,kmax_to_c) &
       bind(C, name="mfix_get_data")

    use get_data_module, only: get_data
    use geometry, only: imax, jmax, kmax

    integer imax_to_c, jmax_to_c, kmax_to_c

    call get_data()

    imax_to_c = imax
    jmax_to_c = jmax
    kmax_to_c = kmax

  end subroutine mfix_get_data

! **************************************************************************

  subroutine mfix_set_domain(flag) &
       bind(C, name="mfix_set_domain")

    use set_domain_module, only: set_domain

    use compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

    integer, intent(inout) :: flag(istart3:iend3,jstart3:jend3,kstart3:kend3,4)

    call set_domain(flag)

  end subroutine mfix_set_domain

! **************************************************************************

  subroutine mfix_init_fvars(u,v,w,ep,p,ro,rop,nlen) &
       bind(C, name="mfix_init_fvars")

    integer          :: nlen
    double precision :: u(nlen), v(nlen), w(nlen), ep(nlen), p(nlen), ro(nlen), rop(nlen)

!    call init_fvars(u,v,w,ep,p,ro,rop,nlen)

  end subroutine mfix_init_fvars

! **************************************************************************

  subroutine mfix_set_constprop(ro_g,mu_g,lambda_g,nlen) &
       bind(C, name="mfix_set_constprop")

    integer          :: nlen
    double precision :: ro_g(nlen), mu_g(nlen), lambda_g(nlen)

!    call set_constprop(ro_g,mu_g,lambda_g)

  end subroutine mfix_set_constprop

! **************************************************************************

end module boxlib_to_mfix_module
