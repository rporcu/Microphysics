module boxlib_to_mfix_module
! _________________________________________________________________

  use iso_c_binding
  use bl_fort_module, only : c_real

  implicit none

contains

!**************************************************************************!
!                                                                          !
!                                                                          !
!**************************************************************************!
  subroutine mfix_MAIN() &
       bind(C, name="mfix_MAIN")

     call mfix()

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

  subroutine mfix_set_domain(icbc_flag,flag,flag_e,flag_n,flag_t,nlen,dx_from_c) &
       bind(C, name="mfix_set_domain")

    use geometry, only: AXY, AYZ, AXZ, VOL

    integer          :: nlen
    integer          :: icbc_flag(nlen)
    integer          :: flag(nlen)
    integer          :: flag_e(nlen)
    integer          :: flag_n(nlen)
    integer          :: flag_t(nlen)
    double precision :: dx_from_c(3)

    AXY = dx_from_c(1)*dx_from_c(2)
    AYZ = dx_from_c(2)*dx_from_c(3)
    AXZ = dx_from_c(1)*dx_from_c(3)

    VOL = dx_from_c(1)*dx_from_c(2)*dx_from_c(3)

    print *,'DX ',dx_from_c(1),dx_from_c(2),dx_from_c(3)

!    call set_domain(icbc_flag,flag,flag_e,flag_n,flag_t,nlen)

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
