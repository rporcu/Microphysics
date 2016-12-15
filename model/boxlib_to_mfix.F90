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
  subroutine mfix_MAIN(flag,vol_surr) &
       bind(C, name="mfix_MAIN")

    use compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

    integer         , intent(inout) :: flag(istart3:iend3,jstart3:jend3,kstart3:kend3,4)
    double precision, intent(inout) :: vol_surr(istart3:iend3,jstart3:jend3,kstart3:kend3,1)

    call mfix(flag,vol_surr)

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
