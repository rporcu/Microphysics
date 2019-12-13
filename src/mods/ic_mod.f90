!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: ic                                                          !
!  Purpose: Global initial conditions variables.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module ic

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int, c_char, c_null_char

  use param, only: dim_ic, dim_m, dim_n_g, dim_n_s

  ! Boundary condition coordinates
  real(rt) :: IC_X_w(dim_ic), IC_X_e(dim_ic)
  real(rt) :: IC_Y_s(dim_ic), IC_Y_n(dim_ic)
  real(rt) :: IC_Z_b(dim_ic), IC_Z_t(dim_ic)

  ! Void fraction in a specified boundary
  real(rt) :: IC_EP_g(dim_ic), IC_EP_s(dim_ic, dim_m)

  ! Initial gas pressure
  real(rt) :: IC_P_g(dim_ic)

  ! Initial velocities in specified region
  real(rt) :: IC_U_g(dim_ic), IC_U_s(dim_ic, dim_m)
  real(rt) :: IC_V_g(dim_ic), IC_V_s(dim_ic, dim_m)
  real(rt) :: IC_W_g(dim_ic), IC_W_s(dim_ic, dim_m)

  ! Heat transfer boundary condition
  real(rt) :: IC_T_g(dim_ic), IC_T_s(dim_ic, dim_m)

  ! Species transfer boundary condition
  real(rt) :: IC_X_g(dim_ic, dim_n_g), IC_X_s(dim_ic, dim_m, dim_n_s)

  ! Particle Size properties
  character(len=16) :: ic_dp_dist(dim_ic, dim_m)
  real(rt) :: ic_dp_mean(dim_ic, dim_m)
  real(rt) :: ic_dp_std(dim_ic, dim_m)
  real(rt) :: ic_dp_min(dim_ic, dim_m)
  real(rt) :: ic_dp_max(dim_ic, dim_m)

  ! Particle density properties
  character(len=16) :: ic_ro_s_dist(dim_ic, dim_m)
  real(rt) :: ic_ro_s_mean(dim_ic, dim_m)
  real(rt) :: ic_ro_s_std(dim_ic, dim_m)
  real(rt) :: ic_ro_s_min(dim_ic, dim_m)
  real(rt) :: ic_ro_s_max(dim_ic, dim_m)

  character(len=16) :: ic_pack_type(dim_ic)

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutines:                                                         !
!                                                                      !
! Purpose: Getters for the initial conditions values                   !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  real(rt) function get_ic_p_g (pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_ic_p_g = ic_p_g(pID)
    return
  end function get_ic_p_g

  real(rt) function get_ic_ep_g (pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_ic_ep_g = ic_ep_g(pID)
    return
  end function get_ic_ep_g

  real(rt) function get_ic_ep_s (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_ep_s = ic_ep_s(pID, solid_type)
    return
  end function get_ic_ep_s

  subroutine get_ic_pack_type (pID, c_string) bind(C)
    integer(c_int), intent(in) :: pID
    character(len=1, kind=c_char), intent(inout) :: c_string(16)
    integer :: N,I
    N = len_trim(ic_pack_type(pID))
    do I=1,N
      c_string(I) = ic_pack_type(pID)(I:I)
    enddo
    c_string(N+1) = c_null_char
  end subroutine get_ic_pack_type

  subroutine get_ic_dp_dist (pID, solid_type, c_string) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    character(len=1, kind=c_char), intent(inout) :: c_string(16)
    integer :: N,I
    N = len_trim(ic_dp_dist(pID, solid_type))
    do I=1,N
      c_string(I) = ic_dp_dist(pID, solid_type)(I:I)
    enddo
    c_string(N+1) = c_null_char
  end subroutine get_ic_dp_dist

  real(rt) function get_ic_dp_mean (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_dp_mean = ic_dp_mean(pID, solid_type)
    return
  end function get_ic_dp_mean

  real(rt) function get_ic_dp_std (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_dp_std = ic_dp_std(pID, solid_type)
    return
  end function get_ic_dp_std

  real(rt) function get_ic_dp_min (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_dp_min = ic_dp_min(pID, solid_type)
    return
  end function get_ic_dp_min

  real(rt) function get_ic_dp_max (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_dp_max = ic_dp_max(pID, solid_type)
    return
  end function get_ic_dp_max

  subroutine get_ic_ro_s_dist (pID, solid_type, c_string) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    character(len=1, kind=c_char), intent(inout) :: c_string(16)
    integer :: N,I
    N = len_trim(ic_ro_s_dist(pID, solid_type))
    do I=1,N
      c_string(I) = ic_ro_s_dist(pID, solid_type)(I:I)
    enddo
    c_string(N+1) = c_null_char
  end subroutine get_ic_ro_s_dist

  real(rt) function get_ic_ro_s_mean (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_ro_s_mean = ic_ro_s_mean(pID, solid_type)
    return
  end function get_ic_ro_s_mean

  real(rt) function get_ic_ro_s_std (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_ro_s_std = ic_ro_s_std(pID, solid_type)
    return
  end function get_ic_ro_s_std

  real(rt) function get_ic_ro_s_min (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_ro_s_min = ic_ro_s_min(pID, solid_type)
    return
  end function get_ic_ro_s_min

  real(rt) function get_ic_ro_s_max (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_ro_s_max = ic_ro_s_max(pID, solid_type)
    return
  end function get_ic_ro_s_max

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutines: get_ic_u_g, get_ic_v_g, get_ic_w_g & many others        !
!                                                                      !
! Purpose: Getters for the initial conditions values                   !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  real(rt) function get_ic_u_g (pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_ic_u_g = ic_u_g(pID)
    return
  end function get_ic_u_g

  real(rt) function get_ic_v_g (pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_ic_v_g = ic_v_g(pID)
    return
  end function get_ic_v_g

  real(rt) function get_ic_w_g (pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_ic_w_g = ic_w_g(pID)
    return
  end function get_ic_w_g

  real(rt) function get_ic_x_w (pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_ic_x_w = ic_x_w(pID)
    return
  end function get_ic_x_w

  real(rt) function get_ic_y_s (pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_ic_y_s = ic_y_s(pID)
    return
  end function get_ic_y_s

  real(rt) function get_ic_z_b (pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_ic_z_b = ic_z_b(pID)
    return
  end function get_ic_z_b

  real(rt) function get_ic_x_e (pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_ic_x_e = ic_x_e(pID)
    return
  end function get_ic_x_e

  real(rt) function get_ic_y_n (pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_ic_y_n = ic_y_n(pID)
    return
  end function get_ic_y_n

  real(rt) function get_ic_z_t (pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_ic_z_t = ic_z_t(pID)
    return
  end function get_ic_z_t

  real(rt) function get_ic_u_s (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_u_s = ic_u_s(pID, solid_type)
    return
  end function get_ic_u_s

  real(rt) function get_ic_v_s (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_v_s = ic_v_s(pID, solid_type)
    return
  end function get_ic_v_s

  real(rt) function get_ic_w_s (pID, solid_type) bind(C)
    integer(c_int), intent(in) :: pID, solid_type
    get_ic_w_s = ic_w_s(pID, solid_type)
    return
  end function get_ic_w_s

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: ic_defined                                               !
!                                                                      !
! Purpose: Return if a IC region has been defined based on coordinates !
! defined in the input deck.                                           !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  logical function ic_defined (icv)
    use param, only: is_defined
    integer, intent(in) :: icv

    ic_defined = &
         is_defined(ic_x_w(icv)) .or. is_defined(ic_x_e(icv)) .or. &
         is_defined(ic_y_s(icv)) .or. is_defined(ic_y_n(icv)) .or. &
         is_defined(ic_z_b(icv)) .or. is_defined(ic_z_t(icv))
   end function ic_defined

  integer(c_int) function ic_defined_cpp (icv) bind(C)
    use param, only: is_defined
    integer(c_int), intent(in) :: icv

    if (ic_defined(icv)) then
      ic_defined_cpp = 1
    else
      ic_defined_cpp = 0
   endif
   
   end function ic_defined_cpp


end module ic
