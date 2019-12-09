!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: bc                                                     !
!  Purpose: Global variables for specifying boundary conditions.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module bc

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int, c_char, c_null_char

  use param, only: undefined
  use param, only: dim_bc, dim_m, dim_n_g, dim_n_s

  ! Type of boundary:
  character(len=16) :: BC_Type(dim_bc)

  ! Flags for periodic boundary conditions
  logical :: cyclic_x = .false.
  logical :: cyclic_y = .false.
  logical :: cyclic_z = .false.

  logical :: bc_defined(1:dim_bc) = .false.

  ! Boundary condition location (EB planes)
  real(rt) :: BC_Normal(1:dim_bc,1:3) = undefined
  real(rt) :: BC_Center(1:dim_bc,1:3) = undefined

  ! Void fraction in a specified boundary
  real(rt) :: BC_EP_g(dim_bc) = 1.0d0
  real(rt) :: BC_EP_s(dim_bc, dim_m) = 0.0d0

  ! Gas phase BC pressure
  real(rt) :: BC_P_g(dim_bc) = undefined

  ! Array containing BC_U_g, BC_V_g, BC_W_g
  real(rt) :: BC_U_g(dim_bc) = 0.0d0
  real(rt) :: BC_V_g(dim_bc) = 0.0d0
  real(rt) :: BC_W_g(dim_bc) = 0.0d0

  ! Heat transfer boundary condition
  real(rt) :: BC_T_g   (dim_bc)        = 293.15d0


  ! External shaking (shaking amplitude vector sets direction of shaking)
  real(rt), dimension(3) :: BC_shaker_A ! shaking amplitude
  real(rt)               :: BC_shaker_F ! shaking frequency

  ! Character variable to determine the flow plane of a flow cell
  character :: BC_Plane(dim_bc)

  ! Cell flag definitions
  integer, parameter :: undef_cell =   0 ! undefined
  integer, parameter :: ignore_    =   9 ! ignore this wall because it is outside an EB boundary
  integer, parameter :: pinf_      =  10 ! pressure inflow cell
  integer, parameter :: pout_      =  11 ! pressure outflow cell
  integer, parameter :: minf_      =  20 ! mass flux inflow cell
  integer, parameter :: nsw_       = 100 ! wall with no-slip b.c.

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutines: getters                                                 !
!                                                                      !
! Purpose: Getters for the boundary conditions values                  !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  integer(c_int) function get_bc_defined(pID) bind(C)
    integer(c_int), intent(in) :: pID
    if(bc_defined(pID)) then
      get_bc_defined = 1
    else
      get_bc_defined = 0
    endif
    return
  end function get_bc_defined

  subroutine get_bc_type(pID, c_string) bind(C)
    integer(c_int), intent(in) :: pID
    character(len=1, kind=c_char), intent(inout) :: c_string(16)
    integer :: N,I
    N = len_trim(BC_Type(pID))
    do I=1,N
      c_string(I) = BC_Type(pID)(I:I)
    enddo
    c_string(N+1) = c_null_char
  end subroutine get_bc_type

  real(rt) function get_bc_u_g(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_u_g = bc_u_g(pID)
    return
  end function get_bc_u_g

  real(rt) function get_bc_v_g(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_v_g = bc_v_g(pID)
    return
  end function get_bc_v_g

  real(rt) function get_bc_w_g(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_w_g = bc_w_g(pID)
    return
  end function get_bc_w_g

  real(rt) function get_bc_t_g(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_t_g = bc_t_g(pID)
    return
  end function get_bc_t_g

  real(rt) function get_bc_ep_g(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_ep_g = bc_ep_g(pID)
    return
  end function get_bc_ep_g

  real(rt) function get_bc_p_g(pID) bind(C)
    integer(c_int), intent(in) :: pID
    get_bc_p_g = bc_p_g(pID)
    return
  end function get_bc_p_g

  integer(c_int) function get_minf() bind(C)
    get_minf = minf_
    return
  end function get_minf

  integer(c_int) function get_pinf() bind(C)
    get_pinf = pinf_
    return
  end function get_pinf

  integer(c_int) function get_pout() bind(C)
    get_pout = pout_
    return
  end function get_pout

  subroutine get_domain_bc (domain_bc_out) bind(C)
    integer(c_int), intent(out)  :: domain_bc_out(6)
    integer :: bcv
    ! Default is that we reflect particles off domain boundaries if not periodic
    domain_bc_out(1:6) = 1
    if (cyclic_x) domain_bc_out(1:2) = 0
    if (cyclic_y) domain_bc_out(3:4) = 0
    if (cyclic_z) domain_bc_out(5:6) = 0

    do bcv = 1,6
       select case (trim(bc_type(bcv)))
         case ('P_OUTFLOW','PO')
            domain_bc_out(bcv) = 0
       end select
    end do
  end subroutine get_domain_bc

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: set_cyclic                                               !
!                                                                      !
! Purpose: Function to set cyclic flags.                               !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  subroutine set_cyclic(cyc_x, cyc_y, cyc_z) &
    bind(C, name="mfix_set_cyclic")

    integer, intent(in) :: cyc_x, cyc_y, cyc_z

    cyclic_x = (cyc_x == 1)
    cyclic_y = (cyc_y == 1)
    cyclic_z = (cyc_z == 1)

  end subroutine set_cyclic


end module bc
