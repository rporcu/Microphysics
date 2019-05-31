!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: bc                                                     !
!  Purpose: Global variables for specifying boundary conditions.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module bc

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

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
  real(rt) :: BC_EP_g(dim_bc), BC_EP_s(dim_bc, dim_m)

  ! Gas phase BC pressure
  real(rt) :: BC_P_g(dim_bc)

  ! Velocities at a specified boundary
  real(rt) :: BC_U_g(dim_bc), BC_U_s(dim_bc, dim_m)
  real(rt) :: BC_V_g(dim_bc), BC_V_s(dim_bc, dim_m)
  real(rt) :: BC_W_g(dim_bc), BC_W_s(dim_bc, dim_m)
  
  ! Array containing BC_U_g, BC_V_g, BC_W_g
  real(rt) :: BC_Vel_g(dim_bc,3)

  ! Volumetric flow rate through a mass inflow boundary
  real(rt) :: BC_VolFlow_g(dim_bc), BC_VolFlow_s(dim_bc, dim_m)

  ! Mass flow rate through a mass inflow boundary
  real(rt) :: BC_MassFlow_g(dim_bc), BC_MassFlow_s(dim_bc, dim_m)

  ! Specified pressure drop cyclic boundary
  real(rt) :: delp_x, delp_y, delp_z

  ! Partial slip wall boundary condition (gas only)
  real(rt) :: BC_hw_g(dim_bc)
  real(rt) :: BC_Uw_g(dim_bc)
  real(rt) :: BC_Vw_g(dim_bc)
  real(rt) :: BC_Ww_g(dim_bc)

  ! Heat transfer boundary condition
  real(rt) :: BC_T_g   (dim_bc), BC_T_s   (dim_bc, dim_m)
  real(rt) :: BC_hw_T_g(dim_bc), BC_hw_T_s(dim_bc, dim_m)
  real(rt) :: BC_Tw_g  (dim_bc), BC_Tw_s  (dim_bc, dim_m)
  real(rt) :: BC_C_T_g (dim_bc), BC_C_T_s (dim_bc, dim_m)

  ! Species transfer boundary condition
  real(rt) :: BC_X_g(dim_bc, dim_n_g), BC_X_s(dim_bc, dim_m, dim_n_s)
  real(rt) :: BC_hw_X_g(dim_bc, dim_n_g)
  real(rt) :: BC_Xw_g  (dim_bc, dim_n_g)
  real(rt) :: BC_C_X_g (dim_bc, dim_n_g)

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
! Subroutines: get_bc_u_g, get_bc_v_g, get_bc_w_g, get_bc_t_g          !
!              get_bc_ep_g                                             !
!                                                                      !
! Purpose: Getters for the boundary conditions values                  !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
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

  real(rt) function get_bc_vel_g(n,pID) bind(C)
    integer(c_int), intent(in) :: n,pID
    get_bc_vel_g = bc_vel_g(n,pID)
    return
  end function get_bc_vel_g

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


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: write_out_bc                                            !
!                                                                      !
!  Purpose: Echo user input for BC regions.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine write_out_bc(unit_out, dx, dy, dz, &
    xlength, ylength, zlength, domlo, domhi)

    use param, only: zero, is_defined

    use calc_cell_module, only: calc_cell_bc_flow
    use calc_cell_module, only: calc_cell_bc_wall

    implicit none

    integer,        intent(in) :: unit_out
    real(rt)  , intent(in) :: dx, dy, dz
    real(rt)  , intent(in) :: xlength, ylength, zlength
    integer(c_int), intent(in) :: domlo(3), domhi(3)

    integer :: bcv, m

    logical :: flow_bc

!-----------------------------------------------


! Boundary Condition Data
    write (unit_out, 1600)
1600  format(//,3x,'7. BOUNDARY CONDITIONS')

    if (cyclic_x .and. abs(delp_x) > epsilon(0.0d0)) then
       write (unit_out, 1602) 'X', ' with pressure drop'
       write (unit_out, 1603) 'X', DELP_X
    else if (cyclic_x) then
       write (unit_out, 1602) 'X'
    endif
    if (cyclic_y .and. abs(delp_y) > epsilon(0.0d0)) then
       write (unit_out, 1602) 'Y', ' with pressure drop'
       write (unit_out, 1603) 'Y', DELP_Y
    else if (cyclic_y) then
       write (unit_out, 1602) 'Y'
    endif
    if (cyclic_z .and. abs(delp_z) > epsilon(0.0d0)) then
       write (unit_out, 1602) 'Z', ' with pressure drop'
       write (unit_out, 1603) 'Z', DELP_Z
    else if (cyclic_z) then
       write (unit_out, 1602) 'Z'
    endif

 1602 format(/7X,'Cyclic boundary conditions in ',A,' direction',A)
 1603 format( 7X,'Pressure drop (DELP_',A,') = ',G12.5)

    do bcv = 1, dim_bc
       if (bc_defined(bcv)) then

          write (unit_out, 1610) bcv, bc_type(bcv)

1610  format(/7x,'Boundary condition no : ',I4,2/&
              9X,'Type of boundary condition : ',A16)

          select case (trim(bc_type(bcv)))
          case ('MASS_INFLOW','MI')
             write (unit_out,"(9x,'Inlet with specified mass flux')")
             flow_bc = .true.
          case ('MASS_OUTFLOW','MO')
             write (unit_out,"(9x,'Outlet with specified mass flux')")
          case ('P_INFLOW','PI')
             write (unit_out,"(9x,'Inlet with specified gas pressure')")
             flow_bc = .true.
          case ('P_OUTFLOW','PO')
             write (unit_out,"(9x,'Outlet with specified gas pressure')")
             flow_bc = .true.
          case ('NO_SLIP_WALL','NSW')
             write (unit_out,"(9x,'Velocity is zero at wall')")
             flow_bc = .false.
          end select

          if(flow_bc) then
             write (unit_out, "(' ')")
             write (unit_out,1640) bc_ep_g(bcv)
             if(is_defined(bc_p_g(bcv))) &
               write (unit_out,1641) bc_p_g(bcv)
             if(is_defined(bc_t_g(bcv))) &
               write (unit_out,1642) bc_t_g(bcv)
             if(is_defined(bc_massflow_g(bcv))) &
               write (unit_out,1648) bc_massflow_g(bcv)
             if(is_defined(bc_volflow_g(bcv)))  &
               write (unit_out,1649) bc_volflow_g(bcv)
             write (unit_out, 1650) bc_u_g(bcv)
             write (unit_out, 1651) bc_v_g(bcv)
             write (unit_out, 1652) bc_w_g(bcv)

1640  format(9X,'Gas phase volume fraction (BC_EP_g) ....... ',g12.5)
1641  format(9X,'Gas pressure (BC_P_g) ..................... ',g12.5)
1642  format(9X,'Gas temperature (BC_T_g) .................. ',g12.5)
1648  format(9X,'Gas mass flow rate (BC_MassFlow_g) ........ ',g12.5)
1649  format(9X,'Gas volumetric flow rate (BC_VOLFLOW_g) ... ',g12.5)
1650  format(9X,'X-component of gas velocity (BC_U_g) ...... ',g12.5)
1651  format(9X,'Y-component of gas velocity (BC_V_g) ...... ',g12.5)
1652  format(9X,'Z-component of gas velocity (BC_W_g) ...... ',g12.5)

             do m = 1, dim_m
                if(bc_ep_s(bcv,m) > zero) then
                   write (unit_out, "(' ')")
                   write (unit_out, 1660) m, bc_ep_s(bcv,m)
                   write (unit_out, "(' ')")
                   if(is_defined(bc_massflow_s(bcv,m))) &
                     write(unit_out, 1668) m, bc_massflow_s(bcv,m)
                   if(is_defined(bc_volflow_s(bcv,m)))  &
                     write(unit_out, 1669) m, bc_volflow_s(bcv,m)
                   write(unit_out,1670)m,bc_u_s(bcv,m)
                   write(unit_out,1671)m,bc_v_s(bcv,m)
                   write(unit_out,1672)m,bc_w_s(bcv,m)
                endif
             enddo

1660  format(9X,'Solids phase-',I2,' Volume fraction (BC_EP_s) ............. ',g12.5)
1668  format(9X,'Solids phase-',I2,' mass flow rate (BC_MASSFLOW_s) ........ ',g12.5)
1669  format(9X,'Solids phase-',I2,' volumetric flow rate (BC_VOLFLOW_s) ... ',g12.5)
1670  format(9X,'X-component of solids phase-',I2,' velocity (BC_U_s) ...... ',g12.5)
1671  format(9X,'Y-component of solids phase-',I2,' velocity (BC_V_s) ...... ',g12.5)
1672  format(9X,'Z-component of solids phase-',I2,' velocity (BC_W_s) ...... ',g12.5)

          endif
       endif
    enddo

    return
  end subroutine write_out_bc

end module bc
