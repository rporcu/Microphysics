!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: bc                                                     !
!  Purpose: Global variables for specifying boundary conditions.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module bc

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int

  use param, only: dim_bc, dim_m, dim_n_g, dim_n_s

  ! Type of boundary:
  character(len=16) :: BC_Type(dim_bc)

  ! Flags for periodic boundary conditions
  logical :: cyclic_x = .false.
  logical :: cyclic_y = .false.
  logical :: cyclic_z = .false.

  ! Boundary condition coordinates
  real(c_real) :: BC_X_w(dim_bc), BC_X_e(dim_bc)
  real(c_real) :: BC_Y_s(dim_bc), BC_Y_n(dim_bc)
  real(c_real) :: BC_Z_b(dim_bc), BC_Z_t(dim_bc)

  real(c_real) :: BC_Normal(dim_bc,3)
  real(c_real) :: BC_Center(dim_bc,3)

  ! Void fraction in a specified boundary
  real(c_real) :: BC_EP_g(dim_bc), BC_EP_s(dim_bc, dim_m)

  ! Gas phase BC pressure
  real(c_real) :: BC_P_g(dim_bc)

  ! Velocities at a specified boundary
  real(c_real) :: BC_U_g(dim_bc), BC_U_s(dim_bc, dim_m)
  real(c_real) :: BC_V_g(dim_bc), BC_V_s(dim_bc, dim_m)
  real(c_real) :: BC_W_g(dim_bc), BC_W_s(dim_bc, dim_m)

  ! Volumetric flow rate through a mass inflow boundary
  real(c_real) :: BC_VolFlow_g(dim_bc), BC_VolFlow_s(dim_bc, dim_m)

  ! Mass flow rate through a mass inflow boundary
  real(c_real) :: BC_MassFlow_g(dim_bc), BC_MassFlow_s(dim_bc, dim_m)

  ! Specified pressure drop cyclic boundary
  real(c_real) :: delp_x, delp_y, delp_z

  ! Partial slip wall boundary condition (gas only)
  real(c_real) :: BC_hw_g(dim_bc)
  real(c_real) :: BC_Uw_g(dim_bc)
  real(c_real) :: BC_Vw_g(dim_bc)
  real(c_real) :: BC_Ww_g(dim_bc)

  ! Heat transfer boundary condition
  real(c_real) :: BC_T_g   (dim_bc), BC_T_s   (dim_bc, dim_m)
  real(c_real) :: BC_hw_T_g(dim_bc), BC_hw_T_s(dim_bc, dim_m)
  real(c_real) :: BC_Tw_g  (dim_bc), BC_Tw_s  (dim_bc, dim_m)
  real(c_real) :: BC_C_T_g (dim_bc), BC_C_T_s (dim_bc, dim_m)

  ! Species transfer boundary condition
  real(c_real) :: BC_X_g(dim_bc, dim_n_g), BC_X_s(dim_bc, dim_m, dim_n_s)
  real(c_real) :: BC_hw_X_g(dim_bc, dim_n_g)
  real(c_real) :: BC_Xw_g  (dim_bc, dim_n_g)
  real(c_real) :: BC_C_X_g (dim_bc, dim_n_g)

  ! Character variable to determine the flow plane of a flow cell
  character :: BC_Plane(dim_bc)

  ! Cell flag definitions
  integer, parameter :: undef_cell =   0 ! undefined
  integer, parameter :: pinf_      =  10 ! pressure inflow cell
  integer, parameter :: pout_      =  11 ! pressure outflow cell
  integer, parameter :: minf_      =  20 ! mass flux inflow cell
  integer, parameter :: nsw_       = 100 ! wall with no-slip b.c.
  integer, parameter :: fsw_       = 101 ! wall with free-slip
  integer, parameter :: psw_       = 102 ! wall with partial-slip b.c.
  integer, parameter :: cycl_      =  50 ! cyclic b.c.
  integer, parameter :: cycp_      =  51 ! cyclic b.c. with pressure drop


contains

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
! Subroutine: bc_defined                                               !
!                                                                      !
! Purpose: Return if a BC region has been defined based on coordinates !
! defined in the input deck.                                           !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  logical function bc_defined(icv)

    use param, only: is_defined

    integer, intent(in) :: icv

    bc_defined = is_defined(bc_x_w(icv)) .or. is_defined(bc_x_e(icv)) .or. &
                 is_defined(bc_y_s(icv)) .or. is_defined(bc_y_n(icv)) .or. &
                 is_defined(bc_z_b(icv)) .or. is_defined(bc_z_t(icv))

! An IC is defined for restart runs only if it is a 'PATCH'.
    if(bc_type(icv) == 'DUMMY') bc_defined = .false.

  end function bc_defined



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
    real(c_real)  , intent(in) :: dx, dy, dz
    real(c_real)  , intent(in) :: xlength, ylength, zlength
    integer(c_int), intent(in) :: domlo(3), domhi(3)

    integer :: bcv, m
    integer :: i_w, j_s, k_b
    integer :: i_e, j_n, k_t

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
             call calc_cell_bc_flow(&
               xlength, ylength, zlength, dx, dy, dz, &
               bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
               bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
               i_w, i_e, j_s, j_n, k_b, k_t)
             flow_bc = .true.
          case ('MASS_OUTFLOW','MO')
             write (unit_out,"(9x,'Outlet with specified mass flux')")
             call calc_cell_bc_flow(&
               xlength, ylength, zlength, dx, dy, dz, &
               bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
               bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
               i_w, i_e, j_s, j_n, k_b, k_t)
          case ('P_INFLOW','PI')
             write (unit_out,"(9x,'Inlet with specified gas pressure')")
             call calc_cell_bc_flow(&
               xlength, ylength, zlength, dx, dy, dz, &
               bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
               bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
               i_w, i_e, j_s, j_n, k_b, k_t)
             flow_bc = .true.
          case ('P_OUTFLOW','PO')
             write (unit_out,"(9x,'Outlet with specified gas pressure')")
             call calc_cell_bc_flow(&
               xlength, ylength, zlength, dx, dy, dz, &
               bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
               bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
               i_w, i_e, j_s, j_n, k_b, k_t)
             flow_bc = .true.
          case ('FREE_SLIP_WALL','FSW')
             write (unit_out,"(9x,'Velocity gradients are zero')")
             call calc_cell_bc_wall(domlo, domhi, &
               xlength, ylength, zlength, dx, dy, dz, &
               bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
               bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
               i_w, i_e, j_s, j_n, k_b, k_t)
             flow_bc = .false.
          case ('NO_SLIP_WALL','NSW')
             write (unit_out,"(9x,'Velocity is zero at wall')")
             call calc_cell_bc_wall(domlo, domhi, &
               xlength, ylength, zlength, dx, dy, dz, &
               bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
               bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
               i_w, i_e, j_s, j_n, k_b, k_t)
             flow_bc = .false.
          case ('PAR_SLIP_WALL','PSW')
             write (unit_out,"(9x,'Partial slip condition at wall')")
             call calc_cell_bc_wall(domlo, domhi, &
               xlength, ylength, zlength, dx, dy, dz, &
               bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
               bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
               i_w, i_e, j_s, j_n, k_b, k_t)
             flow_bc = .false.
          end select

          write (unit_out, 1620) &
            bc_x_w(bcv), dx*dble(i_w-1), bc_x_e(bcv), dx*dble(i_e), &
            bc_y_s(bcv), dy*dble(j_s-1), bc_y_n(bcv), dy*dble(j_n), &
            bc_z_b(bcv), dz*dble(k_b-1), bc_z_t(bcv), dz*dble(k_t)

1620  format(9x,45X,' Specified  ',5X,' Simulated  ',/&
         9X,'X coordinate of west face   (BC_X_w) ...... ',g12.5, 5x, g12.5/,&
         9x,'X coordinate of east face   (BC_X_e) ...... ',g12.5, 5x, g12.5/,&
         9x,'Y coordinate of south face  (BC_Y_s) ...... ',g12.5, 5x, g12.5/,&
         9x,'Y coordinate of north face  (BC_Y_n) ...... ',g12.5, 5x, g12.5/,&
         9x,'Z coordinate of bottom face (BC_Z_b) ...... ',g12.5, 5x, g12.5/,&
         9x,'Z coordinate of top face    (BC_Z_t) ...... ',g12.5, 5x, g12.5/)

          write (unit_out, 1630) i_w, i_e, j_s, j_n, k_b, k_t

1630  format(&
         9X,'I index of cell at west   (BC_I_w) ',24('.'),1x,I4,/,&
         9X,'I index of cell at east   (BC_I_e) ',24('.'),1x,I4,/,&
         9X,'J index of cell at south  (BC_J_s) ',24('.'),1x,I4,/,&
         9X,'J index of cell at north  (BC_J_n) ',24('.'),1x,I4,/,&
         9X,'K index of cell at bottom (BC_K_b) ',24('.'),1x,I4,/,&
         9X,'K index of cell at top    (BC_K_t) ',24('.'),1x,I4)


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

          else
             if (bc_type(bcv) == 'PAR_SLIP_WALL' .or. bc_type(bcv) == 'PSW') &
               write (unit_out, 1675) bc_hw_g(bcv), &
               bc_uw_g(bcv), bc_vw_g(bcv), bc_ww_g(bcv)

1675  format(9X,'Partial slip coefficient (BC_hw_g) .... ',G12.5,/,&
             9X,'Slip velocity U at wall (BC_Uw_g) ..... ',G12.5,/,&
             9X,'Slip velocity V at wall (BC_Vw_g) ..... ',G12.5,/,&
             9X,'Slip velocity W at wall (BC_Ww_g) ..... ',G12.5)

          endif
       endif
    enddo

    return
  end subroutine write_out_bc

end module bc
