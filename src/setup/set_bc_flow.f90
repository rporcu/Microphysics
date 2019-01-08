module set_bc_flow_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

  use bc, only: bc_type, bc_plane

  use param,  only: dim_m
  use param, only: zero, one, equal, is_defined

  use bc, only: bc_u_g, bc_v_g, bc_w_g
  use bc, only: bc_massflow_g, bc_volflow_g

  use bc, only: bc_u_s, bc_v_s, bc_w_s
  use bc, only: bc_massflow_s, bc_volflow_s

  implicit none
  private

  public set_bc_flow

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_FLOW                                             !
!                                                                      !
!  Purpose: Check boundary condition specifications                    !
!     - convert physical locations to i, j, k's                        !
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine set_bc_flow(xlength, ylength, zlength, dx, dy, dz) &
     bind(C,name ="set_bc_flow")

     use param,    only: dim_bc
     use bc,       only: bc_defined, bc_type
     use bc,       only: bc_ep_s

    implicit none

    real(rt)  , intent(in) :: xlength, ylength, zlength
    real(rt)  , intent(in) :: dx, dy, dz

    integer :: bcv, i
    logical :: check(dim_m)  ! Flag to skip checks on indexed solid phase.

    ! Loop over each defined BC and check the user data.
    do bcv = 1, dim_bc

       if(bc_defined(bcv)) then

          ! Determine which solids phases are present.
          do i = 1, dim_m
             check(i) = (bc_ep_s(bcv,i) > zero)
          end do

          select case (trim(bc_type(bcv)))

             case ('MASS_INFLOW','MI','MASS_OUTFLOW','MO')
                call flow_to_vel(bcv, check, &
                                 xlength, ylength, zlength, dx, dy, dz)

          end select
       endif
    enddo

  end subroutine set_bc_flow

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: FLOW_TO_VEL                                             !
!                                                                      !
!  Purpose: Convert volumetric and mass flow rates to velocities       !
!     A specified mass flow rate is first converted to volumetric      !
!     flow rate. The volumetric flow rate is then converted to a       !
!     velocity.                                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine flow_to_vel(bcv, check, xlength, ylength, zlength, &
     dx, dy, dz)

      integer, intent(in) :: bcv
      logical, intent(in) :: check(dim_m)
      real(rt)  , intent(in) :: xlength, ylength, zlength
      real(rt)  , intent(in) :: dx, dy, dz

      integer :: m

      ! mass flows rates are converted to volumetric flow rates.
      if(is_defined(bc_massflow_g(bcv))) &
         call gas_massflow_to_volflow(bcv)

      ! volumetric flow rates are converted to velocities.
      if(is_defined(bc_volflow_g(bcv))) &
         call gas_volflow_to_velocity(bcv, xlength, ylength, zlength, dx, dy, dz)

      do m=1,dim_m
         if(check(m)) then
            if(is_defined(bc_massflow_s(bcv,m))) &
               call solids_massflow_to_volflow(bcv,m)

            if(is_defined(bc_volflow_s(bcv,m))) &
               call solids_volflow_to_velocity(bcv,m, &
               xlength, ylength, zlength, dx, dy, dz)
         endif
      enddo

   end subroutine flow_to_vel

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_MASSFLOW_TO_VOLFLOW                                 !
!                                                                      !
!  Purpose: Convert a gas phase BC input from a mass flow rate to      !
!  a volumetric flow rate.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine gas_massflow_to_volflow(bcv)

    use eos,       only: eosg
    use fld_const, only: ro_g0

    integer, intent(in) :: bcv
    real(rt)        :: volflow

    ! No need to convert if the mass flow is zero.
    if (equal(bc_massflow_g(bcv),zero)) then
       volflow = zero

    ! incompressible gas bc.
    else 

       volflow = bc_massflow_g(bcv)/ro_g0

    endif

! store the calculated volumetric flow rate.
    bc_volflow_g(bcv) = volflow

  end subroutine gas_massflow_to_volflow

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOLIDS_MASSFLOW_TO_VOLFLOW                              !
!                                                                      !
!  Purpose: Convert solids phase BC input from a mass flow rate to     !
!  a volumetric flow rate.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine solids_massflow_to_volflow(bcv,m)


      integer, intent(in) :: bcv, m

      bc_volflow_s(bcv,m) = 0.0d0
      stop 343143


  end subroutine solids_massflow_to_volflow


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_VOLFLOW_TO_VELOCITY                                 !
!                                                                      !
!  Purpose: Convert gas phase volumetric rate to a velocity.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine gas_volflow_to_velocity(bcv, xlength, ylength, zlength, &
     dx, dy, dz)

     use bc, only: bc_ep_g
     use bc, only: bc_x_w, bc_y_s, bc_z_b
     use bc, only: bc_x_e, bc_y_n, bc_z_t
     use calc_cell_module, only: calc_cell_bc_flow

     implicit none

     integer,        intent(in) :: bcv
     real(rt)  , intent(in) :: xlength, ylength, zlength
     real(rt)  , intent(in) :: dx, dy, dz

     real(rt) :: sgn, off, vel, area
     integer      :: i_w, i_e, j_s, j_n, k_b, k_t

    select case (trim(bc_type(bcv)))
    case ('MASS_INFLOW', 'MI'); SGN =  ONE; OFF = ZERO
    case ('MASS_OUTFLOW','MO'); SGN = -ONE; OFF = ONE
    end select

    select case (bc_plane(bcv))
    case ('W'); sgn = -sgn
    case ('S'); sgn = -sgn
    case ('B'); sgn = -sgn
    end select

    call calc_cell_bc_flow(&
       xlength, ylength, zlength, dx, dy, dz, &
       bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
       bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
       i_w, i_e, j_s, j_n, k_b, k_t)

    select case(bc_plane(bcv))
    case('W','E')
       area = dy*dble(j_n-j_s+1)*dz*dble(k_t-k_b+1)
       vel = sgn*bc_volflow_g(bcv)/(area*bc_ep_g(bcv))
       bc_u_g(bcv) = vel
       bc_v_g(bcv) = off * bc_v_g(bcv)
       bc_w_g(bcv) = off * bc_w_g(bcv)
    case('S','N')
       area = dx*dble(i_e-i_w+1)*dz*dble(k_t-k_b+1)
       vel = sgn*bc_volflow_g(bcv)/(area*bc_ep_g(bcv))
       bc_v_g(bcv) = vel
       bc_u_g(bcv) = off * bc_u_g(bcv)
       bc_w_g(bcv) = off * bc_w_g(bcv)
    case('B','T')
       area = dx*dble(i_e-i_w+1)*dy*dble(j_n-j_s+1)
       vel = sgn*bc_volflow_g(bcv)/(area*bc_ep_g(bcv))
       bc_w_g(bcv) = vel
       bc_u_g(bcv) = off * bc_u_g(bcv)
       bc_v_g(bcv) = off * bc_v_g(bcv)
    end select

  end subroutine gas_volflow_to_velocity

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOLIDS_VOLFLOW_TO_VELOCITY                              !
!                                                                      !
!  Purpose: Convert volumetric and mass flow rates to velocities       !
!     A specified mass flow rate is first converted to volumetric      !
!     flow rate. The volumetric flow rate is then converted to a       !
!     velocity.                                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine solids_volflow_to_velocity(bcv, m, &
     xlength, ylength, zlength, dx, dy, dz)


     use bc, only: bc_ep_s
     use bc, only: bc_x_w, bc_y_s, bc_z_b
     use bc, only: bc_x_e, bc_y_n, bc_z_t
     use calc_cell_module, only: calc_cell_bc_flow

     implicit none

     integer,        intent(in) :: bcv, m
     real(rt)  , intent(in) :: xlength, ylength, zlength
     real(rt)  , intent(in) :: dx, dy, dz

     real(rt) :: vel, sgn, off, area
     integer      :: i_w, i_e, j_s, j_n, k_b, k_t

     select case (trim(bc_type(bcv)))
     case ('MASS_INFLOW', 'MI'); sgn =  one; off = zero
     case ('MASS_OUTFLOW','MO'); sgn = -one; off = one
     end select

     select case (bc_plane(BCV))
     case ('W'); sgn = -sgn
     case ('S'); sgn = -sgn
     case ('B'); sgn = -sgn
     end select

    call calc_cell_bc_flow(&
       xlength, ylength, zlength, dx, dy, dz, &
       bc_x_w(bcv), bc_y_s(bcv), bc_z_b(bcv), &
       bc_x_e(bcv), bc_y_n(bcv), bc_z_t(bcv), &
       i_w, i_e, j_s, j_n, k_b, k_t)

     select case(bc_plane(bcv))
     case('W','E')
        area = dy*dble(j_n-j_s+1)*dz*dble(k_t-k_b+1)
        vel = sgn * bc_volflow_s(bcv,m)/(area*bc_ep_s(bcv,m))
        bc_u_s(bcv,m) = vel
        bc_v_s(bcv,m) = off * bc_v_s(bcv,m)
        bc_w_s(bcv,m) = off * bc_w_s(bcv,m)
     case('S','N')
        area = dx*dble(i_e-i_w+1)*dz*dble(k_t-k_b+1)
        vel = sgn * bc_volflow_s(bcv,m)/(area*bc_ep_s(bcv,m))
        bc_v_s(bcv,m) = vel
        bc_u_s(bcv,m) = off * bc_u_s(bcv,m)
        bc_w_s(bcv,m) = off * bc_w_s(bcv,m)
     case('B','Y')
        area = dx*dble(i_e-i_w+1)*dy*dble(j_n-j_s+1)
        vel = sgn * bc_volflow_s(bcv,m)/(area*bc_ep_s(bcv,m))
        bc_w_s(bcv,m) = vel
        bc_u_s(bcv,m) = off * bc_u_s(bcv,m)
        bc_v_s(bcv,m) = off * bc_v_s(bcv,m)
     end select

  end subroutine solids_volflow_to_velocity
end module set_bc_flow_module
