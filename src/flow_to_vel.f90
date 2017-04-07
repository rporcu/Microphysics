module flow_to_vel_new_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use exit_mod,       only: mfix_exit
  use toleranc,       only: compare
  use param,          only: DIM_M
  use run,            only: IFILE_NAME
  use param1,         only: UNDEFINED, ZERO, ONE, IS_DEFINED, IS_UNDEFINED
  use error_manager,  only: finl_err_msg, err_msg, flush_err_msg, &
                          & init_err_msg, ivar

  use bc, only: bc_type, bc_plane

  use bc, only: bc_massflow_g,  bc_volflow_g,  &
              & bc_ep_g, bc_u_g, bc_v_g, bc_w_g

  use bc, only: bc_massflow_s, bc_volflow_s, &
              &  bc_ep_s, bc_u_s, bc_v_s, bc_w_s

  implicit none
  private

  public flow_to_vel

contains

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
   subroutine flow_to_vel(m_tot, skip, bcv, xlength, ylength, zlength, dx, dy, dz)

      logical, intent(in) :: skip(dim_m)
      integer, intent(in) :: m_tot, bcv
      real(c_real)  , intent(in) :: xlength, ylength, zlength
      real(c_real)  , intent(in) :: dx, dy, dz

      integer :: m

      call init_err_msg("FLOW_TO_VEL_NEW")

      ! mass flows rates are converted to volumetric flow rates.
      if(is_defined(bc_massflow_g(bcv))) &
         call gas_massflow_to_volflow(bcv)

      do m=1,m_tot
         if(is_defined(bc_massflow_s(bcv,m))) &
            call solids_massflow_to_volflow(bcv,m,skip(m))
      enddo

      ! volumetric flow rates are converted to velocities.
      if(is_defined(bc_volflow_g(bcv))) &
         call gas_volflow_to_velocity(bcv, xlength, ylength, zlength, dx, dy, dz)

      do m=1,m_tot
         if(is_defined(bc_volflow_s(bcv,m))) &
            call solids_volflow_to_velocity(bcv,m,skip(m), &
               xlength, ylength, zlength, dx, dy, dz)
      enddo

      call finl_err_msg

   end subroutine flow_to_vel

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_MASSFLOW_TO_VOLFLOW                                 !
!  Author: M. Syamlal                                 Date: 28-JUL-92  !
!                                                                      !
!  Purpose: Convert a gas phase BC input from a mass flow rate to      !
!  a volumetric flow rate.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine gas_massflow_to_volflow(bcv)

    use bc,        only: bc_massflow_g, bc_p_g, bc_t_g, bc_volflow_g
    use eos,       only: eosg
    use fld_const, only: mw_avg, ro_g0
    use scales   , only: p_ref

    integer, intent(in) :: bcv
    real(c_real)        :: volflow
    real(c_real)        :: mw

    call init_err_msg("GAS_MASSFLOW_TO_VOLFLOW")

    ! No need to convert if the mass flow is zero.
    if(compare(bc_massflow_g(bcv),zero)) then
       volflow = zero

    ! incompressible gas bc.
    elseif(is_defined(ro_g0)) then
       volflow = bc_massflow_g(bcv)/ro_g0

    ! well-defined compresible gas bc.
    elseif(is_defined(bc_p_g(bcv)) .and. is_defined(bc_t_g(bcv))) then
       mw = mw_avg
       volflow = bc_massflow_g(bcv) / &
          eosg(mw,(bc_p_g(bcv)-p_ref),bc_t_g(bcv))

    ! fails. this shouldn't happen as previous checks should catch any
    ! errors leading to this routine.
    else
       write(err_msg, 1100) bcv
       call flush_err_msg(abort=.true.)

1100   format('Error 1100: Boundary condition ',I3,' failed sanity ',   &
            'check.',/'Please report this to the MFIX mailing list.')

    endif

    ! Check that a specified volumetric flow matches the calculated value.
    if(is_defined(bc_volflow_g(bcv))) then
       if(.not.compare(volflow,bc_volflow_g(bcv))) then
          write(err_msg,1101) trim(ivar('BC_MASSFLOW_g',bcv)), bcv,  &
               volflow, bc_volflow_g(bcv), trim(IFILE_NAME)
          call flush_err_msg(abort=.true.)
       endif
    else

       ! store the calculated volumetric flow rate.
       bc_volflow_g(bcv) = volflow
    endif

1101 FORMAT('Error 1101: Volumetric flow rate calculated from ',A,/   &
         'does NOT equal the specified volumetric flow rate for BC',I3,&
         /3x,'>>> Calculated: ',G14.7,/3x,'>>> Specified:  ',G14.7,/   &
         'Please correct the ',A,' file.')

    ! clean up and return.
    call finl_err_msg

  end subroutine gas_massflow_to_volflow



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOLIDS_MASSFLOW_TO_VOLFLOW                              !
!                                                                      !
!  Purpose: Convert solids phase BC input from a mass flow rate to     !
!  a volumetric flow rate.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine solids_massflow_to_volflow(bcv,m, skip_m)

      use bc,       only: bc_massflow_s, bc_volflow_s
      use constant, only: ro_s0

      integer, intent(in) :: bcv, m
      logical, intent(in) :: skip_m
      real(c_real)        :: volflow

      call init_err_msg("SOLIDS_MASSFLOW_TO_VOLFLOW")


      if(skip_m) then
         write(err_msg,1100) m, bcv, trim(ivar("BC_MASSFLOW_S",bcv,m)), &
            trim(IFILE_NAME)
         call flush_err_msg(abort=.true.)
      endif

1100 FORMAT('Error 1100: Solids phase ',I2,' has a specified mass ',  &
         'flow rate',/'at BC ',I3,', ',A,'. But, both BC_ROP_s and ',&
         'BC_EP_s are zero or undefined.',/'Please correct the ',&
         A,' file.')

      VOLFLOW = BC_MASSFLOW_S(BCV,M)/RO_S0(M)

    ! If volumetric flow is also specified compare both
      if(is_defined(bc_volflow_s(bcv,m))) then
         if(.not.compare(volflow,bc_volflow_s(bcv,m))) then
            write(err_msg,1101) trim(ivar('BC_MASSFLOW_S',bcv,m)), bcv, &
               volflow, bc_volflow_s(bcv,m),  trim(IFILE_NAME)
            call flush_err_msg(abort=.true.)
         endif
      else
         bc_volflow_s(bcv,m) = volflow
      endif

      call finl_err_msg

1101 format('Error 1101: Volumetric flow rate calculated from ',A,/   &
         'does NOT equal the specified volumetric flow rate for BC',I3,&
         /3x,'>>> Calculated: ',G14.7,/3x,'>>> Specified:  ',G14.7,/   &
         'Please correct the ',A,' file.')

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

     use compar, only: mype
     use calc_cell_module, only: calc_cell_bc_flow

     implicit none

     integer,        intent(in) :: bcv
     real(c_real)  , intent(in) :: xlength, ylength, zlength
     real(c_real)  , intent(in) :: dx, dy, dz

     real(c_real) :: sgn, off, vel, area
     integer      :: i_w, i_e, j_s, j_n, k_b, k_t

    call init_err_msg("GAS_VOLFLOW_TO_VELOCITY")

    select case (trim(bc_type(bcv)))
    case ('MASS_INFLOW','MI');  SGN =  ONE; OFF = ZERO
    case ('MASS_OUTFLOW','MO'); SGN = -ONE; OFF = ONE
    case DEFAULT
       write(*,*) 'error in GAS_VOLFLOW_TO_VELOCITY'
       call mfix_exit(myPE)
    end select

    select case (bc_plane(bcv))
    case ('W'); sgn = -sgn
    case ('S'); sgn = -sgn
    case ('B'); sgn = -sgn
    end select

    call calc_cell_bc_flow(bcv, xlength, ylength, zlength, &
       dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)

    select case(bc_plane(bcv))
    case('W','E'); area = dy*dble(j_n-j_s+1)*dz*dble(k_t-k_b+1)
    case('S','N'); area = dx*dble(i_e-i_w+1)*dz*dble(k_t-k_b+1)
    case('B','Y'); area = dx*dble(i_e-i_w+1)*dy*dble(j_n-j_s+1)
    end select

    ! BC area and BC volume fraction.
    vel = sgn*bc_volflow_g(bcv)/(area*bc_ep_g(bcv))

    ! if the user also defined the boundary velocity through the plane, then
    ! check that the calculated value agrees with the specified value. if
    ! the user did not define the boundary velocity through the plane, then
    ! if mass_inflow set the value of the boundary velocity to the
    ! calculated value. otherwise do nothing.
    if(bc_plane(BCV) == 'W' .or. bc_plane(BCV)== 'E') then

       if(is_defined(bc_u_g(bcv))) then
          if(.not.compare(vel,bc_u_g(bcv))) then
             write(err_msg,1100) bcv, vel, 'BC_U_g', bc_u_g(bcv)
             call flush_err_msg(abort=.true.)
          endif
       else
          bc_u_g(bcv) = vel
          bc_v_g(bcv) = off * bc_v_g(bcv)
          bc_w_g(bcv) = off * bc_w_g(bcv)
       endif

    elseif(bc_plane(BCV) == 'S' .or. bc_plane(BCV)== 'N') then
       if(is_defined(bc_v_g(bcv))) then
          if(.not.compare(vel,bc_v_g(bcv))) then
             write(err_msg, 1100) bcv, vel, 'BC_V_g', bc_v_g(bcv)
             call flush_err_msg(abort=.true.)
          endif
       else
          bc_v_g(bcv) = vel
          bc_u_g(bcv) = off * bc_u_g(bcv)
          bc_w_g(bcv) = off * bc_w_g(bcv)
       endif

    elseif(bc_plane(BCV) == 'B' .or. bc_plane(BCV)== 'T') then
       if(is_defined(bc_w_g(bcv))) then
          if(.not.compare(vel, bc_w_g(bcv))) then
             write(err_msg, 1100) bcv, vel, 'BC_W_g', bc_w_g(bcv)
             call flush_err_msg(abort=.true.)
          endif
       else
          bc_w_g(bcv) = vel
          bc_u_g(bcv) = off * bc_u_g(bcv)
          bc_v_g(bcv) = off * bc_v_g(bcv)
       endif

    endif

    call finl_err_msg


1100 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,') = ',G14.7,/1X,70('*')/)

  end subroutine gas_volflow_to_velocity

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: SOLIDS_VOLFLOW_TO_VELOCITY                              !
  !  Author: M. Syamlal                                 Date: 28-JUL-92  !
  !                                                                      !
  !  Purpose: Convert volumetric and mass flow rates to velocities       !
  !     A specified mass flow rate is first converted to volumetric      !
  !     flow rate. The volumetric flow rate is then converted to a       !
  !     velocity.                                                        !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine solids_volflow_to_velocity(bcv, m, skip_m, &
     xlength, ylength, zlength, dx, dy, dz)

     use calc_cell_module, only: calc_cell_bc_flow

     implicit none

     integer,        intent(in) :: bcv, m
     logical,        intent(in) :: skip_m
     real(c_real)  , intent(in) :: xlength, ylength, zlength
     real(c_real)  , intent(in) :: dx, dy, dz

     real(c_real) :: vel, sgn, off, area
     integer      :: i_w, i_e, j_s, j_n, k_b, k_t

     call init_err_msg("SOLIDS_VOLFLOW_TO_VELOCITY")

     if(skip_m) then
        write(err_msg,1100) m, bcv, trim(ivar("BC_VOLFLOW_S",bcv,m)), &
           & trim(IFILE_NAME)
        call flush_err_msg(abort=.true.)
     endif

1100 format('Error 1100: Solids phase ',I2,' has a specified ',       &
         'volumetric flow rate',/'at BC ',I3,', ',A,'. But, both ',&
         'BC_ROP_s and BC_EP_s are zero or undefined.',/'Please ',&
         'the ',A,' file.')

     select case (trim(bc_type(bcv)))
     case ('MASS_INFLOW', 'MI'); sgn =  one; off = zero
     case ('MASS_OUTFLOW','MO'); sgn = -one; off = one
     case DEFAULT
        write(*,*) 'error in SOLIDS_VOLFLOW_TO_VELOCITY'
        call mfix_exit(0)
     end select

     select case (bc_plane(BCV))
     case ('W'); sgn = -sgn
     case ('S'); sgn = -sgn
     case ('B'); sgn = -sgn
     end select

     call calc_cell_bc_flow(bcv, xlength, ylength, zlength, &
        dx, dy, dz, i_w, i_e, j_s, j_n, k_b, k_t)

     select case(bc_plane(bcv))
     case('W','E'); area = dy*dble(j_n-j_s+1)*dz*dble(k_t-k_b+1)
     case('S','N'); area = dx*dble(i_e-i_w+1)*dz*dble(k_t-k_b+1)
     case('B','Y'); area = dx*dble(i_e-i_w+1)*dy*dble(j_n-j_s+1)
     end select

     if(abs(bc_ep_s(bcv,m)) > zero) then
        vel = sgn * bc_volflow_s(bcv,m)/(area*bc_ep_s(bcv,m))
     else
        if(abs(bc_volflow_s(bcv,m)) > zero) then
           vel = zero
        else
           call flush_err_msg(abort=.true.)
        endif
     endif

     if(bc_plane(BCV) == 'W' .or. bc_plane(BCV)== 'E') then
        if(is_defined(bc_u_s(bcv,m))) then
           if(.not.compare(vel, bc_u_s(bcv,m))) then
              write(err_msg, 1300) bcv, (-vel), 'BC_u_s', m, bc_u_s(bcv,m)
              call flush_err_msg(abort=.true.)
           endif
        else
           bc_u_s(bcv,m) = vel
           bc_v_s(bcv,m) = off * bc_v_s(bcv,m)
           bc_w_s(bcv,m) = off * bc_w_s(bcv,m)
        endif

     elseif(bc_plane(BCV) == 'S' .or. bc_plane(BCV)== 'N') then
        if(is_defined(bc_v_s(bcv,m))) then
           if(.not.compare(vel,bc_v_s(bcv,m))) then
              write(err_msg,1300) bcv, vel, 'BC_v_s', m, bc_v_s(bcv,m)
              call flush_err_msg(abort=.true.)
           endif
        else
           bc_v_s(bcv,m) = vel
           bc_u_s(bcv,m) = off * bc_u_s(bcv,m)
           bc_w_s(bcv,m) = off * bc_w_s(bcv,m)
        endif

     elseif(bc_plane(BCV) == 'B' .or. bc_plane(BCV)== 'T') then
        if(is_defined(bc_w_s(bcv,m))) then
           if(.not.compare(vel,bc_w_s(bcv,m))) then
              write(err_msg, 1300) bcv, vel, 'BC_w_s', m, bc_w_s(bcv,m)
              call flush_err_msg(abort=.true.)
           endif
        else
           bc_w_s(bcv,m) = vel
           bc_u_s(bcv,m) = off * bc_u_s(bcv,m)
           bc_v_s(bcv,m) = off * bc_v_s(bcv,m)
        endif
     endif

    call finl_err_msg



1300 format(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed velocity is not equal to specified value',/,&
         ' Value computed from vol. or mass flow  = ',G14.7,/,&
         ' Specified value (',A,I1,') = ',G14.7,/1X,70('*')/)

  end subroutine solids_volflow_to_velocity

end module flow_to_vel_new_module
