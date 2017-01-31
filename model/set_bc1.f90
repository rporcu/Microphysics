module set_bc1_module

  use bl_fort_module, only: c_real
  use iso_c_binding , only: c_int

  implicit none 
  private

  public set_bc1


contains


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !                                                                      
  !  SUBROUTINE:
  !  set_bc1() 
  !                                                
  !  PURPOSE:
  !  set transient flow boundary conditions
  ! 
  !  AUTHOR:                                                                     
  !  M. Syamlal 
  ! 
  !  DATE:
  !  29 January 1992 
  !                                                                      
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine set_bc1 ( time, dt, slo, shi, p_g, ep_g, ro_g, rop_g, u_g, &
       & v_g, w_g, flux_ge, flux_gn, flux_gt, dx, dy, dz) &
       & bind( c, name="set_bc1" )

    use bc,                 only: bc_defined, bc_type
    use param ,             only: dimension_bc
    use set_outflow_module, only: set_outflow


    integer(c_int), intent(in)    :: slo(3),shi(3)
    real(c_real),   intent(in)    :: dt, time, dx, dy, dz

    real(c_real),   intent(inout) :: &
         &     p_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &    ep_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &    ro_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &   rop_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &     u_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &     v_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &     w_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         & flux_ge( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         & flux_gn( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         & flux_gt( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) )

    integer :: l

    ! set the boundary conditions
    do l = 1, dimension_bc

       if ( .not. bc_defined(l) )                  cycle 

       ! Do not anything for mass_inflow BC
       if ( trim( bc_type(l) ) == 'mass_inflow' )  cycle                        

       ! All transient BCs need this except for mass inflow case.
       ! P_inflow needs ONLY this
       call set_outflow( l, slo, shi, p_g, ep_g, ro_g, rop_g, &
            & u_g, v_g, w_g, flux_ge, flux_gn, flux_gt ) 

       select case ( trim( bc_type(l) ) )

       case ('p_outflow')
          call set_bc1_report_outflow(l, time, dt, slo, shi, &
               & u_g, v_g, w_g, rop_g, ep_g, dx, dy, dz )

       case ('mass_outflow')
          call set_bc1_adjust_outflow(l, time, dt, slo, shi, &
               u_g,v_g,w_g,rop_g,ep_g,dx,dy,dz)

       case ('outflow')
          call set_bc1_report_outflow(l,time, dt, slo, shi, &
               u_g,v_g,w_g,rop_g,ep_g,dx,dy,dz)

       case default
          ! Do nothing for p_inflow condition
       end select

    end do

  end subroutine set_bc1


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !                                                                      
  !  SUBROUTINE:
  !  set_bc1_report_outflow() 
  !                                                
  !  PURPOSE:
  !  print out outflow conditions
  !                                                                      
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine set_bc1_report_outflow( bcv, time, dt, slo, shi, &
       &  u_g, v_g, w_g, rop_g, ep_g, dx, dy, dz )

    use bc,                  only: bc_dt_0, bc_time
    use bc,                  only: bc_mout_g
    use bc,                  only: bc_out_n
    use bc,                  only: bc_vout_g
    use calc_outflow_module, only: calc_outflow
    use funits,              only: dmp_log, unit_log
    use param1,              only: is_undefined, zero
    use run,                 only: tstop

    integer,      intent(in) :: bcv, slo(3), shi(3)
    real(c_real), intent(in) :: dt, time, dx, dy, dz

    real(c_real),   intent(inout) :: &
         &     u_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &     v_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &     w_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &   rop_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &    ep_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) )

    character(len=*), parameter :: &
         & str1 = '(/,1x,"average outflow rates at bc no. ",i2,"  at time = ",g12.5)', &
         & str2 = '(3x,"gas : mass flow = ",g12.5,"     volumetric flow = ",g12.5)'


    if ( is_undefined( bc_dt_0(bcv) ) ) return

    call calc_outflow( bcv, slo, shi, u_g, v_g, w_g, rop_g, ep_g, dx, dy, dz )

    ! calculate and accumulate the actual mass and volume outflow
    if ( ( time + 0.1d0 * dt >= bc_time(bcv) ) .or. &
         & ( time + 0.1d0 * dt >= tstop ) ) then

       bc_time(bcv) = time + bc_dt_0(bcv)

       ! average and print out the flow rates
       bc_mout_g(bcv) = abs( bc_mout_g(bcv) ) / bc_out_n(bcv)
       bc_vout_g(bcv) = abs( bc_vout_g(bcv) ) / bc_out_n(bcv)

       if (dmp_log) write (unit_log, str1) bcv, time
       if (dmp_log) write (unit_log, str2) bc_mout_g(bcv), bc_vout_g(bcv)

       bc_mout_g(bcv) = zero
       bc_vout_g(bcv) = zero
       bc_out_n(bcv)  = 0

    end if


  end subroutine set_bc1_report_outflow



  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !                                                                      
  !  SUBROUTINE:
  !  set_bc1_report_outflow() 
  !                                                
  !  PURPOSE:
  !  adjust velocities to get specified mass or volumetric flow rate
  !  on average outflow rate.
  !                                                                      
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine set_bc1_adjust_outflow(bcv, time, dt, slo, shi, &
       u_g, v_g, w_g, rop_g, ep_g, dx, dy, dz)

    use calc_outflow_module, only: calc_outflow
    use functions,           only: iminus, jminus, kminus
    use funits,              only: dmp_log, unit_log
    use param1,              only: is_defined, zero, small_number
    use run,                 only: tstop
    use bc,                  only:                                   & 
         & bc_dt_0, bc_time, bc_i_w, bc_i_e,                         &
         & bc_j_s, bc_j_n, bc_k_b, bc_k_t,  bc_massflow_g,           &
         & bc_mout_g, bc_out_n, bc_plane, bc_u_g, bc_v_g, bc_w_g,    &
         & bc_volflow_g, bc_vout_g


    integer,      intent(in)    :: bcv, slo(3), shi(3)
    real(c_real), intent(in)    :: dt, time, dx, dy, dz

    real(c_real), intent(inout) :: &
         &     u_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &     v_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &     w_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &   rop_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &    ep_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) )

    integer  :: i, j, k

    character(len=*), parameter :: &
         & str1 = '(/,1x,"average outflow rates at bc no. ",i2,"  at time = ",g12.5)', &
         & str2 = '(3x,"gas : mass flow = ",g12.5,"     volumetric flow = ",g12.5)'


    call calc_outflow( bcv, slo, shi, u_g, v_g, w_g, rop_g, ep_g, dx, dy, dz )

    ! calculate and accumulate the actual mass and volume outflow
    time_if: if (    ( time + 0.1d0 * dt >= bc_time(bcv) ) .or. &
         &  ( time + 0.1d0 * dt >= tstop )         ) then

       bc_time(bcv) = time + bc_dt_0(bcv)

       ! average and print out the flow rates
       bc_mout_g(bcv) = abs( bc_mout_g(bcv) ) / bc_out_n(bcv)
       bc_vout_g(bcv) = abs( bc_vout_g(bcv) ) / bc_out_n(bcv)

       if (dmp_log)  write(unit_log, str1) bcv, time
       if (dmp_log)  write(unit_log, str2) bc_mout_g(bcv), bc_vout_g(bcv)

       bc_out_n(bcv) = 0

       ! now that we know the mass and volume outflow update the bc velocities
       ! (gas phase)
       bc_if: if ( is_defined( bc_massflow_g(bcv) ) ) then

          mass_outflow: if ( bc_mout_g(bcv) > small_number ) then

             select case ( trim( bc_plane(bcv) ) )
             case ('w', 'e')
                bc_u_g(bcv) = bc_u_g(bcv) * bc_massflow_g(bcv) / bc_mout_g(bcv)
             case ('s', 'n')
                bc_v_g(bcv) = bc_v_g(bcv) * bc_massflow_g(bcv) / bc_mout_g(bcv)
             case ('b', 't')
                bc_w_g(bcv) = bc_w_g(bcv) * bc_massflow_g(bcv) / bc_mout_g(bcv)
             end select

          end if mass_outflow

       else if ( is_defined( bc_volflow_g(bcv) ) ) then 

          vol_outflow: if (bc_vout_g(bcv) > small_number) then
             select case  ( trim( bc_plane(bcv) ) )
             case ('w', 'e')
                bc_u_g(bcv) = bc_u_g(bcv) * bc_volflow_g(bcv) / bc_vout_g(bcv)
             case ('s', 'n')
                bc_v_g(bcv) = bc_v_g(bcv) * bc_volflow_g(bcv) / bc_vout_g(bcv)
             case ('b', 't')
                bc_w_g(bcv) = bc_w_g(bcv) * bc_volflow_g(bcv) / bc_vout_g(bcv)
             end select
          end if vol_outflow

       end if bc_if


       ! zero out counter for new cycle
       bc_mout_g(bcv) = zero
       bc_vout_g(bcv) = zero

       ! apply updated boundary velocities - define the field variables at the
       ! boundaries according to user specifications with modifications from
       ! the above calculations.
       ! if the boundary plane is w, s, or b (i.e., the fluid cell is on the
       ! west, south or bottom of the boundary cell) then define the velocity
       ! of the adjacent fluid cell according to the boundary velocity rather
       ! than the velocity of the boundary cell.
       ! why not set the velocity in the boundary cell itself?  based on the
       ! momentum bc routine it should not really matter for mo.
       do k = bc_k_b(bcv), bc_k_t(bcv)
          do j = bc_j_s(bcv), bc_j_n(bcv)
             do i = bc_i_w(bcv), bc_i_e(bcv)

                select case ( trim( bc_plane(bcv) ) )
                case ('w'); u_g(iminus(i,j,k),j,k) = bc_u_g(bcv)
                case ('e'); u_g(i,j,k)             = bc_u_g(bcv)
                case ('s'); v_g(i,jminus(i,j,k),k) = bc_v_g(bcv)
                case ('n'); v_g(i,j,k)             = bc_v_g(bcv)
                case ('b'); w_g(i,j,kminus(i,j,k)) = bc_w_g(bcv)
                case ('t'); w_g(i,j,k)             = bc_w_g(bcv)
                end select

             end do
          end do
       end do

    end if  time_if

  end subroutine set_bc1_adjust_outflow

end module set_bc1_module
