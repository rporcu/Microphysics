module set_outflow_module

  use bl_fort_module, only : c_real
  use iso_c_binding , only: c_int


  implicit none
  private

  public set_outflow


  !>>> TODO: move definition like the following to a module
  character(16), parameter :: POUTFLOW = 'P_OUTFLOW'


contains


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
  !                                                                      C
  !  Subroutine: SET_OUTFLOW                                             C
  !  Purpose: Set specified outflow bc for pressure outflow,             C
  !  mass outflow, outflow and now also pressure inflow bc               C
  !                                                                      C
  !  Comments:                                                           C
  !  If the outflow boundary is on the W, S or B side of the domain and  C
  !  the component of velocity through the plane is defined then this    C
  !  routine will NOT modify it (i.e., its value from the momentum       C
  !  solver is maintained).                                              C
  !                                                                      C
  !  In general it would seem this routine does little more than what    C
  !  is already done within the momentum routines in terms of velocity.  C
  !  The normal component is either 1) untouched if the outflow is on    C
  !  the W, S, B side of the domain or 2) is set to value of the         C
  !  adjacent fluid cell if it is on the E, N, T side (similar to the    C
  !  respective momentum bc routines). The primary addition here is      C
  !  that the tangential components of a bc cell are set to that of      C
  !  the adjacent fluid cell. Note the tangential components are not     C
  !  explicitly handled in the momentum _BC_ routines; instead their     C
  !  values are based on solution of the momentum equation which is      C
  !  replaced here                                                       C
  !                                                                      C
  !  Several routines are called which perform the following tasks:      C
  !  set_outflow_misc - several derived quantities are set in the        C
  !      boundary                                                        C
  !  set_outflow_ep - the void/volume fraction and bulk densities are    C
  !      set in the boundary                                             C
  !  set_outflow_fluxes - convective fluxes are set in the boundary      C
  !                                                                      C
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  subroutine set_outflow ( bcv, slo, shi, p_g, ep_g, ro_g, rop_g, &
       & u_g, v_g, w_g, flux_ge, flux_gn, flux_gt) 

    use bc,        only: bc_k_b, bc_k_t
    use bc,        only: bc_j_s, bc_j_n
    use bc,        only: bc_i_w, bc_i_e
    use functions, only: iminus, iplus, jminus, jplus, kminus, kplus
    use param1,    only: is_undefined
    use geometry,  only: domlo, domhi

    integer, intent(in) :: slo(3),shi(3)
    integer, intent(in) :: bcv

    real(c_real), intent(inout) :: &
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

    integer      :: i, j, k


    ! loop over the range of boundary cells
    do k = bc_k_b(bcv), bc_k_t(bcv)
       do j = bc_j_s(bcv), bc_j_n(bcv)
          do i = bc_i_w(bcv), bc_i_e(bcv)

             west: if ( i == domlo(1) ) then
                call set_outflow_face ( bcv, slo, shi, i, j, k, iminus(i,j,k), j, k, &
                     & p_g, ro_g, rop_g, ep_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt)
             end if west

             east: if ( i == domhi(1) ) then
                call set_outflow_face ( bcv, slo, shi, i, j, k, iplus(i,j,k), j, k, &
                     & p_g, ro_g, rop_g, ep_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt)
             end if east

             south: if ( j == domlo(2) ) then
                call set_outflow_face ( bcv, slo, shi, i, j, k, i, jminus(i,j,k), k, &
                     & p_g, ro_g, rop_g, ep_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt)
             end if south

             north:  if ( j == domhi(2) ) then
                call set_outflow_face ( bcv, slo, shi, i, j, k, i, jplus(i,j,k), k, &
                     & p_g, ro_g, rop_g, ep_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt)
             end if north

             bottom: if ( k == domlo(3) ) then 
                call set_outflow_face ( bcv, slo, shi, i, j, k, i, j, kminus(i,j,k), &
                     & p_g, ro_g, rop_g, ep_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt)
             end if bottom

             top:  if ( k == domhi(3) ) then
                call set_outflow_face ( bcv, slo, shi, i, j, k, i, j, kplus(i,j,k), &
                     & p_g, ro_g, rop_g, ep_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt)
             end if top

          end do   
       end do   
    end do   

  end subroutine set_outflow



  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !                                                                      
  !  SUBROUTINE:                                                         
  !  set_outflow_face                                                    
  !                                                                      
  !  PURPOSE:                                                            
  !  Set specified outflow bc for pressure outflow,  mass flow,          
  !  outflow and now also pressure inflow bc ona given cell face         
  !                                                                      
  !  COMMENTS:                                                           
  !  This subroutine does two things:                                    
  ! 
  !      1)  Set the value of certain variables in the specified outflow 
  !          boundary cell that would not otherwise be set according     
  !          to their value in the adjacent fluid cell.                  
  !          Previously, this was accomplished via set_outflow_misc()    
  ! 
  !      2)  Set the volume fraction/bulk density (i.e., ep_g) in the 
  !          specified outflow boundary cell that would not otherwise be
  !          set according to their value in the adjacent fluid cell.
  !          Previously, this was accomplished via set_outflow_eps()    
  ! 
  !      3)  Set the boundary cell value of the normal component of 
  !          velocity according to the value in the adjacent fluid cell.
  !          Note the value of the boundary velocity is a scaled version 
  !          of the value of the adjacent fluid velocity based on the 
  !          concentration ratio of the fluid cell to the boundary cell.
  !          For the gas phase, this ratio is most likely 1 except for
  !          compressible cases with a po/pi boundary where p_g of the 
  !          boundary is set and may differ from the value of the adjacent 
  !          fluid cell.
  !          For the solids phase this seems unnecessary..? 
  !          differences may arise  
  ! 
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine set_outflow_face ( bcv, slo, shi, ic, jc, kc, in, jn, kn, &
       & p_g, ro_g, rop_g, ep_g, u_g, v_g, w_g, flux_ge, flux_gn, flux_gt)

    use bc       , only: bc_type, bc_ep_g
    use fld_const, only: ro_g0, mw_avg
    use eos      , only: eosg
    use param1,    only: is_undefined, zero, one

    integer, intent(in)      :: slo(3), shi(3)
    integer, intent(in)      :: bcv         ! boundary condition number
    integer, intent(in)      :: ic, jc, kc  ! Cell     indeces
    integer, intent(in)      :: in, jn, kn  ! Neighbor indeces

    real(c_real), intent(inout)    ::                              &
         &     p_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &    ro_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &   rop_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &    ep_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &     u_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &     v_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         &     w_g( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         & flux_ge( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         & flux_gn( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) ), &
         & flux_gt( slo(1):shi(1), slo(2):shi(2), slo(3):shi(3) )

    real(c_real)        :: sum_eps     ! sum of solids phases volume fractions
    real(c_real)        :: sum_rops    ! sum of solids phases bulk densities


    if ( bc_type(bcv) /= POUTFLOW ) p_g(ic,jc,kc) = p_g(in,jn,kn)

    ! >>> TODO:  what is 295.15D0? Better define constants as parameters
    if (is_undefined(ro_g0)) ro_g(ic,jc,kc) = eosg(mw_avg,p_g(in,jn,kn),295.15d0)


    sum_rops = zero
    sum_eps = zero

    ! if bc_ep_g undefined, set ep_g accordingly (based on flow condition
    ! or based on bc_rop_s). if bc_ep_g is defined its set value will be
    ! maintained (from set_bc0).
    if (is_undefined(bc_ep_g(bcv))) ep_g(ic,jc,kc) = one 

    ! now that ep_g in the boundary cell is known, define the bulk density
    ! of the gas phase in the boundary cell
    rop_g(ic,jc,kc) = ro_g(ic,jc,kc)*ep_g(ic,jc,kc)

    ! provide an initial value for the velocity component through the domain
    ! otherwise its present value (from solution of the corresponding
    ! momentum eqn) is kept. values for the velocity components in the off
    ! directions are modified (needed for po or o boundaries but not mo or
    ! pi as velocities should be fully specified by this point)
    ! >>> IMPORTANT: the following check ( is_defined ) was not performed    
    !                for the west boundary in the original code: why???
    if ( is_undefined( w_g(ic,jc,kc) ) ) then

       ! the tangential components are not explicitly handled in the boundary
       ! condition routines of the corresponding momentum equation
       if (rop_g(ic,jc,kc) > zero) then
          w_g(ic,jc,kc) = rop_g(in,jn,kn) * w_g(in,jn,kn) / rop_g(ic,jc,kc)
       else
          w_g(ic,jc,kc) = zero
       end if

    end if

    u_g(ic,jc,kc) = u_g(in,jn,kn)
    v_g(ic,jc,kc) = v_g(in,jn,kn)

    flux_ge(ic,jc,kc) = flux_ge(in,jn,kn)
    flux_gn(ic,jc,kc) = flux_gn(in,jn,kn)
    flux_gt(ic,jc,kc) = flux_gt(in,jn,kn)


  end subroutine set_outflow_face


end module set_outflow_module
