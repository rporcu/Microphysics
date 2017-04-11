module check_initial_conditions_module

  use amrex_fort_module, only: c_real => amrex_real
  use iso_c_binding ,    only: c_int

  use param1,  only: zero, one, undefined, undefined_i
  use param1,  only: is_defined, is_undefined

  use error_manager, only: init_err_msg, finl_err_msg, flush_err_msg
  use error_manager, only: err_msg, ivar, ival

  implicit none
  private

  public check_initial_conditions

contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_INITIAL_CONDITIONS                                !
!                                                                      !
!  Purpose: check the initial conditions input section                 !
!     - check geometry of any specified IC region                      !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_initial_conditions(dx,dy,dz,domlo,domhi)

    use ic,                    only: ic_defined
    use run,                   only: dem_solids, run_type
    use param,                 only: dim_ic

    integer(c_int), intent(in) :: domlo(3),domhi(3)
    real(c_real)  , intent(in) :: dx, dy, dz
    integer(c_int)             :: icv

    ! Determine which ICs are DEFINED
    call check_ic_geometry(dx,dy,dz,domlo,domhi)

    ! Loop over all IC arrays.
    do icv=1, dim_ic

       ! Verify user input for defined defined IC.
       if(ic_defined(icv)) then
          ! Gas phase checks.
          call check_ic_gas_phase(ICV)
          ! Generic solids phase checks.
          call check_ic_solids_phases(ICV)

       ! Verify that no data was defined for unspecified IC. ICs are only
       ! defined for new runs, so these checks are restricted to new runs.
       elseif(run_type == 'NEW') then
          call check_ic_overflow(ICV)
       endif
    enddo

    ! Check the initial conditions for the DEM model as well
    if(dem_solids) call check_ic_common_discrete

 end subroutine check_initial_conditions

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_GEOMETRY                                        !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    subroutine check_ic_geometry(dx,dy,dz,domlo,domhi)

      use param, only: dim_ic
      use ic,    only: ic_defined
      use ic,    only: IC_X_e, IC_Y_n, IC_Z_t
      use ic,    only: IC_X_w, IC_Y_s, IC_Z_b

      use calc_cell_module, only: calc_cell_ic

      integer(c_int), intent(in) :: domlo(3),domhi(3)
      real(c_real)  , intent(in) :: dx,dy,dz

      integer :: icv, i_w, i_e, j_s, j_n, k_b, k_t

      ! Initialize the error manager.
      call init_err_msg("CHECK_IC_GEOMETRY")

      ! Check geometry of any specified IC region
      do icv = 1, dim_ic

         if(ic_defined(icv)) then

            if (is_undefined(ic_x_w(icv))) then
               write(err_msg, 1100) icv, 'IC_X_w'
               call flush_err_msg(abort=.true.)
            endif

           if (is_undefined(ic_x_e(icv))) then
              write(err_msg, 1100) icv, 'IC_X_e'
              call flush_err_msg(abort=.true.)
           endif

           if (is_undefined(ic_y_s(icv))) then
              write(err_msg, 1100) icv, 'IC_Y_s'
              call flush_err_msg(abort=.true.)
           endif

           if (is_undefined(ic_y_n(icv))) then
              write(err_msg, 1100) icv, 'IC_Y_n'
              call flush_err_msg(abort=.true.)
           endif

           if (is_undefined(ic_z_b(icv))) then
              write(err_msg, 1100) icv, 'IC_Z_b'
              call flush_err_msg(abort=.true.)
           endif

           if (is_undefined(ic_z_t(icv))) then
              write(err_msg, 1100) icv, 'IC_Z_t'
              call flush_err_msg(abort=.true.)
           endif

1100 FORMAT('Error 1100: Initial condition region ',I3,' is ill-',    &
        'defined.',/' > ',A,' is not specified.',/'Please correct ', &
        'the input deck.')

           call calc_cell_ic(dx, dy, dz, &
                  ic_x_w(icv), ic_y_s(icv), ic_z_b(icv), &
                  ic_x_e(icv), ic_y_n(icv), ic_z_t(icv), &
                  i_w, i_e, j_s, j_n, k_b, k_t)

           ! Report problems with calculated bounds.
           if(i_w > i_e) then
              write(err_msg, 1101) icv, 'IC_I_W > IC_I_E'
              call flush_err_msg(abort=.true.)
           elseif(i_w < domlo(1)) then
              write(err_msg, 1101) icv, 'IC_I_W < domlo(1)'
              call flush_err_msg(abort=.true.)
           elseif(i_w > domhi(1)) then
              write(err_msg, 1101) icv, 'IC_I_W > domhi(1)'
              call flush_err_msg(abort=.true.)
           elseif(i_e < domlo(1)) then
              write(err_msg, 1101) icv, 'IC_I_E < domlo(1)'
              call flush_err_msg(abort=.true.)
           elseif(i_e > domhi(1)) then
              write(err_msg, 1101) icv, 'IC_I_E > domhi(1)'
              call flush_err_msg(abort=.true.)
           endif

           if(j_s > j_n) then
              write(err_msg, 1101) icv, 'IC_J_S > IC_J_N'
              call flush_err_msg(abort=.true.)
           elseif(j_s<domlo(2)) then
              write(err_msg, 1101) icv, 'IC_J_S < domlo(2)'
              call flush_err_msg(abort=.true.)
           elseif(j_s>domhi(2)) then
              write(err_msg, 1101) icv, 'IC_J_S >  domhi(2)'
              call flush_err_msg(abort=.true.)
           elseif(j_n<domlo(2)) then
              write(err_msg, 1101) icv, 'IC_J_N < domlo(2)'
              call flush_err_msg(abort=.true.)
           elseif(j_n>domhi(2)) then
              write(err_msg, 1101) icv, 'IC_J_N > domhi(2)'
              call flush_err_msg(abort=.true.)
           endif

           if(k_b > k_t) then
              write(err_msg, 1101) icv, 'IC_K_B > IC_K_T'
              call flush_err_msg(abort=.true.)
           elseif(k_b < domlo(3)) then
              write(err_msg, 1101) icv, 'IC_K_B < domlo(3)'
              call flush_err_msg(abort=.true.)
           elseif(k_b > domhi(3)) then
              write(err_msg, 1101) icv, 'IC_K_B > domhi(3)'
              call flush_err_msg(abort=.true.)
           elseif(k_t < domlo(3)) then
              write(err_msg, 1101) icv, 'IC_K_T < domlo(3)'
              call flush_err_msg(abort=.true.)
           elseif(k_t > domhi(3)) then
              write(err_msg, 1101) icv, 'IC_K_T > domhi(3)'
              call flush_err_msg(abort=.true.)
           endif

        endif
     enddo   ! end loop over (icv=1,dim_ic)

1101 format('Error 1101: Initial condition region ',I2,' is ill-',    &
        'defined.',/3x,A,/'Please correct the input deck.')

      call finl_err_msg

   end subroutine check_ic_geometry


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_GAS_PHASE                                       !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Verify gas phase input variables in IC region.              !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   subroutine check_ic_gas_phase(ICV)

      use fld_const, only: ro_g0
      use ic,        only: IC_P_g, IC_U_g, IC_V_g, IC_W_g,  IC_TYPE

      integer, intent(in) :: ICV
      logical             :: BASIC_IC

      call init_err_msg("CHECK_IC_GAS_PHASE")

      ! Patch ICs skip various checks.
      basic_ic = (ic_type(icv) /= 'PATCH')

      ! Check that gas phase velocity components are initialized.
      if(basic_ic) then
         if(is_undefined(ic_u_g(icv))) then
            write(err_msg, 1000) trim(ivar('IC_U_g',icv))
            call flush_err_msg(abort=.true.)
         endif

         if(is_undefined(ic_v_g(icv))) then
            write(err_msg, 1000) trim(ivar('IC_V_g',icv))
            call flush_err_msg(abort=.true.)
         endif

         if(is_undefined(ic_w_g(icv))) then
            write(err_msg, 1000) trim(ivar('IC_W_g',icv))
            call flush_err_msg(abort=.true.)
         endif
      endif

1000  format('Error 1000: Required input not specified: ',A,/&
         'Please correct the input deck.')

      ! Check that if the gas phase pressure is initialized and the gas is
      ! compressible that the gas phase pressure is not zero or negative
      if(is_defined(ic_p_g(icv))) then
         if(is_undefined(ro_g0).and. ic_p_g(icv)<=zero) then
            write(err_msg, 1100) trim(iVar('IC_P_g',icv)),             &
                 ival(ic_p_g(icv))
            call flush_err_msg(abort=.true.)
         endif
      endif

1100  format('Error 1100: Pressure must be greater than 0.0 for ',     &
         'compressible flow',/'Illegal value: ',A,' = ',A,/&
         'Please correct the input deck.')

      call finl_err_msg

   end subroutine check_ic_gas_phase


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_IC_SOLIDS_PHASES                                  !
!  Author: P. Nicoletti                               Date: 02-DEC-91  !
!  Author: J.Musser                                   Date: 01-MAR-14  !
!                                                                      !
!  Purpose: Verify solids phase(s) input variables in IC region.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
    subroutine check_ic_solids_phases(ICV)

      use param1, only: equal
      use constant, only: mmax
      use ic,       only: ic_ep_s, ic_u_s, ic_v_s, ic_w_s
      use ic,       only: ic_ep_g, ic_type

      integer, intent(in) :: icv
      integer             :: m
      real(c_real)        :: sum_ep
      logical             :: basic_ic

      ! Initialize the error manager.
      CALL init_err_msg("CHECK_IC_SOLIDS_PHASES")

      ! Patch ICs skip various checks.
      basic_ic = (ic_type(icv) /= 'PATCH')

      ! Initialize the sum of the total volume fraction.
      sum_ep = ic_ep_g(icv)

      ! check that solids phase-m components are initialized
      if(basic_ic) then
         do m=1, mmax

            if(ic_ep_s(icv,m) > epsilon(0.d0)) then

               if(is_undefined(ic_u_s(icv,m))) then
                  write(err_msg, 1000)trim(ivar('IC_U_s',icv,m))
                  call flush_err_msg(abort=.true.)
               endif

               if(is_undefined(ic_v_s(icv,m))) then
                  write(err_msg, 1000)trim(ivar('IC_V_s',icv,m))
                  call flush_err_msg(abort=.true.)
               endif

               if(is_undefined(ic_w_s(icv,m))) then
                  write(err_msg, 1000)trim(ivar('IC_W_s',icv,m))
                  call flush_err_msg(abort=.true.)
               endif

               sum_ep = sum_ep + ic_ep_s(icv,m)

            endif
         enddo
      endif

1000  format('Error 1000: Insufficient solids phase ',I2,' ',          &
           'information for IC',/'region ',I3,'. ',A,' not specified.',/ &
           'Please correct the input deck.')

      ! Verify that the volume fractions sum to one.
      if(basic_ic .and. .not.equal(sum_ep,one)) then
         write(err_msg,1410) icv
         call flush_err_msg(abort=.true.)
      endif

1410  format('Error 1410: Illegal initial condition region : ',I3,/    &
           'Sum of volume fractions does NOT equal ONE. Please correct',/&
           'the input deck.')

      call finl_err_msg

    end subroutine check_ic_solids_phases


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_OVERFLOW                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Verify that no data was defined for unspecified IC.         !
!                                                                      !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    subroutine check_ic_overflow(ICV)

      use param, only: dim_m
      use ic,    only: ic_type
      use ic,    only: ic_ep_g, ic_u_g, ic_v_g, ic_w_g
      use ic,    only: ic_ep_s, ic_u_s, ic_v_s, ic_w_s

      integer, intent(in) :: icv
      integer :: m

      if (ic_type(icv) == 'PATCH') return

      ! Initialize the error manager.
      call init_err_msg("CHECK_IC_OVERFLOW")

      ! gas phase quantities
      if(is_defined(ic_u_g(icv))) then
         write(err_msg, 1010) trim(ivar('IC_U_g',ICV))
         call flush_err_msg(abort=.true.)
      elseif(is_defined(ic_v_g(icv))) then
         write(err_msg, 1010) trim(ivar('IC_V_g',ICV))
         call flush_err_msg(abort=.true.)
      elseif(is_defined(ic_w_g(icv))) then
         write(err_msg, 1010) trim(ivar('IC_W_g',ICV))
         call flush_err_msg(abort=.true.)
      elseif(abs(ic_ep_g(icv)-one) > epsilon(0.)) then
         write(err_msg, 1010) trim(ivar('IC_EP_g',ICV))
         call flush_err_msg(abort=.true.)
      endif

      ! solids phase quantities
      do m=1, dim_m
         if(ic_ep_s(icv,m) > epsilon(0.0d0)) then
            write(err_msg, 1010) trim(ivar('IC_ROP_s',ICV,M))
            call flush_err_msg(abort=.true.)
         elseif(is_defined(ic_u_s(icv,m))) then
            write(err_msg, 1010) trim(ivar('IC_U_s',ICV,M))
            call flush_err_msg(abort=.true.)
         elseif(is_defined(ic_v_s(icv,m))) then
            write(err_msg, 1010) trim(ivar('IC_V_s',ICV,M))
            call flush_err_msg(abort=.true.)
         elseif(is_defined(ic_w_s(icv,m))) then
            write(err_msg, 1010) trim(ivar('IC_W_s',ICV,M))
            call flush_err_msg(abort=.true.)
         endif
      enddo

      call finl_err_msg

1010  format('Error 1010: ',A,' specified in an undefined IC region')

    end subroutine check_ic_overflow


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_IC_COMMON_DISCRETE                                !
!                                                                      !
!  Purpose: check the initial conditions input section for DEM.        !
!     - ensure the ICs with solids are non-overlapping                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine check_ic_common_discrete

      use param1,   only: zero
      use param,    only: dim_ic
      use ic,       only: ic_defined, ic_ep_s
      use ic,       only: ic_x_w, ic_x_e, ic_y_s, ic_y_n, ic_z_b, ic_z_t

      use error_manager,  only: init_err_msg, finl_err_msg
      use error_manager,  only: err_msg, flush_err_msg

      integer      :: icv1, icv2

      ! Initialize the error manager.
      call init_err_msg("CHECK_IC_COMMON_DISCRETE")

      ! Check if the ICs are non-overlapping.
      do icv1 = 1, dim_ic
         if(ic_defined(icv1) .and. sum(ic_ep_s(icv1,:)) > zero) then
            do icv2 = icv1+1, dim_ic
               if(ic_defined(icv2) .and. sum(ic_ep_s(icv1,:)) > zero) then
                  if((ic_x_w(icv1) < ic_x_e(icv2) .and. ic_x_e(icv1) > ic_x_w(icv2)) .or. &
                     (ic_y_s(icv1) < ic_y_n(icv2) .and. ic_y_n(icv1) > ic_y_s(icv2)) .or. &
                     (ic_z_b(icv1) < ic_z_t(icv2) .and. ic_z_t(icv1) > ic_z_b(icv2))) then
                     write(err_msg, 1004) icv1, icv2
                     call flush_err_msg(abort=.true.)
                  endif
               endif
            enddo
         endif
      enddo

1004 format('Error 1004: Overlapping IC regions with nonzero solids ',&
        'volume',/'fraction detected. This is not supported for ',    &
        'discrete solids.',2/'Overlapping ICs: ',2(2x,I4),2/,         &
        'Please correct the input deck.')


      call finl_err_msg

   end subroutine check_ic_common_discrete

end module check_initial_conditions_module
