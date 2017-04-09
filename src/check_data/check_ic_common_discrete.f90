module check_ic_common_discrete_module

  implicit none
  private

  public check_ic_common_discrete

contains

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
      use param,    only: dimension_ic
      use ic,       only: ic_defined, ic_ep_s
      use ic,       only: ic_x_w, ic_x_e, ic_y_s, ic_y_n, ic_z_b, ic_z_t

      use error_manager,  only: init_err_msg, finl_err_msg
      use error_manager,  only: err_msg, flush_err_msg

      integer      :: icv1, icv2

      ! Initialize the error manager.
      call init_err_msg("CHECK_IC_COMMON_DISCRETE")

      ! Check if the ICs are non-overlapping.
      do icv1 = 1, dimension_ic
         if(ic_defined(icv1) .and. sum(ic_ep_s(icv1,:)) > zero) then
            do icv2 = icv1+1, dimension_ic
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

end module check_ic_common_discrete_module
