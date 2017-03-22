module check_ic_common_discrete_module


  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use run,            only: IFILE_NAME
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg, &
       & ivar, ival, err_msg

  implicit none
  private

  public check_ic_common_discrete

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: CHECK_IC_COMMON_DISCRETE                                !
  !  Author:   R.Garg                                   Date: 11-Mar-14  !
  !                                                                      !
  !  Purpose: check the initial conditions input section for DEM.        !
  !     - ensure the first IC is defined over the entire domain with     !
  !        ep_g = 1 when more than one IC has solids                     !
  !     - ensure the ICs are non-overlapping                             !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_ic_common_discrete

 
    use param1,   only: undefined_i, zero, one, equal
    use geometry, only: xlength, ylength, zlength
    use param,    only: dimension_ic
    use ic,       only: ic_defined, ic_ep_g, ic_x_w, ic_x_e, ic_y_s, &
                      &  ic_y_n, ic_z_b, ic_z_t

    integer      :: icv, icv2, idim, count_ic, count_ic_with_sols, first_def_ic
    real(c_real) :: ic_orig(3), ic_end(3), ic2_orig(3) , ic2_end(3)
    real(c_real) :: ic_min, ic_max, ic2_min, ic2_max , tol_ic_reg
    logical      :: sep_axis, first_ic_ok



    ! Initialize the error manager.
    call init_err_msg("CHECK_IC_COMMON_DISCRETE")

    ! First check if multiple IC regions are defined for non-zero solids volume
    ! fraction, then check if the first IC is specified over the whole domain with IC_EP_g = 1
    count_ic           = 0     !total count of defined ics
    count_ic_with_sols = 0     !total count of defined ic's with solids
    first_def_ic       = undefined_i
    
    do icv = 1, dimension_ic
       
       if (ic_defined(icv)) then
          count_ic = count_ic + 1
          first_def_ic = min(first_def_ic, icv)
          
          if(ic_ep_g(icv).lt.one) count_ic_with_sols &
               = count_ic_with_sols  + 1
          
       endif
    end do
    
    if(count_ic_with_sols >= 1 .and. &
         count_ic > count_ic_with_sols+1) then
       
       ! if the number of ic's with solids is greater than one, make sure the
       ! first ic spans the entire domain with voidage of one. this ensures
       ! that the entire domain has valid ics defined.
       icv = first_def_ic
       first_ic_ok = .false.
       if(equal(ic_ep_g(icv), one) &
            .and.ic_x_w(icv).le.zero.and.ic_x_e(icv).ge.xlength         &
            .and.ic_y_s(icv).le.zero.and.ic_y_n(icv).ge.ylength)        &
            first_ic_ok = .true.
       
       if (first_ic_ok .and. ic_z_b(icv) <= zero .and. &
            ic_z_t(icv) >= zlength) first_ic_ok = .true.
       
       if(.not.first_ic_ok) then
          write(err_msg, 1003) trim(IFILE_NAME)
          call flush_err_msg(abort=.true.)
       endif
       
1003   format(' Error 1003: Particle seeding with more than one IC ',   &
            'region requires',/'that IC 1 span the entire domain and ',   &
            'have IC_EP_g(1) = 1.0.',/'Please correct the ',A,' file.')
       
    endif
    
    ! Check if the ICs are non-overlapping.
    tol_ic_reg  = 1e-04
    icvloop : do icv = 1, dimension_ic
       
       if(.not.ic_defined(icv)) cycle icvloop
       if(equal(ic_ep_g(icv), 1.d0)) cycle icvloop
       ic_orig(1) = ic_x_w(icv)
       ic_orig(2) = ic_y_s(icv)
       ic_orig(3) = ic_z_b(icv)
       ic_end(1)  = ic_x_e(icv)
       ic_end(2)  = ic_y_n(icv)
       ic_end(3)  = ic_z_t(icv)
       icvtwoloop : do icv2 = icv+1, dimension_ic
          
          if(.not.ic_defined(icv2)) cycle icvtwoloop
          if(equal(ic_ep_g(icv2), 1.0d0)) cycle icvtwoloop
          
          ic2_orig(1) = ic_x_w(icv2)
          ic2_orig(2) = ic_y_s(icv2)
          ic2_orig(3) = ic_z_b(icv2)
          ic2_end(1)  = ic_x_e(icv2)
          ic2_end(2)  = ic_y_n(icv2)
          ic2_end(3)  = ic_z_t(icv2)
          
          sep_axis  = .false.
          do idim = 1, 3
             
             ic_min = ic_orig(idim)
             ic_max = ic_end(idim)
             ic2_min = ic2_orig(idim)
             ic2_max = ic2_end(idim)
             
             ! Check for separating axis. If the separating axis exists, then the IC
             ! regions can't overlap generally equality implies lack of sep_axis,
             ! and thus, overlapping. However, doing so will flag all IC's as
             ! overlapping since IC's have to share common edges. So here the
             ! equality is considered as existence of a separating axis, and hence,
             ! no overlap equality is also considered as separating axis which is
             if ((ic_min .ge. ic2_max)  .or. (ic_max .le. ic2_min) ) then
                sep_axis = .true.
                exit
             endif
          end do
          
          ! Implies the IC regions could not find a separating axis and are
          ! therefore overlapping.
          if(.not.sep_axis) then
             write(err_msg, 1004) icv, icv2, trim(ifile_name)
             call flush_err_msg(abort=.true.)
          endif
          
1004      format('Error 1004: Overlapping IC regions with nonzero solids ',&
               'volume',/'fraction detected. This is not supported for ',    &
               'discrete solids.',2/'Overlapping ICs: ',2(2x,I4),2/,         &
               'Please correct the ',A,' file.')
          
       end do icvtwoloop
    end do icvloop
    
    
    
    call finl_err_msg
    
  end subroutine check_ic_common_discrete

end module check_ic_common_discrete_module
