module check_solids_common_all_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
  use param1,         only: is_undefined, is_defined
  use run,            only: IFILE_NAME
  use error_manager,  only: finl_err_msg, flush_err_msg, init_err_msg, &
                          & ivar, ival, err_msg

  implicit none
  private 

  public check_solids_common_all



contains
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: CHECK_SOLIDS_COMMON_ALL                                 !
  !  Purpose: Check the solid phase input that is common to all solids   !
  !  phase models.                                                       !
  !                                                                      !
  !    ****** DO NOT PUT MODEL SPECIFIC CHECKS IN THIS ROUTINE ******    !
  !                                                                      !
  !  Use the companion routines for checks specific to a particular      !
  !  solids phase model:                                                 !
  !                                                                      !
  !    > CHECK_SOLIDS_DEM       :: DEM solids phase model                !
  !                                                                      !
  !  Author: J.Musser                                  Date: 03-FEB-14   !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_solids_common_all

    use constant, only: MMAX     ! Number of continuum solids phases.
    use constant, only: D_P0     ! User specified: Initial solids diameter.
    use param,    only: DIM_M    ! Maximum number of solids phases.
    use param1,   only: ZERO     ! Parameter constants

    integer :: M          ! Loop counters.
    integer :: MMAX_L     ! Total number of all solids


    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_COMMON_ALL")

    ! Set the number of solids phases to be checked.
    MMAX_L = MMAX

    ! Check D_p0
    do M = 1, MMAX_L
       if(IS_UNDEFINED(D_P0(M))) then
          write(ERR_MSG, 1000) trim(iVar('D_p0',M)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       elseif(D_P0(M) <= ZERO)then
          write(ERR_MSG, 1001) trim(iVar('D_p0',M)), iVal(D_P0(M)), &
               trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
    enddo

    do M = MMAX_L+1, DIM_M
       if(IS_DEFINED(D_P0(M)))then
          write(ERR_MSG,1002) trim(iVar('D_p0',M)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)
       endif
    enddo

    ! Check solids drag model selection.
    call check_solids_drag

    ! Check the solids density input parameters.
    call check_solids_density(MMAX_L)

    ! Finalize the error messges
    call finl_err_msg


1000 format('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the ',A,' file.')

1001 format('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the ',A,' file.')

1002 format('Error 1002: Illegal input: ',A,' specified out of range.', &
         'Please correct the ',A,' file.')

  end subroutine check_solids_common_all


  !----------------------------------------------------------------------!
  ! Subroutine: CHECK_SOLIDS_DRAG                                        !
  ! Purpose: Check solids species input.                                 !
  !                                                                      !
  ! Author: J. Musser                                  Date: 07-FEB-14   !
  !----------------------------------------------------------------------!
  subroutine check_solids_drag

    use run, only: DRAG_TYPE, DRAG_TYPE_ENUM, SYAM_OBRIEN, GIDASPOW, &
         & GIDASPOW_PCF, GIDASPOW_BLEND, GIDASPOW_BLEND_PCF, WEN_YU, &
         & WEN_YU_PCF, KOCH_HILL, KOCH_HILL_PCF, BVK, HYS, USER_DRAG

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_DRAG")

    select case(trim(adjustl(DRAG_TYPE)))

    case ('SYAM_OBRIEN'); DRAG_TYPE_ENUM = SYAM_OBRIEN
    case ('GIDASPOW'); DRAG_TYPE_ENUM = GIDASPOW
    case ('GIDASPOW_PCF'); DRAG_TYPE_ENUM = GIDASPOW_PCF
    case ('GIDASPOW_BLEND'); DRAG_TYPE_ENUM = GIDASPOW_BLEND
    case ('GIDASPOW_BLEND_PCF'); DRAG_TYPE_ENUM = GIDASPOW_BLEND_PCF
    case ('WEN_YU'); DRAG_TYPE_ENUM = WEN_YU
    case ('WEN_YU_PCF'); DRAG_TYPE_ENUM = WEN_YU_PCF
    CASE ('KOCH_HILL'); DRAG_TYPE_ENUM = KOCH_HILL
    case ('KOCH_HILL_PCF'); DRAG_TYPE_ENUM = KOCH_HILL_PCF
    case ('BVK'); DRAG_TYPE_ENUM = BVK
    case ('HYS'); DRAG_TYPE_ENUM = HYS
    case ('USER_DRAG','USR_DRAG'); DRAG_TYPE_ENUM = USER_DRAG

    case DEFAULT
       write(ERR_MSG,1001)'DRAG_TYPE', trim(adjustl(DRAG_TYPE)), &
            trim(IFILE_NAME)
       call flush_err_msg(ABORT=.true.)
    end select
    
    call finl_err_msg

1001 format('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the ',A,' file.')

  end subroutine check_solids_drag



  !----------------------------------------------------------------------!
  !  Subroutine: CHECK_SOLIDS_DENSITY                                    !
  !  Purpose: check the solid phase density input                        !
  !                                                                      !
  !  Author: J.Musser                                  Date: 03-FEB-14   !
  !----------------------------------------------------------------------!
  subroutine check_solids_density(MMAX_LL)

    use constant, only: RO_s0
    use param,    only: DIM_M

    integer, intent(in) :: MMAX_LL      ! Total number of solids phases
    integer             :: M

    ! Initialize the error manager.
    call init_err_msg("CHECK_SOLIDS_DENSITY")

    ! Check each solids phase.
    do M = 1, MMAX_LL

       ! Verify that one -and only one- solids density model is in use.
       if(IS_UNDEFINED(RO_S0(M))) then
          write(ERR_MSG, 1100) M, trim(IFILE_NAME)
          call flush_err_msg(ABORT=.true.)

1100      format('Error 1101: No solids density information for phase ',  &
               I2,'.',/'Please correct the ',A,' file.')

       endif

    enddo

    ! Check for input overflow.
    do M = MMAX_LL+1, DIM_M
       if(IS_DEFINED(RO_S0(M))) then
          write(ERR_MSG,1002) trim(iVar('RO_s0',M)), trim(IFILE_NAME)
          call flush_err_msg(ABORT=.TRUE.)
       endif
    enddo

    ! Finalize the error messges
    call finl_err_msg


1002 format('Error 1002: Illegal input: ',A,' specified out of range.',&
         'Please correct the ',A,' file.')
    
  end subroutine check_solids_density

end module 
