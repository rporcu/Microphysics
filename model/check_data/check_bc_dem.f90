module check_bc_dem_module

  use bl_fort_module, only: c_real
  use iso_c_binding , only: c_int
  use run,            only: IFILE_NAME
  use error_manager,  only: finl_err_msg, err_msg, flush_err_msg, &
                          & init_err_msg, ivar


  implicit none
  private

  public check_bc_dem


contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  ! minimum amount of geometry data.                                     !
  !                                                                      !
  ! Subroutine: CHECK_BC_DEM                                             !
  ! Author: J.Musser                                    Date: 01-Mar-14  !
  !                                                                      !
  ! Purpose: Determine if BCs are "DEFINED" and that they contain the    !
  ! minimum amount of geometry data.                                     !
  !                                                                      !
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  subroutine check_bc_dem(M_TOT)

    use bc,     only: BC_TYPE, BC_EP_s, BC_PO_APPLY_TO_DES
    use run,    only: SOLIDS_MODEL
    use des_bc, only: DEM_BCMI, DEM_BCMO, DEM_BCMI_MAP, DEM_BCMO_MAP
    use param,  only: DIMENSION_BC
    use param1, only: ZERO, IS_DEFINED

    integer, intent(in) :: M_TOT
    integer             :: BCV, M

    ! Initialize the error manager.
    call init_err_msg("CHECK_BC_DEM")

    ! Initialize
    DEM_BCMI = 0
    DEM_BCMO = 0

    ! Loop over all BCs looking for DEM solids inlets/outlets
    do BCV = 1, DIMENSION_BC
       
       select case (trim(BC_TYPE(BCV)))
       
          ! Determine the number of mass inlets that contain DEM solids.
       case ('MASS_INFLOW')
          M_LP: do M=1,M_TOT
             if(SOLIDS_MODEL(M)=='DEM' .and.                         &
                  BC_EP_s(BCV,M) > ZERO) then
                DEM_BCMI = DEM_BCMI + 1
                DEM_BCMI_MAP(DEM_BCMI) = BCV
                exit M_LP
             endif
          enddo M_LP
          
          ! Count the number of pressure outflows.
       case ('P_OUTFLOW','MASS_OUTFLOW')
          if(BC_PO_APPLY_TO_DES(BCV)) then
             DEM_BCMO = DEM_BCMO + 1
             DEM_BCMO_MAP(DEM_BCMO) = BCV
          endif

          ! Flag CG_MI as an error if DEM solids are present.
       case ('CG_MI')
          do M=1,M_TOT
             if(SOLIDS_MODEL(M)=='DEM') then
                if(IS_DEFINED(BC_EP_s(BCV,M)) .and.                 &
                     BC_EP_s(BCV,M) > ZERO) then
                   write(ERR_MSG,1100) trim(iVar('BC_TYPE',BCV)),    &
                        'GC_MI', trim( IFILE_NAME )
                   call flush_err_msg(ABORT=.true.)
                endif
             endif
          enddo

       case ('CG_PO')
          write(ERR_MSG,1100) trim(iVar('BC_TYPE',BCV)), 'GC_PO', &
               & trim( IFILE_NAME )
          call flush_err_msg(ABORT=.true.)

       case ('OUTFLOW', 'P_INFLOW')
          write(ERR_MSG,1100) trim(iVar('BC_TYPE',BCV)),             &
               trim(BC_TYPE(BCV)), trim( IFILE_NAME )
          call flush_err_msg(ABORT=.true.)

       end select

    enddo

    call finl_err_msg

1100 format('Error 1100: Unsupported boundary condition specified ',  &
          'with',/'DEM simulation: ',A,' = ',A,/'Please correct the ',&
          A,' file.')
    
  end subroutine check_bc_dem

end module check_bc_dem_module
