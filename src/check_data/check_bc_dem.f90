module check_bc_dem_module

  use amrex_fort_module, only : c_real => amrex_real
  use iso_c_binding , only: c_int
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
       case ('MASS_INFLOW', 'MI')
          m_lp: do m=1,m_tot
             if(SOLIDS_MODEL(M)=='DEM' .and.                         &
                  bc_ep_s(bcv,m) > zero) then
                dem_bcmi = dem_bcmi + 1
                dem_bcmi_map(dem_bcmi) = bcv
                exit m_lp
             endif
          enddo m_lp

       ! Count the number of pressure outflows.
       case ('P_OUTFLOW','PO','MASS_OUTFLOW','MO')
          if(bc_po_apply_to_des(bcv)) then
             dem_bcmo = dem_bcmo + 1
             dem_bcmo_map(dem_bcmo) = bcv
          endif

       case ('P_INFLOW','PI')
          write(err_msg,1100) trim(ivar('BC_TYPE',BCV)),             &
               trim(bc_type(bcv))
          call flush_err_msg(abort=.true.)

       end select

    enddo

    call finl_err_msg

1100 format('Error 1100: Unsupported boundary condition specified ',  &
        'with',/'DEM simulation: ',A,' = ',A,/&
        'Please correct the input deck.')

  end subroutine check_bc_dem

end module check_bc_dem_module
