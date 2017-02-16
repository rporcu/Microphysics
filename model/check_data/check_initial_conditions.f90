MODULE CHECK_INITIAL_CONDITIONS_MODULE

  use bl_fort_module,                  only: c_real
  use iso_c_binding ,                  only: c_int
  use check_ic_common_discrete_module, only: check_ic_common_discrete
  use run,                             only: IFILE_NAME
  use param1,                          only: UNDEFINED, UNDEFINED_I,      &
                                           & IS_DEFINED, IS_UNDEFINED,    & 
                                           & ZERO, ONE
  use error_manager,                   only: finl_err_msg, flush_err_msg, & 
                                           & init_err_msg, ivar, ival,    &
                                           & err_msg

  implicit none
  private 

  public check_initial_conditions

contains
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                      !
  !  Subroutine: CHECK_INITIAL_CONDITIONS                                !
  !  Author: P. Nicoletti                               Date: 02-DEC-91  !
  !  Author: J.Musser                                   Date: 01-MAR-14  !
  !                                                                      !
  !  Purpose: check the initial conditions input section                 !
  !     - check geometry of any specified IC region                      !
  !     - check specification of physical quantities                     !
  !                                                                      !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_initial_conditions(dx,dy,dz)

    use ic,                    only: IC_DEFINED
    use run,                   only: DEM_SOLIDS, RUN_TYPE
    use param,                 only: DIMENSION_IC
    use calc_cell_module,      only: calc_cell
    use location_check_module, only: location_check


    real(c_real), intent(in) :: dx,dy,dz
    integer                  :: ICV


    ! Initialize the error manager.
    call init_err_msg("CHECK_INITIAL_CONDITIONS")

    ! Determine which ICs are DEFINED
    call check_ic_geometry(dx,dy,dz)

    ! Loop over all IC arrays.
    do ICV=1, DIMENSION_IC

       ! Verify user input for defined defined IC.
       if(IC_DEFINED(ICV)) then
          ! Gas phase checks.
          call check_ic_gas_phase(ICV)
          ! Generic solids phase checks.
          call check_ic_solids_phases(ICV)
          
          ! Verify that no data was defined for unspecified IC. ICs are only
          ! defined for new runs, so these checks are restricted to new runs.
       elseif(RUN_TYPE == 'NEW') then
          call check_ic_overflow(ICV)
       endif
    enddo
    
    ! Check the initial conditions for the DEM model as well
    if(DEM_SOLIDS) call check_ic_common_discrete
    
    ! Finalize the error manager.
    call finl_err_msg
        
  contains

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    !                                                                      !
    ! Subroutine: CHECK_IC_GEOMETRY                                        !
    ! Author: J.Musser                                    Date: 01-Mar-14  !
    !                                                                      !
    ! Purpose: Provided a detailed error message when the sum of volume    !
    !                                                                      !
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
    subroutine check_ic_geometry(dx,dy,dz)

      use geometry, only: domlo, domhi
      use param,    only: DIMENSION_IC
      use run,      only: RUN_TYPE
      use ic,       only: IC_DEFINED, IC_TYPE,             &
                        & IC_X_e, IC_X_w, IC_I_e, IC_I_w,  &
                        & IC_Y_n, IC_Y_s, IC_J_n, IC_J_s,  &
                        & IC_Z_t, IC_Z_b, IC_K_t, IC_K_b

      real(c_real), intent(in) :: dx,dy,dz
      integer                  :: ICV, I_w, I_e, J_s, J_n, K_b, K_t


      ! Initialize the error manager.
      call init_err_msg("CHECK_IC_GEOMETRY")

      ! Check geometry of any specified IC region
      do ICV = 1, DIMENSION_IC

         IC_DEFINED(ICV) = .false.
         if (IS_DEFINED(IC_X_W(ICV))) IC_DEFINED(ICV) = .true.
         if (IS_DEFINED(IC_X_E(ICV))) IC_DEFINED(ICV) = .true.
         if (IS_DEFINED(IC_Y_S(ICV))) IC_DEFINED(ICV) = .true.
         if (IS_DEFINED(IC_Y_N(ICV))) IC_DEFINED(ICV) = .true.
         if (IS_DEFINED(IC_Z_B(ICV))) IC_DEFINED(ICV) = .true.
         if (IS_DEFINED(IC_Z_T(ICV))) IC_DEFINED(ICV) = .true.
         if (IS_DEFINED(IC_I_W(ICV))) IC_DEFINED(ICV) = .true.
         if (IS_DEFINED(IC_I_E(ICV))) IC_DEFINED(ICV) = .true.
         if (IS_DEFINED(IC_J_S(ICV))) IC_DEFINED(ICV) = .true.
         if (IS_DEFINED(IC_J_N(ICV))) IC_DEFINED(ICV) = .true.
         if (IS_DEFINED(IC_K_B(ICV))) IC_DEFINED(ICV) = .true.
         if (IS_DEFINED(IC_K_T(ICV))) IC_DEFINED(ICV) = .true.

         ! An IC is defined for restart runs only if it is a 'PATCH'.
         if(RUN_TYPE /= 'NEW' .and. IC_TYPE(ICV) /= 'PATCH') &
              IC_DEFINED(ICV) = .false.
         
         ! Ignore patched IC regions for new runs. It may be better to flag this as
         ! and error to avoid user confusion.
         if(RUN_TYPE == 'NEW' .and. IC_TYPE(ICV) == 'PATCH') &
              IC_DEFINED(ICV) = .FALSE.

         IF(.NOT.IC_DEFINED(ICV)) CYCLE

         if (IS_UNDEFINED(IC_X_W(ICV)) .and. IS_UNDEFINED(IC_I_W(ICV))) then
            write(ERR_MSG, 1100) ICV, 'IC_X_w and IC_I_w', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif

         if (IS_UNDEFINED(IC_X_E(ICV)) .and. IS_UNDEFINED(IC_I_E(ICV))) then
            write(ERR_MSG, 1100) ICV, 'IC_X_e and IC_I_e', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif

         if (IS_UNDEFINED(IC_Y_S(ICV)) .and. IS_UNDEFINED(IC_J_S(ICV))) then
            write(ERR_MSG, 1100) ICV, 'IC_Y_s and IC_J_s', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif

         if (IS_UNDEFINED(IC_Y_N(ICV)) .and. IS_UNDEFINED(IC_J_N(ICV))) then
            write(ERR_MSG, 1100) ICV, 'IC_Y_n and IC_J_n', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif

         if (IS_UNDEFINED(IC_Z_B(ICV)) .and. IS_UNDEFINED(IC_K_B(ICV))) then
            write(ERR_MSG, 1100) ICV, 'IC_Z_b and IC_K_b', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif

         if (IS_UNDEFINED(IC_Z_T(ICV)) .and. IS_UNDEFINED(IC_K_T(ICV))) then
            write(ERR_MSG, 1100) ICV, 'IC_Z_t and IC_K_t', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif
         
      enddo   ! end loop over (icv = 1,dimension_ic)

1100  format('Error 1100: Initial condition region ',I3,' is ill-',    &
           'defined.',/' > ',A,' are not specified.',/'Please correct ', &
           'the ',A,' file.')


      do ICV = 1, DIMENSION_IC

         ! Skip this check if the IC region is not specified.
         if(.not.IC_DEFINED(ICV)) cycle
         
         if (IS_DEFINED(IC_X_W(ICV)).and. IS_DEFINED(IC_X_E(ICV))) then
            I_W = CALC_CELL (IC_X_W(ICV), DX)
            I_W = I_W + 1
            I_E = CALC_CELL (IC_X_E(ICV), DX)
            if (IC_I_W(ICV)/=UNDEFINED_I .or. IC_I_E(ICV)/=UNDEFINED_I) then
               call location_check (IC_I_W(ICV), I_W, ICV, 'IC - west')
               call location_check (IC_I_E(ICV), I_E, ICV, 'IC - east')
            else
               IC_I_W(ICV) = I_W
               IC_I_E(ICV) = I_E
            endif
         endif
         
         ! Report problems with calculated bounds.
         if(IC_I_W(ICV) > IC_I_E(ICV)) then
            write(ERR_MSG, 1101) ICV, 'IC_I_W > IC_I_E', trim(IFILE_NAME)
            write(*,*)' dump:',IC_I_W(ICV),IC_I_E(ICV)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_I_W(ICV) < domlo(1)) then
            write(ERR_MSG, 1101) ICV, 'IC_I_W < domlo(1)', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_I_W(ICV) > domhi(1)) then
            write(ERR_MSG, 1101) ICV, 'IC_I_W > domhi(1)', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_I_E(ICV) < domlo(1)) then
            write(ERR_MSG, 1101) ICV, 'IC_I_E < domlo(1)', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_I_E(ICV) > domhi(1)) then
            write(ERR_MSG, 1101) ICV, 'IC_Z_t and IC_K_t', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif
         
         if (IS_DEFINED(IC_Y_S(ICV)) .and. IS_DEFINED(IC_Y_N(ICV))) then
            J_S = CALC_CELL (IC_Y_S(ICV), DY)
            J_S = J_S + 1
            J_N = CALC_CELL (IC_Y_N(ICV), DY)
            if (IC_J_S(ICV)/=UNDEFINED_I .or. IC_J_N(ICV)/=UNDEFINED_I) then
               call location_check (IC_J_S(ICV), J_S, ICV, 'IC - south')
               call location_check (IC_J_N(ICV), J_N, ICV, 'IC - north')
            else
               IC_J_S(ICV) = J_S
               IC_J_N(ICV) = J_N
            endif
         endif
         
         if(IC_J_S(ICV) > IC_J_N(ICV)) then
            write(ERR_MSG, 1101) ICV, 'IC_J_S > IC_J_N', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_J_S(ICV)<domlo(2)) then
            write(ERR_MSG, 1101) ICV, 'IC_J_S < domlo(2)', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_J_S(ICV)>domhi(2)) then
            write(ERR_MSG, 1101) ICV, 'IC_J_S >  domhi(2)', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_J_N(ICV)<domlo(2)) then
            write(ERR_MSG, 1101) ICV, 'IC_J_N < domlo(2)', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_J_N(ICV)>domhi(2)) then
            write(ERR_MSG, 1101) ICV, 'IC_J_N > domhi(2)', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif
         
         
         if (IS_DEFINED(IC_Z_B(ICV)) .and. IS_DEFINED(IC_Z_T(ICV))) then
            K_B = CALC_CELL (IC_Z_B(ICV), DZ)
            K_B = K_B + 1
            K_T = CALC_CELL (IC_Z_T(ICV), DZ)
            if (IC_K_B(ICV)/=UNDEFINED_I .or. IC_K_T(ICV)/=UNDEFINED_I) then
               call location_check (IC_K_B(ICV), K_B, ICV, 'IC - bottom')
               call location_check (IC_K_T(ICV), K_T, ICV, 'IC - top')
            else
               IC_K_B(ICV) = K_B
               IC_K_T(ICV) = K_T
            endif
         endif
         
         if(IC_K_B(ICV) > IC_K_T(ICV)) then
            write(ERR_MSG, 1101) ICV, 'IC_K_B > IC_K_T', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_K_B(ICV) < domlo(3)) then
            write(ERR_MSG, 1101) ICV, 'IC_K_B < domlo(3)', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_K_B(ICV) > domhi(3)) then
            write(ERR_MSG, 1101) ICV, 'IC_K_B > domhi(3)', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_K_T(ICV) < domlo(3)) then
            write(ERR_MSG, 1101) ICV, 'IC_K_T < domlo(3)', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         elseif(IC_K_T(ICV) > domhi(3)) then
            write(ERR_MSG, 1101) ICV, 'IC_K_T > domhi(3)', trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif
         
         
1101     format('Error 1101: Initial condition region ',I2,' is ill-',    &
              'defined.',/3x,A,/'Please correct the ',A,' file.')
         
      enddo   ! end loop over (icv=1,dimension_ic)
      
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
      use ic,        only: IC_EP_g, IC_P_g, IC_U_g, IC_V_g, IC_W_g,  IC_TYPE

      integer, intent(in) :: ICV
      logical             :: BASIC_IC
      
      call init_err_msg("CHECK_IC_GAS_PHASE")
      
      ! Patch ICs skip various checks.
      BASIC_IC = (IC_TYPE(ICV) /= 'PATCH')
      
      ! Check that gas phase velocity components are initialized.
      if(BASIC_IC) then
         if(IS_UNDEFINED(IC_U_G(ICV))) then
            write(ERR_MSG, 1000) trim(iVar('IC_U_g',ICV)), trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif
         
         if(IS_UNDEFINED(IC_V_G(ICV))) then
            write(ERR_MSG, 1000) trim(iVar('IC_V_g',ICV)), trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif
         
         if(IS_UNDEFINED(IC_W_G(ICV))) then
            write(ERR_MSG, 1000) trim(iVar('IC_W_g',ICV)), trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif
      endif
      
      ! Check that gas phase void fraction is initialized. Patched ICs may
      ! have an undefined volume fration. A second check is preformed on
      ! the solids.
      if(IS_UNDEFINED(IC_EP_G(ICV)) .and. BASIC_IC) then
         write(ERR_MSG, 1000) trim(iVar('IC_EP_g',ICV)), trim(IFILE_NAME)
         call flush_err_msg(ABORT=.true.)
      endif
      
      ! Check that if the gas phase pressure is initialized and the gas is
      ! compressible that the gas phase pressure is not zero or negative
      if(IS_DEFINED(IC_P_G(ICV))) then
         if(IS_UNDEFINED(RO_G0).and. IC_P_G(ICV)<=ZERO) then
            write(ERR_MSG, 1100) trim(iVar('IC_P_g',ICV)),             &
                 iVal(IC_P_G(ICV)), trim(IFILE_NAME)
            call flush_err_msg(ABORT=.true.)
         endif
      endif
      
1100  format('Error 1100: Pressure must be greater than 0.0 for ',     &
           'compressible flow',/'Illegal value: ',A,' = ',A,/'Please ',  &
           'correct the ',A,' file.')

      call finl_err_msg


1000  format('Error 1000: Required input not specified: ',A,/'Please ',&
           'correct the ',A,' file.')

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


      use param,    only: DIM_M
      use toleranc, only: compare
      use constant, only: MMAX,ro_s0
      use ic,       only: IC_EP_s, IC_ROP_s, IC_U_s, IC_V_s, IC_W_s, &
                        & IC_EP_g, IC_TYPE

      integer, intent(IN) :: ICV
      integer             :: M
      real(c_real)        :: SUM_EP, IC_ROs(1:DIM_M)
      logical             :: SKIP(1:DIM_M), BASIC_IC

      ! Initialize the error manager.
      CALL init_err_msg("CHECK_IC_SOLIDS_PHASES")

      ! Patch ICs skip various checks.
      BASIC_IC = (IC_TYPE(ICV) /= 'PATCH')

      ! Calculate EP_s from EP_g if there is only one solids phase.
      if(MMAX == 1 .and. IS_DEFINED(IC_EP_S(ICV,1))) then
         if(IS_DEFINED(IC_EP_g(ICV))) IC_EP_S(ICV,1) = ONE-IC_EP_g(ICV)
      endif
      
      ! Bulk density or solids volume fraction must be explicitly defined
      ! if there are more than one solids phase.
      if(MMAX > 1 .and. .not.COMPARE(IC_EP_g(ICV),ONE)) then
         ! IC_EP_g may be undefined for PATCH IC regions.
         if(IS_DEFINED(IC_EP_g(ICV))) then
            do M = 1, MMAX
               if(IS_DEFINED(IC_ROP_S(ICV,M)) .and. &
                    IS_DEFINED(IC_EP_S(ICV,M))) then
                  write(ERR_MSG, 1400) M, ICV, 'IC_ROP_s and IC_EP_s', trim(IFILE_NAME)
                  call flush_err_msg(ABORT=.true.)
               endif
            enddo
            
            ! If IC_EP_G is undefined, then ROP_s and EP_s should be too.
         else
            do M = 1, MMAX
               if(IS_DEFINED(IC_ROP_S(ICV,M)) .and. &
                    IS_DEFINED(IC_EP_S(ICV,M))) then
                  write(ERR_MSG, 1400) M, ICV, 'IC_ROP_s and IC_EP_s', trim(IFILE_NAME)
                  call flush_err_msg(ABORT=.true.)
               endif
            enddo
         endif
      endif
      
1400  format('Error 1400: Insufficient solids phase ',I2,' ',          &
           'information for IC',/'region ',I3,'. ',A,' not specified.',/ &
           'Please correct the ',A,' file.')
      
      ! Determine which solids phases are present.
      do M = 1, MMAX
         SKIP(M)=(IS_UNDEFINED(IC_ROP_S(ICV,M)).or.abs(IC_ROP_S(ICV,M))<epsilon(ZERO)) &
              .and.(IS_UNDEFINED(IC_EP_S(ICV,M)) .or.abs(IC_EP_S(ICV,M))<epsilon(ZERO))
      enddo
      
      if(MMAX == 1 .and. abs(IC_EP_g(ICV)-ONE)>ZERO) SKIP(1) = .false.
      
      do M=1, MMAX
         
         ! check that solids phase m velocity components are initialized
         if(BASIC_IC) then
            if(IS_UNDEFINED(IC_U_S(ICV,M))) then
               if (SKIP(M)) then
                  IC_U_S(ICV,M) = ZERO
               else
                  write(ERR_MSG, 1000)trim(iVar('IC_U_s',ICV,M)), trim(IFILE_NAME)
                  call flush_err_msg(ABORT=.true.)
               endif
            endif
            
            if(IS_UNDEFINED(IC_V_S(ICV,M))) then
               if(SKIP(M)) then
                  IC_V_S(ICV,M) = ZERO
               else
                  write(ERR_MSG, 1000)trim(iVar('IC_V_s',ICV,M)), trim(IFILE_NAME)
                  call flush_err_msg(ABORT=.true.)
               endif
            endif
            
            if(IS_UNDEFINED(IC_W_S(ICV,M))) then
               if(SKIP(M)) then
                  IC_W_S(ICV,M) = ZERO
               else
                  write(ERR_MSG, 1000)trim(iVar('IC_W_s',ICV,M)), trim(IFILE_NAME)
                  call flush_err_msg(ABORT=.true.)
               endif
            endif
            
         endif
         

         IC_ROs(M) = RO_s0(M)
         
      enddo   ! end loop over (m=1)
      
      
      ! Initialize the sum of the total volume fraction.
      SUM_EP = IC_EP_G(ICV)
      
      do M=1, MMAX
         
         ! Clear out both varaibles if this phase is skipped.
         if(BASIC_IC .and. SKIP(M)) then
            IC_EP_S(ICV,M)  = ZERO
            IC_ROP_S(ICV,M) = ZERO
            
            ! Leave everything undefined for PATCH ICs that are not specifed.
         elseif(.not.BASIC_IC .and. (IS_UNDEFINED(IC_ROP_S(ICV,M))      &
              .and. IS_UNDEFINED(IC_EP_S(ICV,M)))) then
            
            ! If both input parameters are defined. Make sure they are equivalent.
         elseif(IS_DEFINED(IC_ROP_S(ICV,M)) .and.                     &
              IS_DEFINED(IC_EP_S(ICV,M))) then
            
            if(.not.COMPARE(IC_EP_S(ICV,M)*IC_ROs(M),                  &
                 IC_ROP_S(ICV,M))) then
               
               ! BASIC_IC regions require that the IC_ROP_s and IC_EP_s specifications
               ! match although it is unlikely that anyone would specify both.
               if(BASIC_IC) then
                  
                  write(ERR_MSG,1406) M, ICV, trim(IFILE_NAME)
                  call flush_err_msg(ABORT=.true.)
                  
1406              format('Error 1406: IC_EP_s and IC_ROP_s are inconsistent for ',&
                       'phase ',I2,/,'in IC region ', I3,'. Please correct the ',&
                        A, ' file.')
                  
                  
                  ! PachedeIC regions defer to IC_EP_s if the values do not match. This
                  ! prevents a dead lock or the need to define both. This case is rather
                  ! common as a defined IC_EP_s is converted to IC_ROP_s. Therefore, if
                  ! a patch region is used more than once, these values may not match.
               else
                  
                  write(ERR_MSG,1407) trim(iVar('IC_ROP_s',ICV,M)), &
                       trim(iVAL(IC_ROP_S(ICV,M))), trim(iVar('IC_EP_s',&
                       ICV,M)), trim(iVAL(IC_EP_S(ICV,M)))
                  call flush_err_msg()
                  
1407              format('Warning 1407: IC_EP_s and IC_ROP_s are inconsistent:',    &
                       2(/3x,A,' = ',A),/'Deferring to IC_EP_s to overcome conflict.')
                  
                  IC_ROP_S(ICV,M) = IC_EP_S(ICV,M)*IC_ROs(M)
                  
               endif
            endif
            
            
            ! Compute IC_EP_s from IC_ROP_s
         elseif(IS_UNDEFINED(IC_EP_S(ICV,M)))then
            IC_EP_S(ICV,M) = IC_ROP_S(ICV,M) / IC_ROs(M)
            
            ! Compute IC_ROP_s from IC_EP_s and IC_ROs
         elseif(IS_UNDEFINED(IC_ROP_S(ICV,M))) then
            IC_ROP_S(ICV,M) = IC_EP_S(ICV,M) * IC_ROs(M)
            ! This is a sanity check.
         else
            
         endif
         ! Add this phase to the total volume fraction.
         SUM_EP = SUM_EP + IC_EP_S(ICV,M)
      enddo
      
      ! Verify that the volume fractions sum to one.
      if(BASIC_IC .and. .not.COMPARE(SUM_EP,ONE)) then
         write(ERR_MSG,1410) ICV, trim(IFILE_NAME)
         call flush_err_msg(ABORT=.true.)
      endif
      
1410  format('Error 1410: Illegal initial condition region : ',I3,/    &
           'Sum of volume fractions does NOT equal ONE. Please correct',/&
           'the ',A,' file.')
      
      
      call finl_err_msg

1000  format('Error 1000: Required input not specified: ',A,/'Please ',&
           'correct the ',A,' file.')
      
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

      use param, only: DIM_M
      use ic,    only: IC_EP_g, IC_U_g, IC_V_g, IC_W_g, IC_ROP_s,  IC_U_s, &
                     & IC_V_s, IC_W_s, IC_TYPE

      integer, intent(in) :: ICV
      integer :: M

      if (IC_TYPE(ICV) == 'PATCH') return

      ! Initialize the error manager.
      call init_err_msg("CHECK_IC_OVERFLOW")

      ! GAS PHASE quantities
      ! -------------------------------------------->>>
      if(IS_DEFINED(IC_U_G(ICV))) then
         write(ERR_MSG, 1010) trim(iVar('IC_U_g',ICV))
         call flush_err_msg(ABORT=.true.)
      elseif(IS_DEFINED(IC_V_G(ICV))) then
         write(ERR_MSG, 1010) trim(iVar('IC_V_g',ICV))
         call flush_err_msg(ABORT=.true.)
      elseif(IS_DEFINED(IC_W_G(ICV))) then
         write(ERR_MSG, 1010) trim(iVar('IC_W_g',ICV))
         call flush_err_msg(ABORT=.true.)
      elseif(IS_DEFINED(IC_EP_G(ICV))) then
         write(ERR_MSG, 1010) trim(iVar('IC_EP_g',ICV))
         call flush_err_msg(ABORT=.true.)
      endif
      ! --------------------------------------------<<<
      
      
      ! SOLIDS PHASE quantities
      ! -------------------------------------------->>>
      do M=1, DIM_M
         if(IS_DEFINED(IC_ROP_S(ICV,M))) then
            write(ERR_MSG, 1010) trim(iVar('IC_ROP_s',ICV,M))
            call flush_err_msg(ABORT=.true.)
         elseif(IS_DEFINED(IC_U_S(ICV,M))) then
            write(ERR_MSG, 1010) trim(iVar('IC_U_s',ICV,M))
            call flush_err_msg(ABORT=.true.)
         elseif(IS_DEFINED(IC_V_S(ICV,M))) then
            write(ERR_MSG, 1010) trim(iVar('IC_V_s',ICV,M))
            call flush_err_msg(ABORT=.true.)
         elseif(IS_DEFINED(IC_W_S(ICV,M))) then
            write(ERR_MSG, 1010) trim(iVar('IC_W_s',ICV,M))
            call flush_err_msg(ABORT=.true.)
         endif
      enddo
      ! --------------------------------------------<<<
      
      call finl_err_msg
      
1010  format('Error 1010: ',A,' specified in an undefined IC region')

    end subroutine check_ic_overflow

  end subroutine check_initial_conditions

end module check_initial_conditions_module
