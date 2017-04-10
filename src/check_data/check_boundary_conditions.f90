module check_boundary_conditions_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

! Parameter constants
   use param1, only: zero, one, undefined, is_defined, is_undefined, equal

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
   use error_manager, only: init_err_msg, finl_err_msg, flush_err_msg
   use error_manager, only: err_msg, ivar, ival

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BOUNDARY_CONDITIONS                               !
!                                                                      !
!  Purpose: Check boundary condition specifications                    !
!     - convert physical locations to i, j, k's                        !
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine check_boundary_conditions(dx, dy, dz, &
         xlength, ylength, zlength, domlo, domhi)

! Global Variables:
!---------------------------------------------------------------------//
! Total number of (actual) continuum solids.
      use constant, only: MMAX
! Flag: BC dimensions or Type is specified
      use bc, only: BC_DEFINED
! Use specified BC type
      use bc, only: BC_TYPE
! Solids volume fraction at BC
      use bc, only: BC_EP_s

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of BCs
      use param, only: DIM_BC
! Maximum number of disperse phases
      use param, only: DIM_M

      use check_bc_geometry_module, only: check_bc_geometry
      use check_bc_geometry_module, only: check_bc_geometry_flow
      use check_bc_geometry_module, only: check_bc_geometry_wall

      use check_bc_walls_module,    only: check_bc_walls

      use check_bc_inflow_module,   only: check_bc_mass_inflow
      use check_bc_inflow_module,   only: check_bc_p_inflow

      use check_bc_outflow_module,  only: check_bc_outflow
      use check_bc_outflow_module,  only: check_bc_p_outflow

      implicit none

      integer(c_int), intent(in) :: domlo(3),domhi(3)
      real(c_real)  , intent(in) :: dx,dy,dz
      real(c_real)  , intent(in) :: xlength,ylength,zlength

! Local Variables:
!---------------------------------------------------------------------//
! Loop counters
      integer :: BCV, I
! Flag to skip checks on indexed solid phase.
      logical :: SKIP(1:DIM_M)
!......................................................................!

! Initialize the error manager.
      call init_err_msg("CHECK_BOUNDARY_CONDITIONS")

! Determine which BCs are DEFINED
      call check_bc_geometry

! Loop over each defined BC and check the user data.
      do bcv = 1, dim_bc

         if (bc_defined(bcv)) then

! Determine which solids phases are present.
            do i = 1, dim_m
               skip(i) = equal(bc_ep_s(bcv,i), zero)
            enddo

            select case (trim(bc_type(bcv)))

            case ('MASS_INFLOW', 'MI')
               call check_bc_geometry_flow(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)
               call check_bc_mass_inflow(mmax, skip, bcv)
               ! call check_bc_inflow(mmax,skip,bcv)

            case ('P_INFLOW', 'PI')
               call check_bc_geometry_flow(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)
               call check_bc_p_inflow(mmax, skip, bcv)
               ! call check_bc_inflow(mmax, skip, bcv)
               call check_bc_outflow(mmax, bcv)

            case ('P_OUTFLOW','PO')
               call check_bc_geometry_flow(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)
               call check_bc_p_outflow(bcv)
               call check_bc_outflow(mmax, bcv)

            case ('FREE_SLIP_WALL','FSW',&
                  'PAR_SLIP_WALL', 'PSW',&
                  'NO_SLIP_WALL',  'NSW')
               call check_bc_geometry_wall(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)
               call check_bc_walls(bcv)

            end select

! Check whether BC values are specified for undefined BC locations
         elseif(bc_type(bcv) /= 'DUMMY') then

            call check_bc_range(bcv)

         endif
      enddo

! Cleanup and exit.
      call finl_err_msg

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_RANGE                                          !
!                                                                      !
!  Purpose: Verify that data was not given for undefined BC regions.   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine check_bc_range(BCV)

      ! Gas phase BC varaibles
      use bc, only: BC_EP_g, BC_T_g, BC_X_g, BC_P_g
      use bc, only: BC_U_g, BC_V_g, BC_W_g

      ! Solids phase BC variables.
      use bc, only: BC_EP_s, BC_T_s, BC_X_s
      use bc, only: BC_U_s, BC_V_s, BC_W_s

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of disperse phases.
      use param, only: DIM_M
! Maximum number of species gas/solids
      use param, only: DIMENSION_N_G, DIMENSION_N_S

      use param1, only: zero, one, equal

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------/
! Boundary condition index.
      integer, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Generic loop varaibles.
      integer :: M, N
!......................................................................!


! Initialize the error manager.
      call INIT_ERR_MSG("CHECK_BC_RANGE")


! Check gas phase variables.
      if(.not.equal(bc_u_g(bcv),zero)) then
         write(err_msg,1100) trim(ivar('BC_U_g',bcv))
         call flush_err_msg(abort=.true.)
      endif
      if(.not.equal(bc_v_g(bcv),zero)) then
         write(err_msg,1100) trim(ivar('BC_V_g',bcv))
         call flush_err_msg(abort=.true.)
      endif
      if (.not.equal(bc_w_g(bcv),zero)) then
         write(err_msg,1100) trim(ivar('BC_W_g',bcv))
         call flush_err_msg(abort=.true.)
      endif
      if (.not.equal(bc_ep_g(bcv),one)) then
         write(err_msg,1100) trim(ivar('BC_EP_g',bcv))
         call flush_err_msg(abort=.true.)
      endif
      if (is_defined(bc_p_g(bcv))) then
         write(err_msg,1100) trim(ivar('BC_P_g',bcv))
         call flush_err_msg(abort=.true.)
      endif
      if (is_defined(bc_t_g(bcv))) then
         write(err_msg,1100) trim(ivar('BC_T_g',bcv))
         call flush_err_msg(abort=.true.)
      ENDIF

      DO N = 1, DIMENSION_N_G
         IF(IS_DEFINED(BC_X_G(BCV,N))) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_X_g',BCV,N))
            call FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

! Check solids phase variables.
      do m = 1, dim_m
         if(.not.equal(bc_ep_s(bcv,m),zero)) then
            write(err_msg,1100) trim(ivar('BC_EP_s',BCV,M))
            call flush_err_msg(abort=.true.)
         endif
         if(.not.equal(bc_u_s(bcv,m),zero)) then
            write(err_msg,1100) trim(ivar('BC_U_s',BCV,M))
            call flush_err_msg(abort=.true.)
         endif
         if(.not.equal(bc_v_s(bcv,m),zero)) then
            write(err_msg,1100) trim(ivar('BC_V_s',BCV,M))
            call flush_err_msg(abort=.true.)
         endif

         if(.not.equal(bc_w_s(bcv,m),zero)) then
            write(err_msg,1100) trim(ivar('BC_W_s',BCV,M))
            call flush_err_msg(abort=.true.)
         endif
         if(is_defined(bc_t_s(bcv,m))) then
            write(err_msg,1100) trim(ivar('BC_T_s',BCV,M))
            call flush_err_msg(abort=.true.)
         endif

         do n = 1, dimension_n_s
            if(is_defined(bc_x_s(bcv,m,n))) then
               write(err_msg,1100) trim(ivar('BC_X_s',bcv,m,n))
               call flush_err_msg(abort=.true.)
            endif
         enddo

      ENDDO

      call finl_err_msg

 1100 FORMAT('Error 1100:',A,' specified for an undefined BC location')

      end subroutine check_bc_range

   end subroutine check_boundary_conditions
end module check_boundary_conditions_module
