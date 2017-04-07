module check_boundary_conditions_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   use check_bc_dem_module, only: check_bc_dem
   use check_bc_inflow_module, only: check_bc_mass_inflow, check_bc_p_inflow
   use check_bc_outflow_module, only: check_bc_outflow, check_bc_mass_outflow, check_bc_p_outflow
   use check_bc_walls_module, only: check_bc_walls

! Parameter constants
   use param1, only: ZERO, ONE, UNDEFINED, IS_DEFINED, IS_UNDEFINED, EQUAL

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BOUNDARY_CONDITIONS                               !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Check boundary condition specifications                    !
!     - convert physical locations to i, j, k's                        !
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
!     - check specification of physical quantities                     !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine check_boundary_conditions(dx,dy,dz,xlength,ylength,zlength,domlo,domhi)

! Global Variables:
!---------------------------------------------------------------------//
! Total number of (actual) continuum solids.
      use constant, only: MMAX
! Flag: BC dimensions or Type is specified
      use bc, only: BC_DEFINED
! Use specified BC type
      use bc, only: BC_TYPE
! User specified BC solids bulk density
      use bc, only: BC_ROP_s
! Solids volume fraction at BC
      use bc, only: BC_EP_s
      use bc, only: BC_EP_g
! Run-time flag for DEM solids
      use run, only: DEM_SOLIDS

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of BCs
      use param, only: DIMENSION_BC
! Maximum number of disperse phases
      use param, only: DIM_M

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival

      use check_bc_geometry_module, only: check_bc_geometry, check_bc_geometry_flow, check_bc_geometry_wall

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
      do bcv = 1, dimension_bc

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

            case ('OUTFLOW', 'OF')
               call check_bc_geometry_flow(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)
               call check_bc_outflow(mmax, bcv)

            case ('MASS_OUTFLOW', 'MO')
               call check_bc_geometry_flow(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)
               call check_bc_mass_outflow(mmax, bcv)
               call check_bc_outflow(mmax, bcv)

            case ('P_OUTFLOW','PO')
               call check_bc_geometry_flow(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)
               call check_bc_p_outflow(bcv)
               call check_bc_outflow(mmax, bcv)

            case ('FREE_SLIP_WALL','FSW','PAR_SLIP_WALL','PSW','NO_SLIP_WALL','NSW')
               call check_bc_geometry_wall(bcv,dx,dy,dz,&
                  xlength,ylength,zlength,domlo,domhi)
               call check_bc_walls(bcv)

            end select

! Check whether BC values are specified for undefined BC locations
         elseif(bc_type(bcv) /= 'DUMMY') then

            call check_bc_range(bcv)

         endif
      enddo
! Additional checks needed for DEM boundaries
      if(dem_solids) call check_bc_dem(mmax)

! Cleanup and exit.
      call finl_err_msg

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_RANGE                                          !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Verify that data was not given for undefined BC regions.   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine check_bc_range(BCV)

      ! Gas phase BC varaibles
      use bc, only: BC_EP_g, BC_T_g, BC_X_g, BC_P_g
      use bc, only: BC_U_g, BC_V_g, BC_W_g

      ! Solids phase BC variables.
      use bc, only: BC_EP_s, BC_ROP_s, BC_T_s, BC_X_s
      use bc, only: BC_U_s, BC_V_s, BC_W_s


! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of disperse phases.
      use param, only: DIM_M
! Maximum number of species gas/solids
      use param, only: DIMENSION_N_G, DIMENSION_N_S


! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar, ival, err_msg


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
      IF(IS_DEFINED(BC_U_G(BCV))) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_U_g',BCV))
         call FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      IF(IS_DEFINED(BC_V_G(BCV))) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_V_g',BCV))
         call FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      IF (IS_DEFINED(BC_W_G(BCV))) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_W_g',BCV))
         call FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      IF (IS_DEFINED(BC_EP_G(BCV))) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_EP_g',BCV))
         call FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      IF (IS_DEFINED(BC_P_G(BCV))) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_P_g',BCV))
         call FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      IF (IS_DEFINED(BC_T_G(BCV))) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_T_g',BCV))
         call FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      DO N = 1, DIMENSION_N_G
         IF(IS_DEFINED(BC_X_G(BCV,N))) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_X_g',BCV,N))
            call FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

! Check solids phase variables.
      DO M = 1, DIM_M
         IF(IS_DEFINED(BC_ROP_S(BCV,M))) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_ROP_s',BCV,M))
            call FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(IS_DEFINED(BC_EP_S(BCV,M))) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_EP_s',BCV,M))
            call FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(IS_DEFINED(BC_U_S(BCV,M))) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_U_s',BCV,M))
            call FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(IS_DEFINED(BC_V_S(BCV,M))) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_V_s',BCV,M))
            call FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(IS_DEFINED(BC_W_S(BCV,M))) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_W_s',BCV,M))
            call FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(IS_DEFINED(BC_T_S(BCV,M))) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_T_s',BCV,M))
            call FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         DO N = 1, DIMENSION_N_S
            IF(IS_DEFINED(BC_X_S(BCV,M,N))) THEN
               WRITE(ERR_MSG,1100) trim(iVar('BC_X_s',BCV,M,N))
               call FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      ENDDO

      call finl_err_msg

 1100 FORMAT('Error 1100:',A,' specified for an undefined BC location')

      end subroutine check_bc_range

   end subroutine check_boundary_conditions
end module check_boundary_conditions_module
