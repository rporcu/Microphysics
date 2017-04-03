MODULE CHECK_BC_GEOMETRY_MODULE

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   use calc_cell_module, only: calc_loc, calc_cell
   use check_plane_module, only: check_plane

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED, UNDEFINED_I, UNDEFINED_C, IS_UNDEFINED, IS_DEFINED, ZERO, EQUAL

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_GEOMETRY                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Determine if BCs are "DEFINED" and that they contain the    !
! minimum amount of geometry data.                                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_GEOMETRY

! Global Variables:
!---------------------------------------------------------------------//
! Flag: BC contains geometric data and/or specified type
      use bc, only: BC_DEFINED
! User specified BC
      use bc, only: BC_TYPE
! User specified: BC geometry
      use bc, only: BC_X_e, BC_X_w, BC_I_e, BC_I_w
      use bc, only: BC_Y_n, BC_Y_s, BC_J_n, BC_J_s
      use bc, only: BC_Z_t, BC_Z_b, BC_K_t, BC_K_b

! Global Parameters:
!---------------------------------------------------------------------//
! The max number of BCs.
      use param, only: DIMENSION_BC

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      use run, only: IFILE_NAME

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! loop/variable indices
      integer :: BCV, I
! Error flag
      logical :: RECOGNIZED_BC_TYPE
! Total number of valid BC types
      integer, PARAMETER :: DIM_BCTYPE = 21
! Valid boundary condition types
      CHARACTER(LEN=16), DIMENSION(1:DIM_BCTYPE) ::VALID_BC_TYPE = (/&
           'MASS_INFLOW     ', 'MI              ',&
           'MASS_OUTFLOW    ', 'MO              ',&
           'P_INFLOW        ', 'PI              ',&
           'P_OUTFLOW       ', 'PO              ',&
           'FREE_SLIP_WALL  ', 'FSW             ',&
           'NO_SLIP_WALL    ', 'NSW             ',&
           'PAR_SLIP_WALL   ', 'PSW             ',&
           'OUTFLOW         ', 'OF              ',&
           'CG_NSW          ', 'CG_FSW          ',&
           'CG_PSW          ', 'CG_MI           ',&
           'CG_PO           '/)
!......................................................................!


      CALL INIT_ERR_MSG("CHECK_BC_GEOMETRY")

      L50: DO BCV = 1, DIMENSION_BC

         BC_DEFINED(BCV) = .FALSE.
         IF(IS_DEFINED(BC_X_W(BCV))) BC_DEFINED(BCV) = .TRUE.
         IF(IS_DEFINED(BC_X_E(BCV))) BC_DEFINED(BCV) = .TRUE.
         IF(IS_DEFINED(BC_Y_S(BCV))) BC_DEFINED(BCV) = .TRUE.
         IF(IS_DEFINED(BC_Y_N(BCV))) BC_DEFINED(BCV) = .TRUE.
         IF(IS_DEFINED(BC_Z_B(BCV))) BC_DEFINED(BCV) = .TRUE.
         IF(IS_DEFINED(BC_Z_T(BCV))) BC_DEFINED(BCV) = .TRUE.
         IF(IS_DEFINED(BC_I_W(BCV))) BC_DEFINED(BCV) = .TRUE.
         IF(IS_DEFINED(BC_I_E(BCV))) BC_DEFINED(BCV) = .TRUE.
         IF(IS_DEFINED(BC_J_S(BCV))) BC_DEFINED(BCV) = .TRUE.
         IF(IS_DEFINED(BC_J_N(BCV))) BC_DEFINED(BCV) = .TRUE.
         IF(IS_DEFINED(BC_K_B(BCV))) BC_DEFINED(BCV) = .TRUE.
         IF(IS_DEFINED(BC_K_T(BCV))) BC_DEFINED(BCV) = .TRUE.

         IF (BC_TYPE(BCV) == 'DUMMY') BC_DEFINED(BCV) = .FALSE.

         IF(BC_TYPE(BCV)/=UNDEFINED_C .AND. BC_TYPE(BCV)/='DUMMY')THEN

            RECOGNIZED_BC_TYPE = .FALSE.
            DO I = 1, DIM_BCTYPE
                VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
                IF(VALID_BC_TYPE(I) == BC_TYPE(BCV)) THEN
                   RECOGNIZED_BC_TYPE = .TRUE.
                   EXIT
                ENDIF
            ENDDO

            IF(.NOT.RECOGNIZED_BC_TYPE) THEN
               WRITE(ERR_MSG, 1100) trim(iVar('BC_TYPE',BCV)), &
                  trim(BC_TYPE(BCV)), VALID_BC_TYPE
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(.NOT.BC_DEFINED(BCV)) CYCLE

         IF(IS_UNDEFINED(BC_X_W(BCV)) .AND. IS_UNDEFINED(BC_I_W(BCV))) THEN
            write(ERR_MSG,1101) BCV, 'BC_X_w and BC_I_w', trim( IFILE_NAME )
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(IS_UNDEFINED(BC_X_E(BCV)) .AND. IS_UNDEFINED(BC_I_E(BCV))) THEN
            write(ERR_MSG, 1101) BCV, 'BC_X_e and BC_I_e', trim( IFILE_NAME )
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(IS_UNDEFINED(BC_Y_S(BCV)) .AND. IS_UNDEFINED(BC_J_S(BCV))) THEN
            write(ERR_MSG, 1101) BCV, 'BC_Y_s and BC_J_s', trim( IFILE_NAME )
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(IS_UNDEFINED(BC_Y_N(BCV)) .AND. IS_UNDEFINED(BC_J_N(BCV))) THEN
            write(ERR_MSG, 1101) BCV, 'BC_Y_n and BC_J_n', trim( IFILE_NAME )
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(IS_UNDEFINED(BC_Z_B(BCV)) .AND. IS_UNDEFINED(BC_K_B(BCV))) THEN
            write(ERR_MSG, 1101) BCV, 'BC_Z_b and BC_K_b', trim( IFILE_NAME )
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(IS_UNDEFINED(BC_Z_T(BCV)) .AND. IS_UNDEFINED(BC_K_T(BCV))) THEN
            write(ERR_MSG, 1101) BCV, 'BC_Z_t and BC_K_t', trim( IFILE_NAME )
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1101 FORMAT('Error 1101: Boundary condition ',I3,' is ill-defined.',/ &
         A,' are not specified.',/'Please correct the ', A, ' file.')

! Swap BC aliases for the "full name" complement.
         DO I = 1, DIM_BCTYPE
            VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
            IF(VALID_BC_TYPE(I) == BC_TYPE(BCV)) THEN
               IF(MOD(I,2) == 0) BC_TYPE(BCV) = VALID_BC_TYPE(I-1)
               CYCLE  L50
            ENDIF
         ENDDO

         WRITE(ERR_MSG, 1100) trim(iVar('BC_TYPE',BCV)),               &
            trim(BC_TYPE(BCV)), VALID_BC_TYPE
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ENDDO L50   ! end loop over (bcv=1,dimension_bc)

      CALL FINL_ERR_MSG

      RETURN


 1100 FORMAT('Error 1100: Illegal entry: ',A,' = ',A,/'Valid entries:',&
         ' ',10(/5X,A,2x,A),/5X,A)

      END SUBROUTINE CHECK_BC_GEOMETRY



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_GEOMETRY_WALL                                  !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for walls BC's         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_GEOMETRY_WALL(BCV,dx,dy,dz)

! Global Variables:
!---------------------------------------------------------------------//
! Boundary condition locations and corresponding grid index
      use bc, only: BC_X_w, BC_X_e, BC_I_w, BC_I_e
      use bc, only: BC_Y_s, BC_Y_n, BC_J_s, BC_J_n
      use bc, only: BC_Z_b, BC_Z_t, BC_K_b, BC_K_t
! Basic grid information
      use geometry, only: domlo, domhi
      use geometry, only: XLENGTH
      use geometry, only: YLENGTH
      use geometry, only: ZLENGTH
! Function to compare two values
      use toleranc, only: COMPARE

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival

      use location_check_module, only: location_check

      use run, only: IFILE_NAME

      IMPLICIT NONE

      real(c_real), intent(in) :: dx, dy, dz

! Dummy Arguments:
!---------------------------------------------------------------------//
! Index of boundary condition.
      integer, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Calculated indices of the wall boundary
      integer :: I_w , I_e , J_s , J_n , K_b , K_t
! Integer error flag
      integer :: IER
!......................................................................!


      CALL INIT_ERR_MSG("CHECK_BC_GEOMETRY_WALL")

      IF(IS_DEFINED(BC_X_W(BCV)) .AND. IS_DEFINED(BC_X_E(BCV))) THEN

! setting indices to 1 if there is no variation in the i (x) direction
         I_W = CALC_CELL (BC_X_W(BCV), DX)
         I_W = I_W + 1
         I_E = CALC_CELL (BC_X_E(BCV), DX)
! BC along zy plane, checking if far west or far east of domain
         IF(EQUAL(BC_X_W(BCV), BC_X_E(BCV))) THEN
            IF(COMPARE(BC_X_W(BCV),0.0d0)) THEN
               I_W = domlo(1)-1
               I_E = domlo(1)-1
            ELSEIF(COMPARE(BC_X_W(BCV),XLENGTH)) THEN
               I_W = domhi(1)+1
               I_E = domhi(1)+1
            ENDIF
         ENDIF

! checking/setting corresponding i indices according to specified x
! coordinates
         IF(BC_I_W(BCV)/=UNDEFINED_I .OR. BC_I_E(BCV)/=UNDEFINED_I) THEN
            CALL LOCATION_CHECK (BC_I_W(BCV), I_W, BCV, 'BC - west')
            CALL LOCATION_CHECK (BC_I_E(BCV), I_E, BCV, 'BC - east')
         ELSE
            BC_I_W(BCV) = I_W
            BC_I_E(BCV) = I_E
         ENDIF
      ENDIF


      IF(IS_DEFINED(BC_Y_S(BCV)) .AND. IS_DEFINED(BC_Y_N(BCV))) THEN
! setting indices to 1 if there is no variation in the j (y) direction
         J_S = CALC_CELL (BC_Y_S(BCV), DY)
         J_S = J_S + 1
         J_N = CALC_CELL (BC_Y_N(BCV), DY)
! BC along xz plane, checking if far south or far north of domain
         IF(EQUAL(BC_Y_S(BCV), BC_Y_N(BCV))) THEN
            IF(COMPARE(BC_Y_S(BCV),ZERO)) THEN
               J_S = domlo(2)-1
               J_N = domlo(2)-1
            ELSE IF (COMPARE(BC_Y_S(BCV),YLENGTH)) THEN
               J_S = domhi(2)+1
               J_N = domhi(2)+1
            ENDIF
         ENDIF
! checking/setting corresponding j indices according to specified y
! coordinates
         IF(BC_J_S(BCV)/=UNDEFINED_I .OR. BC_J_N(BCV)/=UNDEFINED_I) THEN
            CALL LOCATION_CHECK (BC_J_S(BCV), J_S, BCV, 'BC - south')
            CALL LOCATION_CHECK (BC_J_N(BCV), J_N, BCV, 'BC - north')
         ELSE
            BC_J_S(BCV) = J_S
            BC_J_N(BCV) = J_N
         ENDIF
      ENDIF

      IF(IS_DEFINED(BC_Z_B(BCV)) .AND. IS_DEFINED(BC_Z_T(BCV))) THEN
! setting indices to 1 if there is no variation in the k (z) direction
         K_B = CALC_CELL (BC_Z_B(BCV), DZ)
         K_B = K_B + 1
         K_T = CALC_CELL (BC_Z_T(BCV), DZ)
! BC along xy plane, checking if far bottom or far top of domain
         IF(EQUAL(BC_Z_B(BCV), BC_Z_T(BCV))) THEN
            IF(COMPARE(BC_Z_B(BCV),ZERO)) THEN
               K_B = domlo(3)-1
               K_T = domlo(3)-1
            ELSEIF(COMPARE(BC_Z_B(BCV),ZLENGTH)) THEN
               K_B = domhi(3)+1
               K_T = domhi(3)+1
            ENDIF
         ENDIF
! checking/setting corresponding j indices according to specified y
! coordinates
         IF(BC_K_B(BCV)/=UNDEFINED_I .OR.BC_K_T(BCV)/=UNDEFINED_I) THEN
            CALL LOCATION_CHECK (BC_K_B(BCV), K_B, BCV, 'BC - bottom')
            CALL LOCATION_CHECK (BC_K_T(BCV), K_T, BCV, 'BC - top')
         ELSE
            BC_K_B(BCV) = K_B
            BC_K_T(BCV) = K_T
         ENDIF
      ENDIF


! CHECK FOR VALID VALUES
      IER = 0
      IF (BC_K_B(BCV)<domlo(3)-1 .OR. BC_K_B(BCV)>domhi(3)+1) IER = 1
      IF (BC_J_S(BCV)<domlo(3)-1 .OR. BC_J_S(BCV)>domhi(2)+1) IER = 1
      IF (BC_I_W(BCV)<domlo(2)-1 .OR. BC_I_W(BCV)>domhi(1)+1) IER = 1
      IF (BC_K_T(BCV)<domlo(2)-1 .OR. BC_K_T(BCV)>domhi(3)+1) IER = 1
      IF (BC_J_N(BCV)<domlo(1)-1 .OR. BC_J_N(BCV)>domhi(2)+1) IER = 1
      IF (BC_I_E(BCV)<domlo(1)-1 .OR. BC_I_E(BCV)>domhi(1)+1) IER = 1
      IF (BC_K_B(BCV) > BC_K_T(BCV)) IER = 1
      IF (BC_J_S(BCV) > BC_J_N(BCV)) IER = 1
      IF (BC_I_W(BCV) > BC_I_E(BCV)) IER = 1

      IF(IER /= 0)THEN
         WRITE(ERR_MSG,1100) BCV,                                      &
            'X', BC_X_W(BCV), BC_X_E(BCV),'I',BC_I_W(BCV),BC_I_E(BCV), &
            'Y', BC_Y_S(BCV), BC_Y_N(BCV),'J',BC_J_S(BCV),BC_J_N(BCV), &
            'Z', BC_Z_B(BCV), BC_Z_T(BCV),'K',BC_K_B(BCV),BC_K_T(BCV), &
            trim( IFILE_NAME )
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Invalid location specified for BC ',I3,'.',  &
         3(/3x,A1,': ',g12.5,',',g12.5,8x,A1,': ',I8,',',I8),/         &
         'Please correct the ', A, ' file.')

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_BC_GEOMETRY_WALL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_GEOMETRY_FLOW                                  !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for flow BC's. Also    !
!           set value of bc_plane for flow BC's.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_GEOMETRY_FLOW(BCV,dx,dy,dz)

! Global Variables:
!---------------------------------------------------------------------//
! Boundary condition locations and corresponding grid index
      use bc, only: BC_X_w, BC_X_e, BC_I_w, BC_I_e
      use bc, only: BC_Y_s, BC_Y_n, BC_J_s, BC_J_n
      use bc, only: BC_Z_b, BC_Z_t, BC_K_b, BC_K_t
! Basic grid information
      use geometry, only: domlo,domhi
      use geometry, only: xlength, ylength, zlength

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar, ival

      use location_check_module, only: location_check

      use run, only: IFILE_NAME

      IMPLICIT NONE

      real(c_real), intent(in) :: dx, dy, dz

! Dummy Arguments:
!---------------------------------------------------------------------//
! Index of boundary condition.
      integer, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Calculated indices of the wall boundary
      integer :: I_w, I_e, J_s, J_n, K_b, K_t
! Indices for error checking
      integer :: IER

! surface indictors:
! a value of T indicates that the defined boundary region does not
! vary in indicated coordinate direction. that is, if bc_x_w is
! equal to bc_x_e then the boundary region must be in the yz plane
      logical :: X_CONSTANT, Y_CONSTANT, Z_CONSTANT
!......................................................................!


      CALL INIT_ERR_MSG("CHECK_BC_GEOMETRY_FLOW")

      X_CONSTANT = .TRUE.
      Y_CONSTANT = .TRUE.
      Z_CONSTANT = .TRUE.

      if (is_defined(bc_x_w(bcv)) .and. is_defined(bc_x_e(bcv))) then
         i_w = calc_cell (bc_x_w(bcv), dx)
         i_e = calc_cell (bc_x_e(bcv), dx)
         if (.not.equal(bc_x_w(bcv), bc_x_e(bcv))) then
            x_constant = .false.
            i_w = i_w + 1
         else if(equal(bc_x_w(bcv),xlength)) then
            i_w = i_w + 1
            i_e = i_w
         endif
         bc_i_w(bcv) = i_w
         bc_i_e(bcv) = i_e
      endif

      if (is_defined(bc_y_s(bcv)) .and. is_defined(bc_y_n(bcv))) then
         j_s = calc_cell (bc_y_s(bcv), dy)
         j_n = calc_cell (bc_y_n(bcv), dy)
         if(.not.equal(bc_y_s(bcv), bc_y_n(bcv))) then
            y_constant = .false.
            j_s = j_s + 1
         else if(equal(bc_y_s(bcv),ylength)) then
            j_s = j_s + 1
            j_n = j_s
         endif
         bc_j_s(bcv) = j_s
         bc_j_n(bcv) = j_n
      endif

      if(is_defined(bc_z_b(bcv)) .and. is_defined(bc_z_t(bcv))) then
         k_b = calc_cell (bc_z_b(bcv), dz)
         k_t = calc_cell (bc_z_t(bcv), dz)
         if(.not.equal(bc_z_b(bcv), bc_z_t(bcv))) then
            z_constant = .false.
            k_b = k_b + 1
         else if(equal(bc_z_b(bcv),zlength)) then
            k_b = k_b + 1
            k_t = k_b
         endif
         bc_k_b(bcv) = k_b
         bc_k_t(bcv) = k_t
      endif

! Check whether the boundary is a plane parallel to one of the three
! coordinate planes
      IF(IS_DEFINED(BC_X_W(BCV)) .AND. IS_DEFINED(BC_Y_S(BCV)) .AND. &
         IS_DEFINED(BC_Z_B(BCV))) CALL CHECK_PLANE (X_CONSTANT, &
         Y_CONSTANT, Z_CONSTANT, BCV, 'BC')


! CHECK FOR VALID VALUES
      IER = 0
      IF(BC_I_W(BCV)<domlo(1)-1 .OR. BC_I_W(BCV)>domhi(1)+1) IER = 1
      IF(BC_I_E(BCV)<domlo(1)-1 .OR. BC_I_E(BCV)>domhi(1)+1) IER = 1
      IF(BC_J_S(BCV)<domlo(2)-1 .OR. BC_J_S(BCV)>domhi(2)+1) IER = 1
      IF(BC_J_N(BCV)<domlo(2)-1 .OR. BC_J_N(BCV)>domhi(2)+1) IER = 1
      IF(BC_K_B(BCV)<domlo(3)-1 .OR. BC_K_B(BCV)>domhi(3)+1) IER = 1
      IF(BC_K_T(BCV)<domlo(3)-1 .OR. BC_K_T(BCV)>domhi(3)+1) IER = 1
      IF(BC_K_B(BCV) > BC_K_T(BCV)) IER = 1
      IF(BC_J_S(BCV) > BC_J_N(BCV)) IER = 1
      IF(BC_I_W(BCV) > BC_I_E(BCV)) IER = 1

      IF(IER /= 0)THEN
         WRITE(ERR_MSG,1100) BCV,                                      &
            'X', BC_X_W(BCV), BC_X_E(BCV),'I',BC_I_W(BCV),BC_I_E(BCV), &
            'Y', BC_Y_S(BCV), BC_Y_N(BCV),'J',BC_J_S(BCV),BC_J_N(BCV), &
            'Z', BC_Z_B(BCV), BC_Z_T(BCV),'K',BC_K_B(BCV),BC_K_T(BCV), &
            trim( IFILE_NAME )
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Invalid location specified for BC ',I3,'.',  &
         3(/3x,A1,': ',g12.5,',',g12.5,8x,A1,': ',I8,',',I8),/         &
         'Please correct the ',A,' file.')

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_BC_GEOMETRY_FLOW
END MODULE CHECK_BC_GEOMETRY_MODULE
