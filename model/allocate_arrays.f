!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutinee: ALLOCATE_ARRAYS                                        C
!  Purpose: allocate arrays                                            C
!                                                                      C
!  Author: M. Syamlal                                Date: 17-DEC-98   C
!  Reviewer:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ALLOCATE_ARRAYS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      use matrix
      Use drag
      Use fldvar
      Use geometry
      Use physprop
      Use residual
      Use run
      Use xsi_array

      IMPLICIT NONE

!-----------------------------------------------
! Variables
!-----------------------------------------------

!ambm
      Allocate( A_m(DIMENSION_3, -3:3) )
      Allocate( B_m(DIMENSION_3) )

!fldvar
      Allocate( EP_g (DIMENSION_3) )
      Allocate( EP_go (DIMENSION_3p) )
      Allocate( P_g (DIMENSION_3) )
      Allocate( P_go (DIMENSION_3p) )
      Allocate( RO_g (DIMENSION_3) )
      Allocate( RO_go (DIMENSION_3p) )
      Allocate( ROP_g (DIMENSION_3) )
      Allocate( ROP_go (DIMENSION_3p) )

      Allocate( U_g (DIMENSION_3) )
      Allocate( U_go (DIMENSION_3p) )
      Allocate( V_g (DIMENSION_3) )
      Allocate( V_go (DIMENSION_3p) )
      Allocate( W_g (DIMENSION_3) )
      Allocate( W_go (DIMENSION_3p) )


      Allocate( d_e(DIMENSION_3p) )
      Allocate( d_n(DIMENSION_3p) )
      Allocate( d_t(DIMENSION_3p) )
      Allocate( Pp_g(DIMENSION_3p) )

!physprop
      Allocate( MU_g (DIMENSION_3) )

!visc_g
      Allocate( trD_g(DIMENSION_3) )
      Allocate( LAMBDA_g (DIMENSION_3p) )
      Allocate( TAU_U_g(DIMENSION_3p) )
      Allocate( TAU_V_g(DIMENSION_3p) )
      Allocate( TAU_W_g(DIMENSION_3p) )

!xsi_array
      Allocate( Xsi_e(DIMENSION_3) )
      Allocate( Xsi_n(DIMENSION_3) )
      Allocate( Xsi_t(DIMENSION_3) )

!mflux
      Allocate( Flux_gE(DIMENSION_3p) )
      Allocate( Flux_gN(DIMENSION_3p) )
      Allocate( Flux_gT(DIMENSION_3p) )

      Allocate( ROP_gE(DIMENSION_3p) )
      Allocate( ROP_gN(DIMENSION_3p) )
      Allocate( ROP_gT(DIMENSION_3p) )



      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_ARRAYS_GEOMETRY                               !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Calculate X, X_E,  oX, oX_E                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_ARRAYS_GEOMETRY

! Global Variables:
!---------------------------------------------------------------------//
! Domain decomposition and dimensions
      use geometry, only: oDX, oDX_E
      use geometry, only: oDZ, oDZ_T
      use geometry, only: oDY, oDY_N
      use geometry, only: X, X_E, oX, oX_E
      use geometry, only: Z, Z_T
! Averaging factors.
      use geometry, only: FX_E, FX_E_bar, FX, FX_bar
      use geometry, only: FY_N, FY_N_bar
      use geometry, only: FZ_T, FZ_T_bar
! Domain flags.
      use geometry, only: ICBC_FLAG
      use geometry, only: FLAG, FLAG3
      use geometry, only: FLAG_E, FLAG_N, FLAG_T
! Domain volumes and areas.
      use geometry, only: VOL, VOL_SURR, AYZ, AXZ, AXY! Scalar grid
      use geometry, only: VOL_U, AYZ_U, AXZ_U, AXY_U  ! X-Momentum
      use geometry, only: VOL_V, AYZ_V, AXZ_V, AXY_V  ! Y-Momentum
      use geometry, only: VOL_W, AYZ_W, AXZ_W, AXY_W  ! Z-Momentum
! Axis decomposition
      USE param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
      USE param, only: DIMENSION_3, DIMENSION_4
      USE param, only: DIMENSION_3L, DIMENSION_3P

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Error Flag
      INTEGER :: IER
! Flag indicating that the arrays were previously allocated.
      INTEGER, SAVE :: CALLED = -1
!......................................................................!

      CALLED = CALLED + 1

      IF(CALLED > 0 .and. mod(CALLED,2) /= 0) RETURN

! Initialize the error manager.
      CALL INIT_ERR_MSG("ALLOCATE_ARRAYS_GEOMETRY")

! Allocate geometry components related to the mesh. Check the
! allocation error status and abort if any failure is detected.
      ALLOCATE( X     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( X_E   (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oX    (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oX_E  (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oDX   (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oDX_E (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( oDY   (0:DIMENSION_J), STAT=IER )
      ALLOCATE( oDY_N (0:DIMENSION_J), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( Z     (0:DIMENSION_K), STAT=IER )
      ALLOCATE( Z_T   (0:DIMENSION_K), STAT=IER )
      ALLOCATE( oDZ   (0:DIMENSION_K), STAT=IER )
      ALLOCATE( oDZ_T (0:DIMENSION_K), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( FX     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( FX_bar (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( FX_E     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( FX_E_bar (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( FY_N     (0:DIMENSION_J), STAT=IER )
      ALLOCATE( FY_N_bar (0:DIMENSION_J), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( FZ_T     (0:DIMENSION_K), STAT=IER )
      ALLOCATE( FZ_T_bar (0:DIMENSION_K), STAT=IER )
      IF(IER /= 0) goto 500

! Flags for the scalar grid.
      Allocate( FLAG  (DIMENSION_3), STAT=IER )
      Allocate( FLAG3 (DIMENSION_4), STAT=IER )
      IF(IER /= 0) goto 500

! Flags for the momentum grids.
      Allocate( FLAG_E (DIMENSION_3), STAT=IER )
      Allocate( FLAG_N (DIMENSION_3), STAT=IER )
      Allocate( FLAG_T (DIMENSION_3), STAT=IER )
      IF(IER /= 0) goto 500

! Text flags for scalar grid.
      Allocate( ICBC_FLAG (DIMENSION_3L), STAT=IER )
      IF(IER /= 0) goto 500

! Volume and face-areas of scalar grid.
      Allocate( VOL (DIMENSION_3),  STAT=IER )
      Allocate( AYZ (DIMENSION_3P), STAT=IER )
      Allocate( AXZ (DIMENSION_3P), STAT=IER )
      Allocate( AXY (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

      ! total volume of each cell's surrounding stencil cells
      Allocate( VOL_SURR (DIMENSION_3), STAT=IER )

! Volume and face-areas of X-Momentumn grid.
      Allocate( VOL_U (DIMENSION_3),  STAT=IER )
      Allocate( AYZ_U (DIMENSION_3P), STAT=IER )
      Allocate( AXZ_U (DIMENSION_3P), STAT=IER )
      Allocate( AXY_U (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

! Volume and face-areas of Y-Momentum grid.
      Allocate( VOL_V (DIMENSION_3),  STAT=IER )
      Allocate( AYZ_V (DIMENSION_3P), STAT=IER )
      Allocate( AXZ_V (DIMENSION_3P), STAT=IER )
      Allocate( AXY_V (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

! Volume and face-areas of Z-Momentum grid.
      Allocate( VOL_W (DIMENSION_3),  STAT=IER )
      Allocate( AYZ_W (DIMENSION_3P), STAT=IER )
      Allocate( AXZ_W (DIMENSION_3P), STAT=IER )
      Allocate( AXY_W (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

! Collect the error flags from all ranks. If all allocaitons were
! successfull, do nothing. Otherwise, flag the error and abort.
! Note that the allocation status is checked in groups. This can
! be increase if tracking the source of an allocation failure.
  500 CALL GLOBAL_ALL_SUM(IER)

      IF(IER /= 0) THEN
         WRITE(ERR_MSG,1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Failure during array allocation.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS_GEOMETRY
