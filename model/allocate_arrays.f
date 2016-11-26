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
      Allocate( P_g (DIMENSION_3) )
      Allocate( RO_g (DIMENSION_3) )
      Allocate( ROP_g (DIMENSION_3) )

      Allocate( U_g (DIMENSION_3) )
      Allocate( V_g (DIMENSION_3) )
      Allocate( W_g (DIMENSION_3) )

      Allocate( EP_go  (istart3:iend3,jstart3:jend3,kstart3:kend3))
      Allocate( P_go   (istart3:iend3,jstart3:jend3,kstart3:kend3))
      Allocate( RO_go  (istart3:iend3,jstart3:jend3,kstart3:kend3))
      Allocate( ROP_go (istart3:iend3,jstart3:jend3,kstart3:kend3))

      Allocate( U_go (istart3:iend3,jstart3:jend3,kstart3:kend3))
      Allocate( V_go (istart3:iend3,jstart3:jend3,kstart3:kend3))
      Allocate( W_go (istart3:iend3,jstart3:jend3,kstart3:kend3))

      Allocate( d_e(istart3:iend3,jstart3:jend3,kstart3:kend3))
      Allocate( d_n(istart3:iend3,jstart3:jend3,kstart3:kend3))
      Allocate( d_t(istart3:iend3,jstart3:jend3,kstart3:kend3))

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
      use compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3

! Domain decomposition and dimensions
      use geometry, only: oDX
      use geometry, only: oDZ
      use geometry, only: oDY
! Domain flags.
      use geometry, only: ICBC_FLAG
      use geometry, only: FLAG
      use geometry, only: FLAG_E, FLAG_N, FLAG_T
! Domain volumes and areas.
      use geometry, only: VOL_SURR
! Axis decomposition
      USE param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
      USE param, only: DIMENSION_3, DIMENSION_4
      USE param, only: DIMENSION_3L

! Module procedures
!---------------------------------------------------------------------//
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

! Flags for the scalar grid.
      Allocate( FLAG  (istart3:iend3,jstart3:jend3,kstart3:kend3))
      IF(IER /= 0) goto 500

! Flags for the momentum grids.
      Allocate( FLAG_E (istart3:iend3,jstart3:jend3,kstart3:kend3))
      Allocate( FLAG_N (istart3:iend3,jstart3:jend3,kstart3:kend3))
      Allocate( FLAG_T (istart3:iend3,jstart3:jend3,kstart3:kend3))
      IF(IER /= 0) goto 500

! Text flags for scalar grid.
      Allocate( ICBC_FLAG (DIMENSION_3L), STAT=IER )
      IF(IER /= 0) goto 500

! total volume of each cell's surrounding stencil cells
      Allocate( VOL_SURR (istart3:iend3,jstart3:jend3,kstart3:kend3))

! Collect the error flags from all ranks. If all allocaitons were
! successfull, do nothing. Otherwise, flag the error and abort.
! Note that the allocation status is checked in groups. This can
! be increase if tracking the source of an allocation failure.
  ! 500 CALL GLOBAL_ALL_SUM(IER)

      500 IF(IER /= 0) THEN
         WRITE(ERR_MSG,1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Failure during array allocation.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS_GEOMETRY
