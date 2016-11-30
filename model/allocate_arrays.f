module allocate_mod
   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutinee: ALLOCATE_ARRAYS                                        C
!  Purpose: allocate arrays                                            C
!                                                                      C
!  Author: M. Syamlal                                Date: 17-DEC-98   C
!  Reviewer:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ALLOCATE_ARRAYS(A_m, B_m)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      Use fldvar  , only: ep_g , p_g , ro_g , rop_g , u_g , v_g , w_g
      Use fldvar  , only: ep_go, p_go, ro_go, rop_go, u_go, v_go, w_go
      Use fldvar  , only: d_e, d_n, d_t
      Use fldvar  , only: pp_g
      Use fldvar  , only: mu_g, lambda_g, trD_g, tau_u_g, tau_v_g, tau_w_g
      Use fldvar  , only: flux_ge, flux_gn, flux_gt
      Use fldvar  , only:  rop_ge,  rop_gn,  rop_gt
      Use xsi_array,only: xsi_n, xsi_e, xsi_t

      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3

      IMPLICIT NONE

!-----------------------------------------------
! Variables
!-----------------------------------------------
      DOUBLE PRECISION, allocatable, intent(inOUT) :: A_m(:,:,:,:)
      DOUBLE PRECISION, allocatable, intent(inOUT) :: B_m(:,:,:)

      integer :: is3, ie3
      integer :: js3, je3
      integer :: ks3, ke3

      is3 = istart3;   ie3 = iend3
      js3 = jstart3;   je3 = jend3
      ks3 = kstart3;   ke3 = kend3

      Allocate( A_m(is3:ie3,js3:je3,ks3:ke3, -3:3) )
      Allocate( B_m(is3:ie3,js3:je3,ks3:ke3) )

!fldvar
      Allocate( EP_g (is3:ie3,js3:je3,ks3:ke3) )
      Allocate( P_g (is3:ie3,js3:je3,ks3:ke3) )
      Allocate( RO_g (is3:ie3,js3:je3,ks3:ke3) )
      Allocate( ROP_g (is3:ie3,js3:je3,ks3:ke3) )

      Allocate( U_g (is3:ie3,js3:je3,ks3:ke3) )
      Allocate( V_g (is3:ie3,js3:je3,ks3:ke3) )
      Allocate( W_g (is3:ie3,js3:je3,ks3:ke3) )

      Allocate( EP_go  (is3:ie3,js3:je3,ks3:ke3))
      Allocate( P_go   (is3:ie3,js3:je3,ks3:ke3))
      Allocate( RO_go  (is3:ie3,js3:je3,ks3:ke3))
      Allocate( ROP_go (is3:ie3,js3:je3,ks3:ke3))

      Allocate( U_go (is3:ie3,js3:je3,ks3:ke3))
      Allocate( V_go (is3:ie3,js3:je3,ks3:ke3))
      Allocate( W_go (is3:ie3,js3:je3,ks3:ke3))

      Allocate( d_e(is3:ie3,js3:je3,ks3:ke3))
      Allocate( d_n(is3:ie3,js3:je3,ks3:ke3))
      Allocate( d_t(is3:ie3,js3:je3,ks3:ke3))

      Allocate( Pp_g(is3:ie3,js3:je3,ks3:ke3) )

!physprop
      Allocate( MU_g(is3:ie3,js3:je3,ks3:ke3) )

!visc_g
      Allocate( trD_g(is3:ie3,js3:je3,ks3:ke3))
      Allocate( LAMBDA_g(is3:ie3,js3:je3,ks3:ke3))
      Allocate( TAU_U_g(is3:ie3,js3:je3,ks3:ke3))
      Allocate( TAU_V_g(is3:ie3,js3:je3,ks3:ke3))
      Allocate( TAU_W_g(is3:ie3,js3:je3,ks3:ke3))

!xsi_array
      Allocate( xsi_e(is3:ie3,js3:je3,ks3:ke3))
      Allocate( xsi_n(is3:ie3,js3:je3,ks3:ke3))
      Allocate( xsi_t(is3:ie3,js3:je3,ks3:ke3))


!mflux
      Allocate( Flux_gE(is3:ie3,js3:je3,ks3:ke3))
      Allocate( Flux_gN(is3:ie3,js3:je3,ks3:ke3))
      Allocate( Flux_gT(is3:ie3,js3:je3,ks3:ke3))

      Allocate( ROP_gE(is3:ie3,js3:je3,ks3:ke3))
      Allocate( ROP_gN(is3:ie3,js3:je3,ks3:ke3))
      Allocate( ROP_gT(is3:ie3,js3:je3,ks3:ke3))

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

! Domain flags.
      use geometry, only: FLAG
! Domain volumes and areas.
      use geometry, only: VOL_SURR

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
      Allocate( FLAG  (istart3:iend3,jstart3:jend3,kstart3:kend3,0:4), STAT=IER)
      IF(IER /= 0) goto 500

! total volume of each cell's surrounding stencil cells
      Allocate( VOL_SURR (istart3:iend3,jstart3:jend3,kstart3:kend3))
      IF(IER /= 0) goto 500

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
end module allocate_mod
