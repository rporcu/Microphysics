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
      SUBROUTINE ALLOCATE_ARRAYS(A_m, B_m,ep_g,p_g,ro_g,rop_g,u_g,v_g,w_g,&
         ep_go,p_go,ro_go,rop_go,u_go,v_go,w_go,d_e,d_n,d_t,pp_g,&
         mu_g,lambda_g,trD_g,tau_u_g,tau_v_g,tau_w_g,flux_ge,&
         flux_gn,flux_gt,rop_ge,rop_gn,rop_gt,f_gds, drag_am, drag_bm)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar, only: istart3, iend3
      use compar, only: jstart3, jend3
      use compar, only: kstart3, kend3

      use param1, only: undefined

      IMPLICIT NONE

      double precision, allocatable, intent(inout) :: A_m(:,:,:,:)
      double precision, allocatable, intent(inout) :: B_m(:,:,:)
      double precision, allocatable, intent(inout) :: ep_g(:,:,:)
      double precision, allocatable, intent(inout) :: p_g(:,:,:)
      double precision, allocatable, intent(inout) :: ro_g(:,:,:)
      double precision, allocatable, intent(inout) :: rop_g(:,:,:)
      double precision, allocatable, intent(inout) :: u_g(:,:,:)
      double precision, allocatable, intent(inout) :: v_g(:,:,:)
      double precision, allocatable, intent(inout) :: w_g(:,:,:)
      double precision, allocatable, intent(inout) :: ep_go(:,:,:)
      double precision, allocatable, intent(inout) :: p_go(:,:,:)
      double precision, allocatable, intent(inout) :: ro_go(:,:,:)
      double precision, allocatable, intent(inout) :: rop_go(:,:,:)
      double precision, allocatable, intent(inout) :: u_go(:,:,:)
      double precision, allocatable, intent(inout) :: v_go(:,:,:)
      double precision, allocatable, intent(inout) :: w_go(:,:,:)
      double precision, allocatable, intent(inout) :: d_e(:,:,:)
      double precision, allocatable, intent(inout) :: d_n(:,:,:)
      double precision, allocatable, intent(inout) :: d_t(:,:,:)
      double precision, allocatable, intent(inout) :: pp_g(:,:,:)
      double precision, allocatable, intent(inout) :: mu_g(:,:,:)
      double precision, allocatable, intent(inout) :: lambda_g(:,:,:)
      double precision, allocatable, intent(inout) :: trD_g(:,:,:)
      double precision, allocatable, intent(inout) :: tau_u_g(:,:,:)
      double precision, allocatable, intent(inout) :: tau_v_g(:,:,:)
      double precision, allocatable, intent(inout) :: tau_w_g(:,:,:)
      double precision, allocatable, intent(inout) :: flux_ge(:,:,:)
      double precision, allocatable, intent(inout) :: flux_gn(:,:,:)
      double precision, allocatable, intent(inout) :: flux_gt(:,:,:)
      double precision, allocatable, intent(inout) :: rop_ge(:,:,:)
      double precision, allocatable, intent(inout) :: rop_gn(:,:,:)
      double precision, allocatable, intent(inout) :: rop_gt(:,:,:)
      double precision, allocatable, intent(inout) :: f_gds(:,:,:)
      double precision, allocatable, intent(inout) :: drag_am(:,:,:)
      double precision, allocatable, intent(inout) :: drag_bm(:,:,:,:)

      integer :: is3, ie3
      integer :: js3, je3
      integer :: ks3, ke3

      is3 = istart3;   ie3 = iend3
      js3 = jstart3;   je3 = jend3
      ks3 = kstart3;   ke3 = kend3

      Allocate( A_m(is3:ie3,js3:je3,ks3:ke3, -3:3) )
      Allocate( B_m(is3:ie3,js3:je3,ks3:ke3) )

      Allocate( EP_g (is3:ie3,js3:je3,ks3:ke3) )
      Allocate( P_g (is3:ie3,js3:je3,ks3:ke3) )
      Allocate( RO_g (is3:ie3,js3:je3,ks3:ke3) )
      Allocate( ROP_g (is3:ie3,js3:je3,ks3:ke3) )

      Allocate( U_g (is3:ie3,js3:je3,ks3:ke3) )
      Allocate( V_g (is3:ie3,js3:je3,ks3:ke3) )
      Allocate( W_g (is3:ie3,js3:je3,ks3:ke3) )

      ! Initialize arrays.
      EP_G = UNDEFINED
      P_G = UNDEFINED
      RO_G = UNDEFINED
      ROP_G = UNDEFINED
      U_G = UNDEFINED
      V_G = UNDEFINED
      W_G = UNDEFINED

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

      Allocate( MU_g(is3:ie3,js3:je3,ks3:ke3) )

      Allocate( trD_g(is3:ie3,js3:je3,ks3:ke3))
      Allocate( LAMBDA_g(is3:ie3,js3:je3,ks3:ke3))
      Allocate( TAU_U_g(is3:ie3,js3:je3,ks3:ke3))
      Allocate( TAU_V_g(is3:ie3,js3:je3,ks3:ke3))
      Allocate( TAU_W_g(is3:ie3,js3:je3,ks3:ke3))

      Allocate( Flux_gE(is3:ie3,js3:je3,ks3:ke3))
      Allocate( Flux_gN(is3:ie3,js3:je3,ks3:ke3))
      Allocate( Flux_gT(is3:ie3,js3:je3,ks3:ke3))

      Allocate( ROP_gE(is3:ie3,js3:je3,ks3:ke3))
      Allocate( ROP_gN(is3:ie3,js3:je3,ks3:ke3))
      Allocate( ROP_gT(is3:ie3,js3:je3,ks3:ke3))

      Allocate( DRAG_AM (is3:ie3,js3:je3,ks3:ke3))
      Allocate( F_GDS   (is3:ie3,js3:je3,ks3:ke3))
      Allocate( DRAG_BM (is3:ie3,js3:je3,ks3:ke3, 3))

      DRAG_AM = 0.0d0
      F_GDS   = 0.0d0
      DRAG_BM = 0.0d0

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
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

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
