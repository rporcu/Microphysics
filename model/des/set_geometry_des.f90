MODULE SET_GEOMETRY_DES_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_GEOMETRY_DES                                        !
!  Author:   R.Garg                                   Date: 19-Mar-14  !
!                                                                      !
!  Purpose: Allocate des arrays that are based on Eulerian grid.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_GEOMETRY_DES(dx,dy,dz)


! Global Variables:
!---------------------------------------------------------------------//
! Arrays for DEM simulations delineating cell edges.
      use discretelement, only: XE, YN, ZT, DIMN
! Fluid grid cell dimensions and mesh size
      USE geometry, only: IMIN2, IMAX2
      USE geometry, only: JMIN2, JMAX2
      USE geometry, only: KMIN2, KMAX2

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO

! Module proceedures.
!---------------------------------------------------------------------//
      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE

      real(c_real), intent(in) :: dx, dy, dz

! Local Variables:
!---------------------------------------------------------------------//
! Generic loop indices
      INTEGER :: I, J, K
! Error Flag
      INTEGER :: IER
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("SET_GEOMETRY_DES")

      Allocate( XE (0:imax2), STAT=IER )
      Allocate( YN (0:jmax2), STAT=IER )
      Allocate( ZT (0:kmax2), STAT=IER )

! Collect the error flags from all ranks. If all allocaitons were
! successfull, do nothing. Otherwise, flag the error and abort.
      ! CALL GLOBAL_ALL_SUM(IER)

! Initialize arrays.
      XE(:) = ZERO
      YN(:) = ZERO
      ZT(:) = ZERO

! Each loop starts at 2 and goes to max+2 (i.e., imin1=2, imax2=imax+2)
! However, the indices range to include ghost cells (0-imax2) to avoid

      XE(IMIN2-1) = ZERO-DX
      DO I = IMIN2, IMAX2
         XE(I) = XE(I-1) + DX
      ENDDO

      YN(JMIN2-1) = ZERO-DY
      DO J  = JMIN2, JMAX2
         YN(J) = YN(J-1) + DY
      ENDDO

      IF(DIMN.EQ.3) THEN
         ZT(KMIN2-1) = ZERO-DZ
         DO K = KMIN2, KMAX2
            ZT(K) = ZT(K-1) + DZ
         ENDDO
      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_GEOMETRY_DES
END MODULE SET_GEOMETRY_DES_MODULE
