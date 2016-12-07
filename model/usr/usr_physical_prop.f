!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PHYSICAL_PROP_ROg                                       !
!  Purpose: User hook for calculating the gas phase density.           !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PHYSICAL_PROP_ROg(flag)

! Fluid grid loop bounds.
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3
! Function to identify wall cells
      use functions, only: funijk

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      implicit none

      integer, intent(in   ) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

! Local Variables:
!---------------------------------------------------------------------//
! Loop indicies
      INTEGER :: I, J, K, IJK   ! Computational cell
!......................................................................!


! The following error message is used to make sure that if a user
! defined gas density is invoked, that this routine has been modified.

!- REMOVE THE FOLLOWING ---------------------------------------------->>

      CALL INIT_ERR_MSG('USR_PHYSICAL_PROP_ROg')
      WRITE(ERR_MSG,9999)
      CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 9999 FORMAT('ERROR 9999: The user-defined drag routine was invoked ', &
         'but this',/'generic error message exits. Either choose a ',  &
         'different drag law',/'or correct mfix/model/usr_drag.f')

!- END REMOVE --------------------------------------------------------<<

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
         IJK = FUNIJK(i,j,k)

         if (flag(i,j,k,1) >= 100) cycle

! Calculate the fluid density and bulk density
!         RO_G(IJK) =
!         ROP_G(IJK) =

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE USR_PHYSICAL_PROP_ROg
