
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR1 (L)                                         C
!  Purpose: Write user-defined output                                  C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine write_usr1(l, np, time, dt, particles)

   use amrex_fort_module, only: c_real => amrex_real
   use usr,               only: UPDATE_RK4_SOL
   use particle_mod,      only: particle_t

   implicit none

   integer,          intent(in   ) :: l, np
   real(c_real),     intent(in   ) :: time, dt
   type(particle_t), intent(in   ) :: particles(np)

   
   SELECT CASE(L)
   CASE(1)
      CALL UPDATE_RK4_SOL(TIME)
      CALL WRITE_DES_OUT(TIME, np, particles)
   END SELECT

END SUBROUTINE WRITE_USR1


!......................................................................!
!  Subroutine: WRITE_DES_Out                                           !
!                                                                      !
!  Purpose: Calculate the position and velocity (1D Y-axis) of a free  !
!  falling particle. Compare the results to the MFIX-DEM solultion.    !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!......................................................................!
subroutine WRITE_DES_Out(lTime, np, particles)

   use amrex_fort_module, only: c_real => amrex_real
   use usr,               only: RK4_POS, RK4_VEL
   use particle_mod,      only: particle_t
   
   implicit none

   ! Dummy Arguments
   !---------------------------------------------------------------------//
   real(c_real),     intent(in   ) :: ltime
   integer,          intent(in   ) :: np
   type(particle_t), intent(in   ) :: particles(np)   

   ! Local variables
   !---------------------------------------------------------------------//
   ! file unit for heat transfer data
   integer, parameter :: lUNIT = 2030


   ! Open the file.
   OPEN(UNIT=lUNIT, FILE='POST_POS.dat',                            &
        POSITION="APPEND", STATUS='OLD')
   ! Write the results to file.
   WRITE(lUNIT,1000) lTime, RK4_POS(3), particles(1) % pos(3),        &
        ABS_ERR(RK4_POS(3), particles(1) % pos(3))
   ! Close the output file.
   CLOSE(lUNIT)


   ! Open the file.
   OPEN(UNIT=lUNIT, FILE='POST_VEL.dat',                            &
        POSITION="APPEND", STATUS='OLD')
   ! Write the results to file.
   WRITE(lUNIT,1000) lTime, RK4_VEL(3), particles(1) % vel(3),        &
        ABS_ERR(RK4_VEL(3), particles(1) % vel(3))
   ! Close the output file.
   CLOSE(lUNIT)

   RETURN

1000 FORMAT(3x,F15.6,5X,F15.6,3x,F15.6,3x,F15.4)

CONTAINS


   !......................................................................!
   !                                                                      !
   !  Subroutine name: ABS_ERR                                            !
   !  Author: J.Musser                                   Date:  May-2014  !
   !                                                                      !
   !  Purpose: Calculate either the absolute percent relative error or    !
   !  the absolute error.                                                 !
   !                                                                      !
   !......................................................................!
   REAL(C_REAL) FUNCTION ABS_ERR(EXT, NUM)

      use param, only: SMALL_NUMBER

      IMPLICIT NONE

      REAL(C_REAL), INTENT(IN) :: EXT, NUM

      IF(ABS(EXT) > SMALL_NUMBER) THEN
         ABS_ERR = ABS((EXT - NUM)/EXT)*1.0d2
      ELSE
         ABS_ERR = ABS(EXT - NUM)*1.0d2
      ENDIF
   END FUNCTION ABS_ERR

end subroutine WRITE_DES_Out


   ! !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
! !                                                                      C
! !  Module name: WRITE_USR1 (L)                                         C
! !  Purpose: Write user-defined output                                  C
! !                                                                      C
! !  Author:                                            Date: dd-mmm-yy  C
! !  Reviewer:                                          Date: dd-mmm-yy  C
! !                                                                      C
! !  Revision Number:                                                    C
! !  Purpose:                                                            C
! !  Author:                                            Date: dd-mmm-yy  C
! !  Reviewer:                                          Date: dd-mmm-yy  C
! !                                                                      C
! !  Literature/Document References:                                     C
! !                                                                      C
! !  Variables referenced:                                               C
! !  Variables modified:                                                 C
! !                                                                      C
! !  Local variables:                                                    C
! !                                                                      C
! !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!       subroutine write_usr1(l, time, dt, max_pip, des_pos_new, des_vel_new, omega_new)

!       use amrex_fort_module, only : c_real => amrex_real
!       use usr, only: UPDATE_RK4_SOL

!       IMPLICIT NONE

!       integer,      intent(in   ) :: l, max_pip
!       real(c_real), intent(in   ) :: time, dt
!       real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
!       real(c_real), intent(in   ) :: des_vel_new(max_pip,3)
!       real(c_real), intent(in   ) :: omega_new(max_pip,3)

!       SELECT CASE(L)
!       CASE(1)
!          CALL UPDATE_RK4_SOL(TIME)
!          CALL WRITE_DES_OUT(TIME, max_pip, des_pos_new, des_vel_new)
!       END SELECT

!       RETURN
!       END SUBROUTINE WRITE_USR1


! !......................................................................!
! !  Subroutine: WRITE_DES_Out                                           !
! !                                                                      !
! !  Purpose: Calculate the position and velocity (1D Y-axis) of a free  !
! !  falling particle. Compare the results to the MFIX-DEM solultion.    !
! !                                                                      !
! !  Author: J.Musser                                   Date:  Jan-13    !
! !                                                                      !
! !  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
! !  open-source MFIX-DEM software for gas-solids flows," from URL:      !
! !  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
! !......................................................................!
!       SUBROUTINE WRITE_DES_Out(lTime, max_pip, des_pos_new, des_vel_new)

!       use amrex_fort_module, only : c_real => amrex_real
!       Use usr, only: RK4_POS, RK4_VEL

!       IMPLICIT NONE

! ! Dummy Arguments
! !---------------------------------------------------------------------//
!       real(c_real), intent(in   ) :: ltime
!       integer,      intent(in   ) :: max_pip
!       real(c_real), intent(in   ) :: des_pos_new(max_pip,3)
!       real(c_real), intent(in   ) :: des_vel_new(max_pip,3)

! ! Local variables
! !---------------------------------------------------------------------//
! ! file unit for heat transfer data
!       INTEGER, PARAMETER :: lUNIT = 2030


! ! Open the file.
!       OPEN(UNIT=lUNIT, FILE='POST_POS.dat',                            &
!          POSITION="APPEND", STATUS='OLD')
! ! Write the results to file.
!       WRITE(lUNIT,1000) lTime, RK4_POS(3), DES_POS_NEW(1,3),           &
!                        ABS_ERR(RK4_POS(3), DES_POS_NEW(1,3))
! ! Close the output file.
!       CLOSE(lUNIT)


! ! Open the file.
!       OPEN(UNIT=lUNIT, FILE='POST_VEL.dat',                            &
!          POSITION="APPEND", STATUS='OLD')
! ! Write the results to file.
!       WRITE(lUNIT,1000) lTime, RK4_VEL(3), DES_VEL_NEW(1,3),           &
!                        ABS_ERR(RK4_VEL(3), DES_VEL_NEW(1,3))
! ! Close the output file.
!       CLOSE(lUNIT)

!       RETURN

!  1000 FORMAT(3x,F15.6,5X,F15.6,3x,F15.6,3x,F15.4)

!       CONTAINS


! !......................................................................!
! !                                                                      !
! !  Subroutine name: ABS_ERR                                            !
! !  Author: J.Musser                                   Date:  May-2014  !
! !                                                                      !
! !  Purpose: Calculate either the absolute percent relative error or    !
! !  the absolute error.                                                 !
! !                                                                      !
! !......................................................................!
!       REAL(C_REAL) FUNCTION ABS_ERR(EXT, NUM)

!       use param, only: SMALL_NUMBER

!       IMPLICIT NONE

!       REAL(C_REAL), INTENT(IN) :: EXT, NUM

!       IF(ABS(EXT) > SMALL_NUMBER) THEN
!          ABS_ERR = ABS((EXT - NUM)/EXT)*1.0d2
!       ELSE
!          ABS_ERR = ABS(EXT - NUM)*1.0d2
!       ENDIF
!       END FUNCTION ABS_ERR

!       END SUBROUTINE WRITE_DES_Out
