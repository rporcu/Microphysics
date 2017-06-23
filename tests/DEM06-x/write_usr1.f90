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
   WRITE(lUNIT,1000) lTime, RK4_POS(1), particles(1) % pos(1)
   ! Close the output file.
   CLOSE(lUNIT)


   ! Open the file.
   OPEN(UNIT=lUNIT, FILE='POST_VEL.dat',                            &
        POSITION="APPEND", STATUS='OLD')
   ! Write the results to file.
   WRITE(lUNIT,2000) lTime, RK4_VEL(1), particles(1) % vel(1)
   ! Close the output file.
   CLOSE(lUNIT)

   RETURN

1000 FORMAT(3x,F15.6,5X,F15.6,3x,F15.5)
2000 FORMAT(3x,F15.6,5X,F15.6,3x,F15.4)

end subroutine WRITE_DES_Out
