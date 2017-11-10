!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR1 (L)                                         C
!  Purpose: Write user-defined output                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine write_usr1(l, np, time, dt, particles )

   use amrex_fort_module, only: c_real => amrex_real
   use particle_mod,      only: particle_t

   implicit none

   integer,          intent(in) :: l, np
   real(c_real),     intent(in) :: time, dt
   type(particle_t), intent(in) :: particles(np)

   select case(L)
   case(1); call WRITE_DES_OUT(TIME, np, particles)
   end select

end subroutine write_usr1


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

   use usr, only: b_r, h0, time_c, time_r, w0_r, y_s1, dydt_s1, &
        y_s3, dydt_s3, y_s2, dydt_s2

   use constant,          only: gravity
   use amrex_fort_module, only: c_real => amrex_real
   use particle_mod,      only: particle_t

   implicit none

   ! Passed variables
   !---------------------------------------------------------------------//
   integer,          intent(in) :: np
   real(c_real),     intent(in) :: ltime
   type(particle_t), intent(in) :: particles(np)


   ! Local variables
   !---------------------------------------------------------------------//
   ! file unit for heat transfer data
   integer, parameter :: upos = 2030
   integer, parameter :: uvel = 2031

   ! Analytic position and velocity
   real(c_real) :: lPos_Y, lVel_Y
   ! Absolute relative error between MFIX solution and analytic solution.
   real(c_real) :: Pos_aErr, Pos_rErr
   real(c_real) :: Vel_aErr, Vel_rErr

   real(c_real) :: lGrav
   real(c_real) :: lRad
   integer      :: lStage



   ! Open the files.
   OPEN(UNIT=uPOS,FILE='POST_POS.dat',POSITION="APPEND",STATUS='OLD')
   OPEN(UNIT=uVEL,FILE='POST_VEL.dat',POSITION="APPEND",STATUS='OLD')

   ! Set local variables.
   lGrav  = -gravity(2)
   lRad   = 0.1d0
   lStage = 0

   ! Calculate the position and velocity of the particle

   ! Stage 1: Free fall
   if(lTime < time_c) then
      lStage = 1
      lPos_Y = y_s1(h0, lGrav, lTime)
      lVel_Y = dydt_s1(lGrav, lTime)
      ! Stage 2: Contact
   elseif( lTime < time_r) then
      lStage = 2
      lPos_Y = y_s2(h0, lRad, b_r, w0_r, lGrav, lTime)
      lVel_Y = dydt_s2(h0, lRad, b_r, w0_r, lGrav, lTime)
      ! Stage 3: Rebound
   else
      lStage = 3
      lPos_Y = y_s3(lRad, lGrav, lTime)
      lVel_Y = dydt_s3(lGrav, lTime)
   endif

   ! Calculate the absolute and absolute relative errors.
   Pos_aErr = abs(lPos_Y - particles(1) % pos(2) )
   if(lPos_Y > 0.0d0) then
      Pos_rErr = Pos_aErr/lPos_Y*100
   else
      Pos_rErr = Pos_aErr*100
   endif

   Vel_aErr = abs( lVel_Y - particles(1) % vel(2) )
   if(lVel_Y > 0.0d0) then
      Vel_rErr = Vel_aErr/lVel_Y*100
   else
      Vel_rErr = Vel_aErr*100
   endif

   ! Write the results to a file.
   write(uPos,"(F15.8,5x,I1,5X,F15.8,3x,F15.8)") lTime, lStage, lPos_Y, particles(1)%pos(2)
   close(uPos)

   write(uVel,"(F15.8,5x,I1,5X,F15.8,3x,F15.8)") lTime, lStage, lVel_Y, particles(1)%vel(2)
   close(uVel)


end subroutine WRITE_DES_Out
