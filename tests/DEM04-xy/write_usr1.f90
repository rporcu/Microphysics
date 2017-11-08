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

   integer,          intent(in   ) :: l, np
   real(c_real),     intent(in   ) :: time, dt
   type(particle_t), intent(in   ) :: particles(np)

   select case(L)
   case(1); call write_des_out(TIME, np, particles)
   end SELECT

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
subroutine write_des_out(ltime, np, particles)

   use discretelement, only: mew, mew_w
   use constant, only: gravity
   use param, only: zero
   use param, only: equal
   use usr, only: u0
   use amrex_fort_module, only : c_real => amrex_real
   use particle_mod,      only: particle_t

   implicit none

   ! Passed variables
   !---------------------------------------------------------------------//
   integer,          intent(in   ) ::  np
   real(c_real),     intent(in   ) :: ltime
   type(particle_t), intent(in   ) :: particles(np)

   ! Local variables
   !---------------------------------------------------------------------//
   ! file unit for heat transfer data
   integer, parameter :: lunit = 2030
   ! Slip velocity at contact, error, non-dimensional values.
   real(c_real) :: slip, err, anl_nd, dem_nd
   ! Flag that rolling friction already ended.
   logical, save :: rollfric_end = .false.

   real(c_real) :: lVel
   real(c_real), parameter :: lRad = 0.00050

   ! Return: Rolling friction already ended.
   if(rollfric_end) return

   lVel = dsqrt(particles(1)%vel(1)**2 + particles(1)%vel(2)**2)

   ! Calculate the slip velocity.
   slip = lvel + particles(1) % omega(3)*lrad

   ! Check for a sign flip or a small difference.
   if(equal(abs(slip),1.0d-6) .or. slip < zero) then
      rollfric_end = .true.

      ! Open the files.
      open(unit=lunit,file='POST_TIME.dat', &
           position="append",status='old')

      ! Calculate the non-dimensional end slip times
      anl_nd = 2.0d0/7.0d0
      dem_nd = abs(mew*9.80665/u0) * ltime

      err = (abs(anl_nd-dem_nd)/abs(anl_nd) )*100.

      ! Write the results to a file.
      write(lunit,1000) mew_w, anl_nd, dem_nd
      close(lunit)

      ! Open the files.
      open(unit=lunit,file='POST_TVEL.dat', &
           position="append",status='old')

      ! Calculate the non-dimensional translational velocity.
      ANL_ND = 5.0d0/7.0d0
      DEM_ND = abs(lVel/u0)

      err = (abs(anl_nd-dem_nd)/abs(anl_nd) )*100.

      ! Write the results to a file.
      write(lunit,1000) mew_w, anl_nd, dem_nd
      close(lunit)

      ! Open the files.
      open(unit=lunit,file='POST_AVEL.dat', &
           position="append",status='old')

      ! Calculate the non-dimensional angular velocity.
      anl_nd = 5.0d0/7.0d0
      dem_nd = abs(particles(1) % omega(3)*lrad/u0)

      err = (abs(anl_nd-dem_nd)/abs(anl_nd) )*100.

      ! Write the results to a file.
      write(lunit,1000) mew_w, anl_nd, dem_nd
      close(lunit)

   endif

   return

1000 FORMAT(3x,F15.8,5X,F15.8,3x,F15.8)

END SUBROUTINE WRITE_DES_Out
