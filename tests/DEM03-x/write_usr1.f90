!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: write_usr1 (L)                                         C
!  Purpose: Write user-defined output                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
subroutine write_usr1(l, np, time, dt, particles, xlength, ylength, zlength )

   use amrex_fort_module, only: c_real => amrex_real
   use particle_mod,      only: particle_t
   
   implicit none

      integer,          intent(in   ) :: l, np
      real(c_real),     intent(in   ) :: time, dt, xlength, ylength, zlength
      type(particle_t), intent(in   ) :: particles(np)
      
      select case(L)
         case(1); call write_des_out(time, np, particles, xlength)
      end select

      end subroutine write_usr1

!......................................................................!
!  Subroutine: write_des_out                                           !
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

      subroutine write_des_out(lTime, np, particles, length)

      use amrex_fort_module, only: c_real => amrex_real
      use usr,               only: gx1, gx2, gz1, gz2, rk4_v4
      use param,             only: undefined, is_defined
      use particle_mod,      only: particle_t

      implicit none

! Passed variables
!---------------------------------------------------------------------//
      integer,          intent(in   ) ::  np
      real(c_real),     intent(in   ) :: ltime, length
      type(particle_t), intent(in   ) :: particles(np)
      
! Local variables
!---------------------------------------------------------------------//
! file unit for heat transfer data
   integer, parameter :: upos1 = 2030
   integer, parameter :: upos2 = 2031

   real(c_real), save :: rk4_time = 0.0d0
   real(c_real) :: rk4_dt, rk4_dt_last
   real(c_real) time_interval
   real(c_real) rel_diff
   real(c_real), parameter :: rk4_dt_default = 1.0d-6
   integer :: i, rk4_steps

   ! Open the files.
   open(unit=uPos1,FILE='POST_POS1.dat', &
        position="APPEND",STATUS='OLD')

   open(unit=uPos2,FILE='POST_POS2.dat', &
        position="APPEND",STATUS='OLD')

   ! Calculate the value for the RK4 solutions.
   time_interval = lTime - rk4_time

   rk4_dt = rk4_dt_default

   if (time_interval .le. rk4_dt) then
      rk4_steps = 1
      rk4_dt = time_interval
      rk4_dt_last = undefined
   else
      rk4_steps = floor(real(time_interval/rk4_dt_default))
      rk4_dt = rk4_dt_default
      rk4_dt_last = ltime - (rk4_time + rk4_steps*rk4_dt)
   endif

   do i = 1, rk4_steps
      call rk4_v4(rk4_dt, gx1, gz1, gx2, gz2, length)
      rk4_time = rk4_time + rk4_dt
   enddo

   if (is_defined(rk4_dt_last)) then
      call rk4_v4(rk4_dt_last, gx1, gz1, gx2, gz2, length)
      rk4_time = rk4_time + rk4_dt_last
   endif

   rel_diff = (gx1 - particles(1) % pos(1)) / gx1
   rel_diff = dabs(rel_diff) * 100.d0

   ! Write the results to a file.
   write(uPos1,"(3x,F15.8,5X,3(3x,F15.8))") lTime, gx1, particles(1) % pos(1), rel_diff

   rel_diff = (gx2 - particles(2) % pos(1)) / gx2
   rel_diff = dabs(rel_diff) * 100.d0

   write(uPos2,"(3x,F15.8,5X,3(3x,F15.8))") lTime, gx2, particles(2) % pos(1), rel_diff

   close(uPos1)
   close(uPos2)


   end subroutine write_des_out
