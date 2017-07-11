module particle_mod

  use amrex_fort_module, only: c_real => amrex_real
  use iso_c_binding ,    only: c_int

  implicit none
  private

  public  particle_t
  public  print_particles

  type, bind(C)  :: particle_t
     real(c_real)    :: pos(3)     !< Position -- fortran components 1,2,3
     real(c_real)    :: radius     !< Radius   -- fortran component  4
     real(c_real)    :: volume     !< Volume   -- fortran component  5
     real(c_real)    :: mass       !< Mass     -- fortran component  6
     real(c_real)    :: density    !< Density  -- fortran component  7
     real(c_real)    :: omoi       !< One over momentum of inertia -- fortran component 8
     real(c_real)    :: vel(3)     !< Linear velocity              -- fortran components 9,10,11
     real(c_real)    :: omega(3)   !< Angular velocity             -- fortran components 12,13,14
     real(c_real)    :: drag(3)    !< Drag                         -- fortran components 15,16,17
     integer(c_int)  :: id
     integer(c_int)  :: cpu
     integer(c_int)  :: phase
     integer(c_int)  :: state
  end type particle_t

contains

   subroutine print_particles( particles )
      type(particle_t), intent(in) :: particles(:)
      integer                      :: p

      do p = 1, size(particles)

         write(*,'(2/,A,I0)')        "Particle ID = ", particles(p) % id
         write(*,'(A,3(es15.6,1X))')  "Position    = ", particles(p) % pos
         write(*,'(A,3(es15.6,1X))')  "Velocity    = ", particles(p) % vel
         write(*,'(A,es15.6)')         "Radius      = ", particles(p) % radius
         write(*,'(A,es15.6)')         "Mass        = ", particles(p) % mass
         write(*,'(A,es15.6)')         "Volume      = ", particles(p) % volume

      end do

   end subroutine print_particles

end module
