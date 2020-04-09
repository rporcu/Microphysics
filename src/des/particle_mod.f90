module particle_mod

  use amrex_fort_module, only: rt => amrex_real
  use iso_c_binding ,    only: c_int

  implicit none
  private

  public  particle_t

  type, bind(C)  :: particle_t
     real(rt)    :: pos(3)     !< Position -- fortran components 1,2,3
     real(rt)    :: radius     !< Radius   -- fortran component  4
     real(rt)    :: volume     !< Volume   -- fortran component  5
     real(rt)    :: mass       !< Mass     -- fortran component  6
     real(rt)    :: density    !< Density  -- fortran component  7
     real(rt)    :: omoi       !< One over momentum of inertia -- fortran component 8
     real(rt)    :: vel(3)     !< Linear velocity              -- fortran components 9,10,11
     real(rt)    :: omega(3)   !< Angular velocity             -- fortran components 12,13,14
     real(rt)    :: drag(3)    !< Drag                         -- fortran components 15,16,17
     integer(c_int)  :: id
     integer(c_int)  :: cpu
     integer(c_int)  :: phase
     integer(c_int)  :: state
  end type particle_t

end module
