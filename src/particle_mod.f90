module particle_mod

  use amrex_fort_module, only: c_real => amrex_real
  use iso_c_binding ,    only: c_int

  implicit none
  private

  public  particle_t
  
  type, bind(C)  :: particle_t
     real(c_real)    :: radius
     real(c_real)    :: volume
     real(c_real)    :: mass
     real(c_real)    :: density
     real(c_real)    :: omoi       !< One over momentum of inertia
     real(c_real)    :: pos(3)     !< Position
     real(c_real)    :: vel(3)     !< Linear velocity
     real(c_real)    :: acc(3)     !< Linear acceleration
     real(c_real)    :: omega(3)   !< Angular velocity
     real(c_real)    :: alpha(3)   !< Angular acceleration
     real(c_real)    :: drag(3)    !< Drag
     integer(c_int)  :: id
     integer(c_int)  :: cpu
     integer(c_int)  :: phase
     integer(c_int)  :: state
  end type particle_t
  



end module
