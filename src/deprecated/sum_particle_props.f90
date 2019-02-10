!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Procedure Name: sum_particle_props                                   !
!                                                                      !
! Purpose: Sum diameter, density and count number of particles. This   !
! is used to calculate the average particle diameter and density which !
! go into calculating the particle collision time scale.               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
subroutine sum_particle_props ( np, particles, sum_np, sum_dp, sum_ro) &
     bind(C, name="sum_particle_props")

  use amrex_fort_module, only: rt => amrex_real
  use iso_c_binding ,    only: c_int

  use param, only: dim_m
  use particle_mod

  integer(c_int),   intent(in   ) :: np
  type(particle_t), intent(in   ) :: particles(np)

  real(rt), intent(inout) :: sum_np(dim_m)
  real(rt), intent(inout) :: sum_dp(dim_m)
  real(rt), intent(inout) :: sum_ro(dim_m)

  integer :: p, m

  do p=1, np
     m = particles(p) % phase
     sum_np(m) = sum_np(m) + 1.0d0
     sum_dp(m) = sum_dp(m) + particles(p) % radius * 2.0d0
     sum_ro(m) = sum_ro(m) + particles(p) % density
  enddo

end subroutine sum_particle_props
