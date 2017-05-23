!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: calc_solids_volume                                      !
!                                                                      !
!  Purpose:                                                            !
!    1) Map particle volume to Eulerian grid.                          !
!    2) Convert the solids volume to the fluid volume fraction.        !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
subroutine calc_volume_fraction(lo, hi,  slo, shi, np, particles, &
                                dx, dy, dz, ep_g, rop_g, ro_g) &
     bind(C, name="calc_volume_fraction")

   use amrex_fort_module, only: c_real => amrex_real
   use iso_c_binding ,    only: c_int
   use particle_mod,      only: particle_t
   use discretelement,    only: normal_particle

   implicit none

   integer(c_int),   intent(in   ) :: lo(3), hi(3), slo(3), shi(3), np
   type(particle_t), intent(in   ) :: particles(np)
   real(c_real),     intent(in   ) :: dx, dy, dz
   real(c_real),     intent(in   ) :: ro_g&
        (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

   real(c_real),     intent(inout) :: ep_g&
        (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
   real(c_real),     intent(inout) :: rop_g&
        (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

   ! Local variables
   !-----------------------------------------------
   ! Fluid cell index
   integer :: i, j, k, n

   ! One over cell width
   real(c_real) :: odx, ody, odz

   ! One over cell volume
   real(c_real) :: oovol

   odx = 1.0d0 / dx
   ody = 1.0d0 / dy
   odz = 1.0d0 / dz

   oovol = odx*ody*odz

   ep_g(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = 1.d0

   ! Calculate the gas phase forces acting on each particle.
   do n = 1, np
     
      if ( particles(n) % state  == normal_particle) then

         ! Fluid cell containing the particle
         i = floor( particles(n) % pos(1) * odx )
         j = floor( particles(n) % pos(2) * ody )
         k = floor( particles(n) % pos(3) * odz)

         ep_g(i,j,k) = ep_g(i,j,k) - oovol * particles(n) % volume

      end if
   end do

   rop_g(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = &
    ro_g(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) * &
    ep_g(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

end subroutine calc_volume_fraction

