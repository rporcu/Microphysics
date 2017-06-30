module comp_mean_fields_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   implicit none
   private

   public comp_mean_fields
   
contains


   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !                                                                      !
   !  Subroutine: COMP_MEAN_FIELDS                                        !
   !  Author: J.Musser                                   Date: 11-NOV-14  !
   !                                                                      !
   !  Purpose: Driver routine for calculating field variables (DES_ROP_s) !
   !  from particle data.                                                 !
   !                                                                      !
   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   subroutine comp_mean_fields ( slo, shi, np, particles, dx, dy, dz, ep_g) &
        bind(C, name="comp_mean_fields")

      use param,          only: zero
      use discretelement, only: nonexistent, normal_ghost
      use particle_mod,   only: particle_t
      
      implicit none

      integer(c_int),   intent(in   ) :: slo(3), shi(3), np
      type(particle_t), intent(in   ) :: particles(np)
      real(c_real),     intent(in   ) :: dx, dy, dz
      real(c_real),     intent(inout) :: ep_g&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      !-----------------------------------------------
      ! Local variables
      !-----------------------------------------------
      ! Loop counters: particles, filter cells, phases
      integer n
      ! Fluid cell index
      integer :: I,J,K
      ! Total Mth solids phase volume in IJK
      real(c_real) :: SOLVOLINC&
           (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      ! One over cell volume
      real(c_real) :: OoVol, odx, ody, odz, vol


      SOLVOLINC(:,:,:) = ZERO

      odx = 1.0d0/dx
      ody = 1.0d0/dy
      odz = 1.0d0/dz

      vol = dx*dy*dz

      ! Calculate the gas phae forces acting on each particle.
      do n = 1, np
         if ( NONEXISTENT  == particles(n) % state  ) cycle
         if ( NORMAL_GHOST == particles(n) % state) cycle

         ! Fluid cell containing the particle
         i = floor( particles(n) % pos(1) * odx )
         j = floor( particles(n) % pos(2) * ody )
         k = floor( particles(n) % pos(3) * odz )

         ! Particle phase for data binning.
         ! Accumulate total solids volume (by phase)
         SOLVOLINC(I,J,K) = SOLVOLINC(I,J,K) + particles(n) % volume
         ! Accumulate total solids momenum-ish (by phase)
      end do

      ! Calculate the cell average solids velocity, the bulk density,
      ! and the void fraction.
      !----------------------------------------------------------------//
      OoVol = 1.0d0/VOL
      do K = slo(3),shi(3)
         do J = slo(2),shi(2)
            do I = slo(1),shi(1)
               ep_g(i,j,k) = 1.0d0 - solvolinc(i,j,k)*OoVol
            enddo
         enddo
      enddo
      
   end subroutine comp_mean_fields


end module comp_mean_fields_module
