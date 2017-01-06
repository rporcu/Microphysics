module comp_mean_fields_module

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

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
     subroutine comp_mean_fields(slo, shi, lo, hi, &
                                 ep_g, particle_state, des_pos_new, pvol, flag, nparticles, &
                                 dx, dy, dz) &
         bind(C, name="comp_mean_fields")

      use discretelement, only: max_pip
      use param1, only: zero

      use discretelement, only: nonexistent, normal_ghost
      use discretelement, only: entering_ghost, exiting_ghost

      IMPLICIT NONE

      integer, intent(in   ) :: slo(3),shi(3),lo(3),hi(3)

      integer(c_int), intent(in   ) :: nparticles
      integer(c_int), intent(in   ) :: particle_state(nparticles)
      real(c_real), intent(in   ) :: des_pos_new(nparticles,3)
      real(c_real), intent(in   ) :: pvol(nparticles)

      real(c_real), intent(in   ) :: dx, dy, dz

      real(c_real), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters: particles, filter cells, phases
      INTEGER NP
! Fluid cell index
      INTEGER :: I,J,K
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
      DO NP=1,MAX_PIP
         IF(NONEXISTENT==PARTICLE_STATE(NP)) CYCLE
         IF(NORMAL_GHOST==PARTICLE_STATE(NP) .or. &
            ENTERING_GHOST==PARTICLE_STATE(NP) .or. &
            EXITING_GHOST==PARTICLE_STATE(NP)) CYCLE

! Fluid cell containing the particle
         i = floor(des_pos_new(np,1)*odx) + 1
         j = floor(des_pos_new(np,2)*ody) + 1
         k = floor(des_pos_new(np,3)*odz) + 1

! Particle phase for data binning.
! Accumulate total solids volume (by phase)
         SOLVOLINC(I,J,K) = SOLVOLINC(I,J,K) + PVOL(NP)
! Accumulate total solids momenum-ish (by phase)
      ENDDO

! Calculate the cell average solids velocity, the bulk density,
! and the void fraction.
!----------------------------------------------------------------//
      OoVol = 1.0d0/VOL
      DO K = slo(3),shi(3)
         DO J = slo(2),shi(2)
            DO I = slo(1),shi(1)
               IF(flag(i,j,k,1) == 1) then
                  ep_g(i,j,k) = 1.0d0 - solvolinc(i,j,k)*OoVol
               endif
            ENDDO
         ENDDO
      ENDDO

     end subroutine comp_mean_fields
end module comp_mean_fields_module
