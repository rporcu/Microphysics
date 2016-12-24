module comp_mean_fields_module

   use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
   use discretelement, only: max_pip
   use geometry, only: vol
   use param1, only: zero

   use constant, only:MMAX, RO_S0

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
     SUBROUTINE COMP_MEAN_FIELDS(ep_g, pijk, particle_state, pvol, flag)


      use discretelement, only: nonexistent, normal_ghost
      use discretelement, only: entering_ghost, exiting_ghost

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: pvol(:)
      integer, intent(in) :: particle_state(:)
      integer, intent(in) :: pijk(:,:)

      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER         , INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters: partciles, filter cells, phases
      INTEGER NP
! Fluid cell index
      INTEGER :: I,J,K
! Total Mth solids phase volume in IJK
      DOUBLE PRECISION :: SOLVOLINC&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
! One over cell volume
      double precision :: OoVol

      SOLVOLINC(:,:,:) = ZERO

! Calculate the gas phae forces acting on each particle.
      DO NP=1,MAX_PIP
         IF(NONEXISTENT==PARTICLE_STATE(NP)) CYCLE
         IF(NORMAL_GHOST==PARTICLE_STATE(NP) .or. &
            ENTERING_GHOST==PARTICLE_STATE(NP) .or. &
            EXITING_GHOST==PARTICLE_STATE(NP)) CYCLE

! Fluid cell containing the particle
         I = PIJK(NP,1)
         J = PIJK(NP,2)
         K = PIJK(NP,3)
! Particle phase for data binning.
! Accumulate total solids volume (by phase)
         SOLVOLINC(I,J,K) = SOLVOLINC(I,J,K) + PVOL(NP)
! Accumulate total solids momenum-ish (by phase)
      ENDDO

! Calculate the cell average solids velocity, the bulk density,
! and the void fraction.
!----------------------------------------------------------------//
      OoVol = 1.0d0/VOL
      DO K = kstart3, kend3
         DO J = jstart3, jend3
            DO I = istart3, iend3
               IF(flag(i,j,k,1) == 1) then
                  ep_g(i,j,k) = 1.0d0 - solvolinc(i,j,k)*OoVol
               endif
            ENDDO
         ENDDO
      ENDDO

! Halo exchange of solids volume fraction data.
      ! CALL SEND_RECV(EP_G,2)

      RETURN

      END SUBROUTINE COMP_MEAN_FIELDS

end module comp_mean_fields_module
