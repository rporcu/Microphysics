module comp_mean_fields_module

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
     SUBROUTINE COMP_MEAN_FIELDS(ep_g,ro_g,rop_g,particle_state,&
        particle_phase, pmass, pvol, des_pos_new,des_vel_new, &
        des_radius,des_usr_var,flag,vol_surr)

      use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
      use discretelement, only: max_pip
      use geometry, only: dx, dy, dz, vol
      use param1, only: zero

      use constant, only:MMAX, RO_S0


      use discretelement, only: nonexistent, normal_ghost
      use discretelement, only: entering_ghost, exiting_ghost

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: des_radius(:)
      DOUBLE PRECISION, INTENT(IN) :: pmass(:), pvol(:)
      DOUBLE PRECISION, INTENT(IN) :: des_vel_new(:,:)
      DOUBLE PRECISION, INTENT(IN) :: des_pos_new(:,:)
      DOUBLE PRECISION, INTENT(IN) :: des_usr_var(:,:)
      integer, intent(in) :: particle_state(:)
      integer, intent(in) :: particle_phase(:)


      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER         , INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)
      DOUBLE PRECISION, INTENT(IN   ) :: vol_surr&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

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
      double precision :: OoVol, Oodx, Oody, Oodz

      SOLVOLINC(:,:,:) = ZERO

      Oodx = 1.0d0/dx
      Oody = 1.0d0/dy
      Oodz = 1.0d0/dz

! Calculate the gas phae forces acting on each particle.
      DO NP=1,MAX_PIP
         IF(NONEXISTENT==PARTICLE_STATE(NP)) CYCLE
         IF(NORMAL_GHOST==PARTICLE_STATE(NP) .or. &
            ENTERING_GHOST==PARTICLE_STATE(NP) .or. &
            EXITING_GHOST==PARTICLE_STATE(NP)) CYCLE

! Fluid cell containing the particle
         i = floor(des_pos_new(np,1)*Oodx) + 1
         j = floor(des_pos_new(np,2)*Oody) + 1
         k = floor(des_pos_new(np,3)*Oodz) + 1

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
