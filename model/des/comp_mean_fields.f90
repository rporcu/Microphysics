module comp_mean_fields_module

   use comp_mean_fields0_module, only: comp_mean_fields0
   use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3
   use discretelement, only: max_pip
   use geometry, only: vol
   use param1, only: zero
   use particle_filter, only: DES_INTERP_GARG
   use particle_filter, only: DES_INTERP_MEAN_FIELDS
   use particle_filter, only: DES_INTERP_NONE
   use particle_filter, only: DES_INTERP_SCHEME_ENUM
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
     SUBROUTINE COMP_MEAN_FIELDS(ep_g,ro_g,rop_g,pijk,particle_state,&
        particle_phase, pmass, pvol, des_pos_new,des_vel_new, &
        des_radius,des_usr_var,flag,vol_surr,iglobal_id,pinc)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: des_radius
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pmass, pvol
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new, des_usr_var
      INTEGER, DIMENSION(:), INTENT(IN) :: particle_state
      INTEGER, DIMENSION(:), INTENT(IN) :: particle_phase
      INTEGER, DIMENSION(:), INTENT(OUT) :: iglobal_id
      INTEGER, DIMENSION(:,:), INTENT(IN) :: pijk

      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER         , INTENT(IN   ) :: flag&
         (istart3:iend3, jstart3:jend3, kstart3:kend3, 4)
      INTEGER         , INTENT(INOUT) :: pinc&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: vol_surr&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

!......................................................................!

! Calculate field variables from particle data:
      IF(DES_INTERP_SCHEME_ENUM == DES_INTERP_GARG) then
         CALL COMP_MEAN_FIELDS0(ep_g, ro_g, rop_g, particle_phase, pmass,&
            pvol, des_pos_new, des_vel_new, des_radius, des_usr_var, &
            vol_surr, iglobal_id, flag, pinc)
      ELSE
         CALL COMP_MEAN_FIELDS_ZERO_ORDER
      ENDIF

      RETURN

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER

      use discretelement, only: nonexistent, normal_ghost
      use discretelement, only: entering_ghost, exiting_ghost

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters: partciles, filter cells, phases
      INTEGER NP, M
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
      END SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER

      END SUBROUTINE COMP_MEAN_FIELDS

end module comp_mean_fields_module
