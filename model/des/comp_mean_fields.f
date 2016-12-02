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
      SUBROUTINE COMP_MEAN_FIELDS(ep_g,ro_g,rop_g)

      use compar, only:  istart3, iend3, jstart3, jend3, kstart3, kend3

      use particle_filter, only: DES_INTERP_MEAN_FIELDS
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use particle_filter, only: DES_INTERP_GARG

      use comp_mean_fields0_module
      use comp_mean_fields1_module
      use calc_epg_des_module

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

!......................................................................!

! Calculate field variables from particle data:
      IF(DES_INTERP_MEAN_FIELDS) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE(DES_INTERP_NONE) ; CALL COMP_MEAN_FIELDS_ZERO_ORDER
         CASE(DES_INTERP_GARG) ; CALL COMP_MEAN_FIELDS0(ep_g,ro_g,rop_g)
         CASE DEFAULT          ; CALL COMP_MEAN_FIELDS1
         END SELECT
      ELSE
         CALL COMP_MEAN_FIELDS_ZERO_ORDER
      ENDIF

! Calculate the gas phase volume fraction from DES_ROP_s.
      CALL CALC_EPG_DES(ep_g,ro_g,rop_g)

      RETURN
      END SUBROUTINE COMP_MEAN_FIELDS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE discretelement
      use desgrid
      use desmpi
      USE functions
      use physprop, only:MMAX, RO_S0

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters: partciles, filter cells, phases
      INTEGER NP, M
! Fluid cell index
      INTEGER :: I,J,K, IJK
! Total Mth solids phase volume in IJK
      DOUBLE PRECISION :: SOLVOLINC(istart3:iend3, jstart3:jend3, kstart3:kend3, MMAX)
! PVOL times statistical weight
      DOUBLE PRECISION :: VOL_WT

      SOLVOLINC(:,:,:,:) = ZERO

! Calculate the gas phae forces acting on each particle.
      DO NP=1,MAX_PIP
         IF(NONEXISTENT==PARTICLE_STATE(NP)) CYCLE
         IF(NORMAL_GHOST==PARTICLE_STATE(NP) .or. &
            ENTERING_GHOST==PARTICLE_STATE(NP) .or. &
            EXITING_GHOST==PARTICLE_STATE(NP)) CYCLE

         VOL_WT = PVOL(NP)
! Fluid cell containing the particle
         I = PIJK(NP,1)
         J = PIJK(NP,2)
         K = PIJK(NP,3)
         IJK = PIJK(NP,4)
! Particle phase for data binning.
         M = PIJK(NP,5)
! Accumulate total solids volume (by phase)
         SOLVOLINC(I,J,K,M) = SOLVOLINC(I,J,K,M) + VOL_WT
! Accumulate total solids momenum-ish (by phase)
      ENDDO

! Calculate the cell average solids velocity, the bulk density,
! and the void fraction.
!----------------------------------------------------------------//
        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
         IF(.NOT.fluid_at(i,j,k)) CYCLE

! calculating the cell average solids velocity for each solids phase
         DO M = 1, MMAX

! calculating the bulk density of solids phase m based on the total
! number of particles having their center in the cell
            DES_ROP_S(i,j,k,M) = RO_S0(M)*SOLVOLINC(I,J,K,M)/VOL

         ENDDO   ! end loop over M=1,MMAX

      ENDDO
      ENDDO
      ENDDO

! Halo exchange of solids volume fraction data.
      ! CALL SEND_RECV(DES_ROP_S,2)

      RETURN
      END SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER

end module comp_mean_fields_module
