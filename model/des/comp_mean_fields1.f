module comp_mean_fields1_module

  contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS1

      USE compar, only: iend3, jend3, kend3
      USE compar, only: istart3, jstart3, kstart3
      USE discretelement, only: max_pip, pijk, pvol, des_rop_s
      USE discretelement, only: normal_particle, normal_particle, particle_state, normal_ghost
      USE functions, only: fluid_at, funijk
      USE geometry, only: vol
      USE param1, only: zero
      USE constant, only:MMAX, RO_S0

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters: partciles, filter cells, phases
      INTEGER NP, M
! Fluid cell index
      INTEGER I,J,K,IJK
! Total Mth solids phase volume in IJK
      DOUBLE PRECISION :: SOLVOLINC(istart3:iend3, jstart3:jend3, kstart3:kend3, MMAX)
!-----------------------------------------------

      SOLVOLINC(:,:,:,:) = ZERO
! Calculate the gas phase forces acting on each particle.
      do NP=1,MAX_PIP
         IF(.NOT.NORMAL_PARTICLE==PARTICLE_STATE(NP) .and. .NOT.NORMAL_GHOST==PARTICLE_STATE(NP)) CYCLE

! Particle phase for data binning.
         M = PIJK(NP,5)
! Accumulate total solids volume (by phase)
         SOLVOLINC(I,J,K,M) = SOLVOLINC(I,J,K,M) + PVOL(NP)
      ENDDO

! Calculate the cell average solids velocity, the bulk density,
! and the void fraction.
!---------------------------------------------------------------------//
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

         ENDDO

      ENDDO
      ENDDO
      ENDDO

! Halo exchange of solids volume fraction data.
      ! calL SEND_RECV(DES_ROP_S,2)

      end SUBROUTINE COMP_MEAN_FIELDS1

end module comp_mean_fields1_module
