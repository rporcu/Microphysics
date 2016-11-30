!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS1

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
      INTEGER I,J,K,IJK
! Total Mth solids phase volume in IJK
      DOUBLE PRECISION :: SOLVOLINC(DIMENSION_3,MMAX)

!-----------------------------------------------

      SOLVOLINC(:,:) = ZERO
! Calculate the gas phase forces acting on each particle.
      do NP=1,MAX_PIP
         IF(.NOT.IS_NORMAL(NP) .and. .NOT.IS_GHOST(NP)) CYCLE

! Particle phase for data binning.
         M = PIJK(NP,5)
! Accumulate total solids volume (by phase)
         SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) + PVOL(NP)
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
            DES_ROP_S(IJK,M) = RO_S0(M)*SOLVOLINC(IJK,M)/VOL

         ENDDO

      ENDDO
      ENDDO
      ENDDO

! Halo exchange of solids volume fraction data.
      ! calL SEND_RECV(DES_ROP_S,2)

      end SUBROUTINE COMP_MEAN_FIELDS1
