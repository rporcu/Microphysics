!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS1

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE compar
      USE sendrecv
      USE discretelement
      use desgrid
      use desmpi
      USE functions
      use particle_filter, only: FILTER_WEIGHT
      use particle_filter, only: FILTER_CELL
      use physprop, only:MMAX, RO_S0

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters: partciles, filter cells, phases
      INTEGER NP, LC, M
! Fluid cell index
      INTEGER I,J,K,IJK
! Total Mth solids phase volume in IJK
      DOUBLE PRECISION :: SOLVOLINC(DIMENSION_3,MMAX)
! One divided by the total solids volume.
      DOUBLE PRECISION :: OoSOLVOL
! PVOL times statistical weight, and times filter weight
      DOUBLE PRECISION :: VOL_WT, VOLxWEIGHT
! Loop bound for filter
      INTEGER :: LP_BND


!-----------------------------------------------

      SOLVOLINC(:,:) = ZERO

! Loop bounds for interpolation.
      LP_BND = merge(27,9,DO_K)

! Calculate the gas phase forces acting on each particle.
      do NP=1,MAX_PIP
         IF(.NOT.IS_NORMAL(NP) .and. .NOT.IS_GHOST(NP)) CYCLE

         VOL_WT = PVOL(NP)
! Particle phase for data binning.
         M = PIJK(NP,5)

         DO LC=1,LP_BND
            IJK = FILTER_CELL(LC,NP)
! Particle volume times the weight for this cell.
            VOLxWEIGHT = VOL_WT*FILTER_WEIGHT(LC,NP)
! Accumulate total solids volume (by phase)
            SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) + VOLxWEIGHT
         ENDDO
      ENDDO

! Calculate the cell average solids velocity, the bulk density,
! and the void fraction.
!---------------------------------------------------------------------//
        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
         IF(.NOT.fluid_cell(i,j,k)) CYCLE

! calculating the cell average solids velocity for each solids phase
         DO M = 1, MMAX

! calculating the bulk density of solids phase m based on the total
! number of particles having their center in the cell
            DES_ROP_S(IJK,M) = RO_S0(M)*SOLVOLINC(IJK,M)/VOL(IJK)

         ENDDO

      ENDDO
      ENDDO
      ENDDO


! Halo exchange of solids volume fraction data.
      calL SEND_RECV(DES_ROP_S,2)

      end SUBROUTINE COMP_MEAN_FIELDS1
