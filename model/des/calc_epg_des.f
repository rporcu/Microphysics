!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_EPG_DES                                            !
!  Author: R.Garg                                     Date: ??-???-??  !
!                                                                      !
!  Purpose: Calculate the gas phase volume fraction (and in turn the   !
!  gas phase bulk density) from the sum of the solids volume fractions.!
!                                                                      !
!  NOTE: This routine uses a global communication to notify all ranks  !
!  of potential errors. Therefore all ranks can call MFIX_EXIT and     !
!  prevent dead-lock. Communications may be reduced by passing the     !
!  flag back to the caller and combining with other error checks.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_EPG_DES

! Global Variables:
!---------------------------------------------------------------------//
! Flag: Discrete and continuum solids co-exist
      use discretelement, only: DES_CONTINUUM_COUPLED
! Global ID of particles
      use discretelement, only: iGLOBAL_ID
! Particle positions
      use discretelement, only: DES_POS_NEW
! Number of continuum solids phases
      use physprop, only: MMAX, RO_S0
! Discrete particle material and bulk densities
      use discretelement, only: DES_ROP_s
! Number of particles in indexed fluid cell
      use discretelement, only: PINC
! List of particles in each cell.
      use discretelement, only: PIC
! Gas phae volume fraction, density, and build density
      use fldvar, only: EP_G, RO_G, ROP_G
! Volume of scalar grid cell.
      use geometry, only: VOL
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
! Flag: Fluid exists at indexed cell
      use functions, only: fluid_at
      use functions, only: FUNIJK

! Rank ID of current process
      use compar, only: myPE
! Global communication function to sum to all ranks.

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: ZERO, ONE

      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop indices
      INTEGER :: I,j,k,IJK, M, LC
! Total solids volume fraction
      DOUBLE PRECISION SUM_EPS
! Integer Error Flag
      INTEGER :: IER
!......................................................................!

! Initialize error flag.
      IER = 0

! Calculate gas volume fraction from solids volume fraction:
!---------------------------------------------------------------------//
        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
! Skip wall cells.
         IF(.NOT.fluid_at(i,j,k)) CYCLE
! Initialize EP_g and the accumulator.
         EP_G(I,J,K) = ONE
         SUM_EPS = ZERO
! Sum the DES solids volume fraction.
         DO M = 1, MMAX
            SUM_EPS = SUM_EPS + DES_ROP_S(IJK,M)/RO_S0(M)
         ENDDO
! Calculate the gas phase volume fraction and bulk density.
         EP_G(I,J,K) = ONE - SUM_EPS
         ROP_G(i,j,k) = RO_G(I,J,K) * EP_G(I,J,K)
! Flag an error if gas volume fraction is unphysical.
         IF(DES_CONTINUUM_COUPLED) THEN
            IF(EP_G(I,J,K) <= ZERO .OR. EP_G(I,J,K) > ONE) IER = IER + 1
         ENDIF
      ENDDO
      ENDDO
      ENDDO


      ! CALL GLOBAL_ALL_SUM(IER)
      IF(IER == 0) RETURN


! Report any errors. Volume fraction errors are fatal.
!---------------------------------------------------------------------//
      CALL INIT_ERR_MSG("CALC_EPG_DES")
      CALL OPEN_PE_LOG(IER)

      WRITE(ERR_MSG, 1100)
      CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 1100 FORMAT('Error 1100: Unphysical gas phase volume fraction ',      &
         'calculated. A .vtp',/'file will be written and the code ',   &
         'will exit. Fluid cell details:')

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
            IF(.NOT.fluid_at(i,j,k)) CYCLE
            IF(EP_G(I,J,K) > ZERO .AND. EP_G(I,J,K) <= ONE) CYCLE

            WRITE(ERR_MSG,1101) trim(iVal(IJK)), trim(iVal(I)),&
               trim(iVal(J)), trim(iVal(K)),EP_G(I,J,K), &
               trim(iVal(PINC(IJK))), VOL
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

            WRITE(ERR_MSG,1102)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            DO LC=1,PINC(IJK)
               M=PIC(IJK)%P(LC)
               WRITE(ERR_MSG,1103) iGlobal_ID(M), trim(iVal(           &
                  DES_POS_NEW(M,1))), trim(iVal(DES_POS_NEW(M,2))),    &
                  trim(iVal(DES_POS_NEW(M,3)))
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDDO
         ENDDO
         ENDDO
         ENDDO

 1101 FORMAT(/3x,'Fluid Cell IJK: ',A,6x,'I/J/K: (',A,',',A,',',A,')',/&
         T6,'EP_G = ',g11.4,T30,/T6,'PINC: ',A,T30,&
         'VOL = ',g11.4)

 1102 FORMAT(/T6,'Global ID',T30,'Position')

 1103 FORMAT(T6,I9,3x,'(',A,', ',A,', ',A,')')

      WRITE(ERR_MSG, 1104)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
 1104 FORMAT('This is a fatal error. A particle output file (vtp) ',   &
         'will be written',/'to aid debugging.')

      CALL WRITE_DES_DATA
      CALL MFIX_EXIT(myPE)

      END SUBROUTINE CALC_EPG_DES
