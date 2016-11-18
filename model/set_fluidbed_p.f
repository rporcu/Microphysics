!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_FLUIDBED_P                                          C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Purpose: Set the pressure field inside the bed assuming a fluidized C
!           bed with gravity acting the -ve y-direction                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_FLUIDBED_P

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE constant
      USE discretelement
      USE eos, ONLY: EOSG
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE ic
      USE indices
      USE mpi_utility
      USE param
      USE param1
      USE physprop
      USE scales
      USE sendrecv
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: I, J, K, IJK
! Local loop counter
      INTEGER :: L
! Gas pressure at the axial location j
      DOUBLE PRECISION :: PJ
! Bed weight per unit area
      DOUBLE PRECISION :: BED_WEIGHT
! Total area of a x-z plane
      DOUBLE PRECISION :: AREA
! x-z plane area of one cell
      DOUBLE PRECISION :: dAREA
! Average pressure drop per unit length
      DOUBLE PRECISION :: DPoDX, DPoDY, DPoDZ
!-----------------------------------------------

! If any initial pressures are unspecified skip next section
! calculations.
      DO L = 1, DIMENSION_IC
         IF (IC_DEFINED(L)) THEN
            IF (IC_P_G(L) == UNDEFINED) GOTO 60
            PJ = IC_P_G(L)
         ENDIF
      ENDDO

! Here the pressure in each cell is determined from a specified pressure
! drop across the domain length. This section requires that the pressure
! is already defined in all initial condition regions (otherwise this
! section would be skipped)
! ---------------------------------------------------------------->>>
      IF (DO_I .AND. DELP_X/=UNDEFINED) THEN
         DPODX = DELP_X/XLENGTH
         PJ = PJ - DPODX*HALF*(DX(IMAX1)+DX(IMAX2))
         DO I = IMAX1, IMIN1, -1
            PJ = PJ + DPODX*HALF*(DX(I)+DX(I+1))
            DO K = KMIN1, KMAX1
               DO J = JMIN1, JMAX1
! Bound Checking
                  IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE_PRESSURE(PJ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (DO_J .AND. DELP_Y/=UNDEFINED) THEN
         DPODY = DELP_Y/YLENGTH
         PJ = PJ - DPODY*HALF*(DY(JMAX1)+DY(JMAX2))
         DO J = JMAX1, JMIN1, -1
            PJ = PJ + DPODY*HALF*(DY(J)+DY(J+1))
            DO K = KMIN1, KMAX1
               DO I = IMIN1, IMAX1
! Bound Checking
                  IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE_PRESSURE(PJ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (DO_K .AND. DELP_Z/=UNDEFINED) THEN
         DPODZ = DELP_Z/ZLENGTH
         PJ = PJ - DPODZ*HALF*(DZ(KMAX1)+DZ(KMAX2))
         DO K = KMAX1, KMIN1, -1
            PJ = PJ + DPODZ*HALF*(DZ(K)+DZ(K+1))
            DO J = JMIN1, JMAX1
               DO I = IMIN1, IMAX1
! Bound Checking
                  IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE_PRESSURE(PJ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
! ----------------------------------------------------------------<<<
      GOTO 100   ! pressure in all intial condition region cells was defined

   60 CONTINUE   ! pressure in an initial condition region cell was undefined


! ---------------------------------------------------------------->>>
! Search for an outflow boundary condition where pressure is specified
      PJ = UNDEFINED
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L) .AND. BC_TYPE(L)=='P_OUTFLOW') PJ = BC_P_G(L)
      ENDDO

      IF (PJ == UNDEFINED) THEN
! either a PO was not specified and/or a PO was specified but not the
! pressure at the outlet
         IF (RO_G0 /= UNDEFINED) THEN
! If incompressible flow set P_g to zero
            DO K = kstart3, kend3
            DO J = jstart3, jend3
            DO I = istart3, iend3
               IJK = FUNIJK(i,j,k)
               IF (FLUID_AT(IJK)) P_G(IJK) = ZERO
            ENDDO
            ENDDO
            ENDDO
            GOTO 100

         ELSE   ! compressible case

! Error condition -- no pressure outflow boundary condition is specified
! if a case is compressible and pressure in any of the initial
! conditions regions is unspecified, then a PO is effectively required
! (i.e., is specifies a bc_p_g).
            IF(DMP_LOG)WRITE (UNIT_LOG, 1000)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF


! Set an approximate pressure field assuming that the pressure drop
! balances the weight of the bed, if the initial pressure-field is not
! specified
      DO J = JMAX2, JMIN1, -1

! Find the average weight per unit area over an x-z slice
         BED_WEIGHT = 0.0
         AREA = 0.0
         DO K = KMIN1, KMAX1
            DO I = IMIN1, IMAX1
               IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I,J,K)
               IF (FLUID_AT(IJK)) THEN
                  DAREA = DX(I)*DZ(K)
                  AREA = AREA + DAREA
                  IF (RO_G0 == UNDEFINED) THEN
                     BED_WEIGHT = BED_WEIGHT - DY(J)*GRAVITY_Y*EP_G(IJK)*EOSG(&
                        MW_AVG,PJ,295.15d0)*DAREA
                  ELSE
                     BED_WEIGHT = BED_WEIGHT - DY(J)*GRAVITY_Y*EP_G(IJK)*RO_G0&
                        *DAREA
                  ENDIF
               ENDIF  ! end if (fluid_at(ijk))
            ENDDO    ! end do loop (i=imin1,imax1)
         ENDDO    ! end do loop (k=kmin1,kmax1)

! Global Sum
         call global_all_sum(bed_weight)
         call global_all_sum(area)
         IF (AREA /= 0.0) BED_WEIGHT = BED_WEIGHT/AREA

         PJ = PJ + BED_WEIGHT
         DO K = KMIN1, KMAX1
            DO I = IMIN1, IMAX1
               IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I,J,K)
               IF(FLUID_AT(IJK).AND.P_G(IJK)==UNDEFINED)P_G(IJK)=SCALE_PRESSURE(PJ)
            ENDDO    ! end do (i=imin1,imax1)
         ENDDO   ! end do (k = kmin1,kmax1)
      ENDDO   ! end do (j=jmax2,jimn1, -1)
! end setting an undefined pressure in an initial condition region
! ----------------------------------------------------------------<<<

  100 CONTINUE

      call send_recv(P_G,2)

      RETURN

 1000 FORMAT(/1X,70('*')//' From: SET_FLUIDBED_P'/' Message: Outflow ',&
         'pressure boundary condition (P_OUTFLOW) not found.',/&
         'All the initial pressures (IC_P_g) or at least one P_OUTFLOW',/&
         'condition need to be specified',/1X,70('*')/)

      END SUBROUTINE SET_FLUIDBED_P
