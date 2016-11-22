!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_DRAG_W                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_DRAG_U(A_M, B_M, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use particle_filter, only: DES_INTERP_SCHEME_ENUM, DES_INTERP_GARG
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Coefficient at cell corners added to the gas momentum A matix.
      use discretelement, only: DRAG_AM, F_GDS
! Coefficient at cell corners added to gas momentum B vector.
      use discretelement, only: DRAG_BM
! Volume of X-momentum cell
      use geometry, only: VOL
! Flag to calculate Z direction
      use geometry, only: DO_K

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and solids phase..
      use param, only: DIMENSION_3
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
! Flag: Fluid exists at indexed cell
      USE functions, ONLY: funijk
      use functions, only: fluid_cell
! IJK of cell to east.
      use functions, only: ieast
! IJK function for I,J,K that includes mapped indices.
      use compar, only: FUNIJK_MAP_C
! Function for averaging to a scalar cell's east face.
      use fun_avg, only: AVG
! Domain index bounds.
      use compar, only: ISTART2, JSTART2, KSTART2
      use compar, only: IEND2, JEND2, KEND2

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_M(DIMENSION_3,-3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_M(DIMENSION_3)
! Error index
      INTEGER, INTENT(INOUT) :: IER


! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K, IJK, IJMK, IJKM, IJMKM, IJKE
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION :: tmp_A, tmp_B
! Averaging factor
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN


! Average the interpoalted drag force from the cell corners to the cell face.
      IF(DES_INTERP_SCHEME_ENUM == DES_INTERP_GARG)THEN

         AVG_FACTOR = merge(0.25d0, 0.5d0, DO_K)

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)

            IF(.NOT.fluid_cell(i,j,k)) CYCLE

            IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
            IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
            IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

            IJMK = FUNIJK_MAP_C(I, J-1, K)

            tmp_A = -AVG_FACTOR*(DRAG_AM(IJK) + DRAG_AM(IJMK))
            tmp_B = -AVG_FACTOR*(DRAG_BM(IJK,1) + DRAG_BM(IJMK,1))

            IF(DO_K) THEN
               IJKM = FUNIJK_MAP_C(I, J, K-1)
               IJMKM = FUNIJK_MAP_C(I, J-1, K-1)
               tmp_A = tmp_A - AVG_FACTOR*                             &
                  (DRAG_AM(IJKM) + DRAG_AM(IJMKM))
               tmp_B = tmp_B - AVG_FACTOR*                             &
                  (DRAG_BM(IJKM,1) + DRAG_BM(IJMKM,1))
            ENDIF

            A_M(IJK,0) = A_M(IJK,0) + tmp_A*VOL
            B_M(IJK) = B_M(IJK) + tmp_B*VOL

         ENDDO
         ENDDO
         ENDDO

      ELSE

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

            IJK = FUNIJK(i,j,k)
            IF(fluid_cell(i,j,k)) THEN
               IJKE = FUNIJK(ieast(i,j,k),j,k)

               tmp_A = AVG(F_GDS(IJK), F_GDS(IJKE))
               tmp_B = AVG(DRAG_BM(IJK,1), DRAG_BM(IJKE,1))

               A_M(IJK,0) = A_M(IJK,0) - VOL * tmp_A
               B_M(IJK) = B_M(IJK) - VOL * tmp_B
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDIF



      END SUBROUTINE GAS_DRAG_U

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_DRAG_V                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_DRAG_V(A_M, B_M, IER)


! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use particle_filter, only: DES_INTERP_SCHEME_ENUM, DES_INTERP_GARG
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Coefficient at cell corners added to the gas momentum A matix.
      use discretelement, only: DRAG_AM, F_GDS
! Coefficient at cell corners added to gas momentum B vector.
      use discretelement, only: DRAG_BM
! Volume of Y-momentum cell
      use geometry, only: VOL
! Flag to calculate Z direction
      use geometry, only: DO_K

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and solids phase..
      use param, only: DIMENSION_3
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
! Flag: Fluid exists at indexed cell
      USE functions, ONLY: funijk
      use functions, only: fluid_cell
! IJK of cell to north.
      use functions, only: jnorth
! IJK function for I,J,K that includes mapped indices.
      use compar, only: FUNIJK_MAP_C
! Function for averaging to a scalar cell's north face.
      use fun_avg, only: AVG
! Domain index bounds.
      use compar, only: ISTART2, JSTART2, KSTART2
      use compar, only: IEND2, JEND2, KEND2


      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_M(DIMENSION_3, -3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_M(DIMENSION_3)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K, IJK, IMJK, IJKM, IMJKM, IJKN
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION tmp_A, tmp_B
! Averaging factor
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      IF(DES_INTERP_SCHEME_ENUM == DES_INTERP_GARG)THEN

         AVG_FACTOR = merge(0.25d0, 0.5d0, DO_K)

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
            IF(.NOT.fluid_cell(i,j,k)) CYCLE

            IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
            IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
            IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

            IMJK = FUNIJK_MAP_C(I-1,J,K)

            tmp_A = -AVG_FACTOR*(DRAG_AM(IJK) + DRAG_AM(IMJK))
            tmp_B = -AVG_FACTOR*(DRAG_BM(IJK,2) + DRAG_BM(IMJK,2))

            IF(DO_K) THEN

               IJKM = FUNIJK_MAP_C(I,J,K-1)
               IMJKM = FUNIJK_MAP_C(I-1,J,K-1)

               tmp_A = tmp_A - AVG_FACTOR*                             &
                  (DRAG_AM(IJKM) + DRAG_AM(IMJKM))
               tmp_B = tmp_B - AVG_FACTOR*                             &
                  (DRAG_BM(IJKM,2) + DRAG_BM(IMJKM,2))
            ENDIF

            A_M(IJK,0) = A_M(IJK,0) + tmp_A*VOL
            B_M(IJK) = B_M(IJK) + tmp_B*VOL

         ENDDO
         ENDDO
         ENDDO

      ELSE

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
            IF(fluid_cell(i,j,k)) THEN
               IJKN = FUNIJK(i,jnorth(i,j,k),k)

               tmp_A = AVG(F_GDS(IJK), F_GDS(IJKN))
               tmp_B = AVG(DRAG_BM(IJK,2), DRAG_BM(IJKN,2))

               A_M(IJK,0) = A_M(IJK,0) - VOL * tmp_A
               B_M(IJK) = B_M(IJK) - VOL * tmp_B
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE GAS_DRAG_V

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_DRAG_W                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_DRAG_W(A_M, B_M, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use particle_filter, only: DES_INTERP_SCHEME_ENUM, DES_INTERP_GARG
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Coefficient at cell corners added to the gas momentum A matix.
      use discretelement, only: DRAG_AM, F_GDS
! Coefficient at cell corners added to gas momentum B vector.
      use discretelement, only: DRAG_BM
! Volume of Z-momentum cell
      use geometry, only: VOL

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and solids phase..
      use param, only: DIMENSION_3
! Fluid grid loop bounds.
      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
! Flag: Fluid exists at indexed cell
      USE functions, ONLY: funijk
      use functions, only: fluid_cell
! IJK of cell to top.
      use functions, only: ktop
! IJK function for I,J,K that includes mapped indices.
      use compar, only: FUNIJK_MAP_C
! Function for averaging to a scalar cell's north face.
      use fun_avg, only: AVG
! Domain index bounds.
      use compar, only: ISTART2, JSTART2, KSTART2
      use compar, only: IEND2, JEND2, KEND2

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_M(DIMENSION_3,-3:3)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_M(DIMENSION_3)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K, IJK, IMJK, IJMK, IMJMK, IJKT
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION tmp_A, tmp_B
! Averaging factor
! (=0.25 in 3D and =0.5 in 2D)
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      IF(DES_INTERP_SCHEME_ENUM == DES_INTERP_GARG)THEN

         AVG_FACTOR = 0.25d0

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
            IF(.NOT.fluid_cell(i,j,k)) CYCLE

            IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
            IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
            IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

            IMJK = FUNIJK_MAP_C(I-1,J,K)
            IJMK = FUNIJK_MAP_C(I,J-1,K)
            IMJMK = FUNIJK_MAP_C(I-1,J-1,K)

            tmp_A = -AVG_FACTOR*(DRAG_AM(IJK) + DRAG_AM(IMJK) +        &
               DRAG_AM(IJMK) + DRAG_AM(IMJMK))

            tmp_B = -AVG_FACTOR*(DRAG_BM(IJK,3) + DRAG_BM(IMJK,3) +    &
               DRAG_BM(IJMK,3) + DRAG_BM(IMJMK,3))

            A_M(IJK,0) = A_M(IJK,0) + tmp_A*VOL
            B_M(IJK) = B_M(IJK) + tmp_B*VOL

         ENDDO
         ENDDO
         ENDDO

      ELSE

        DO K = kstart3, kend3
        DO J = jstart3, jend3
        DO I = istart3, iend3

         IJK = FUNIJK(i,j,k)
            IF(fluid_cell(i,j,k)) THEN
               IJKT = FUNIJK(i,j,ktop(i,j,k))

               tmp_A = AVG(F_GDS(IJK), F_GDS(IJKT))
               tmp_B = AVG(DRAG_BM(IJK,3), DRAG_BM(IJKT,3))

               A_M(IJK,0) = A_M(IJK,0) - VOL * tmp_A
               B_M(IJK) = B_M(IJK) - VOL * tmp_B
            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ENDIF

      RETURN
      END SUBROUTINE GAS_DRAG_W
