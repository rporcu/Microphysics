!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_INTERP_WEIGHTS                                     !
!                                                                      !
!  Purpose: Calculate weights used for interpolation.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_INTERP_WEIGHTS

      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_DPVM
      use particle_filter, only: DES_INTERP_GAUSS

      SELECT CASE(DES_INTERP_SCHEME_ENUM)
      CASE(DES_INTERP_DPVM);  CALL CALC_INTERP_WEIGHTS1
      CASE(DES_INTERP_GAUSS); CALL CALC_INTERP_WEIGHTS1
      END SELECT

      RETURN
      END SUBROUTINE CALC_INTERP_WEIGHTS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_INTERP_WEIGHTS0                                    !
!                                                                      !
!  Purpose: Calculate weights used for interpolation.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_INTERP_WEIGHTS1

      use discretelement, only: MAX_PIP
      use discretelement, only: PIJK
      use discretelement, only: DES_POS_NEW
      use discretelement, only: XE, YN, ZT
      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_WEIGHT
      use geometry, only: DO_K

      use param1, only: ZERO, ONE
! compar is needed by functions.inc
      use compar

      IMPLICIT NONE

      INTEGER :: L, IDX
      INTEGER :: IJK

      INTEGER :: I, J, K
      INTEGER :: IC, JC, KC
      INTEGER :: Km, Kp

      DOUBLE PRECISION :: WEIGHT
      DOUBLE PRECISION :: WEIGHT_I(-1:1)
      DOUBLE PRECISION :: WEIGHT_J(-1:1)
      DOUBLE PRECISION :: WEIGHT_K(-1:1)

      DO L = 1, MAX_PIP

         IF(IS_NONEXISTENT(L).or.IS_ENTERING(L).or.IS_EXITING(L).or.IS_ENTERING_GHOST(L).or.IS_EXITING_GHOST(L)) CYCLE

         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

! Tentative weights for I indices to the West and East.
         WEIGHT_I(-1) = CALC_FILTER_WEIGHTS(DES_POS_NEW(L,1), XE(I-1))
         WEIGHT_I( 1) = CALC_FILTER_WEIGHTS(XE(I), DES_POS_NEW(L,1))
         WEIGHT_I( 0) = ONE - WEIGHT_I(-1) - WEIGHT_I(1)

! Tentative weights for J indices to the South and North.
         WEIGHT_J(-1) = CALC_FILTER_WEIGHTS(DES_POS_NEW(L,2), YN(J-1))
         WEIGHT_J( 1) = CALC_FILTER_WEIGHTS(YN(J), DES_POS_NEW(L,2))
         WEIGHT_J( 0) = ONE - WEIGHT_J(-1) - WEIGHT_J(1)

! Tentative weights for K indices to the Top and Bottom.
         IF(DO_K) THEN
            Km=-1;  Kp=1
            WEIGHT_K(-1) = CALC_FILTER_WEIGHTS(DES_POS_NEW(L,3),ZT(K-1))
            WEIGHT_K( 1) = CALC_FILTER_WEIGHTS(ZT(K), DES_POS_NEW(L,3))
            WEIGHT_K( 0) = ONE - WEIGHT_K(-1) - WEIGHT_K(1)
         ELSE
            Km= 0; Kp=0
            WEIGHT_K( 0) = ONE
         ENDIF

! Set the default fluid cell index and loop counter.
         IJK = PIJK(L,4)
         IDX=0

! Calculate weights for ghost particles. Only store weights that the
! current process owns.
         DO KC=Km,Kp
         DO IC=-1,+1
         DO JC=-1,+1
            IDX=IDX+1
            WEIGHT = WEIGHT_I(IC)*WEIGHT_J(JC)*WEIGHT_K(KC)
            IF(IS_ON_MYPE_PLUS2LAYERS(I+IC,J+JC,K+KC)) THEN
               IF(fluid_at(i+ic,j+jc,k+kc) .OR. CYCLIC_AT(i+ic,j+jc,k+kc)) THEN
                  FILTER_CELL(IDX,L) = FUNIJK(i+ic,j+jc,k+kc)
                  FILTER_WEIGHT(IDX,L) = WEIGHT
               ELSE
                  FILTER_CELL(IDX,L) = IJK
                  FILTER_WEIGHT(IDX,L) = WEIGHT
               ENDIF
            ELSE
               FILTER_CELL(IDX,L) = IJK
               FILTER_WEIGHT(IDX,L) = ZERO
            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ENDDO

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      DOUBLE PRECISION FUNCTION CALC_FILTER_WEIGHTS(POS1, POS2)

      use particle_filter, only: FILTER_WIDTH_INTERP
      use particle_filter, only: OoFILTER_VOL
      use particle_filter, only: FILTER_WIDTH_INTERPx3

      DOUBLE PRECISION, INTENT(IN) :: POS1, POS2
      DOUBLE PRECISION :: OVERLAP, HEIGHT

      OVERLAP = POS1 - POS2
      IF(OVERLAP < FILTER_WIDTH_INTERP) THEN
         HEIGHT = FILTER_WIDTH_INTERP - OVERLAP
         CALC_FILTER_WEIGHTS = HEIGHT**2 * &
            (FILTER_WIDTH_INTERPx3 - HEIGHT)*OoFILTER_VOL
      ELSE
         CALC_FILTER_WEIGHTS = ZERO
      ENDIF

      END FUNCTION CALC_FILTER_WEIGHTS

      INCLUDE 'functions.inc'

      END SUBROUTINE CALC_INTERP_WEIGHTS1
