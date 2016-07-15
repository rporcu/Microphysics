!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      !
!  Module name: SET_INDEX1A                                            !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_INDEX1A(I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
         IJKW, IJKE, IJKS, IJKN, IJKB, IJKT)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE geometry
      USE compar
      USE fldvar
      USE indices
      USE boundfunijk
      USE functions

      IMPLICIT NONE

!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, IJKW, IJKE, &
         IJKS, IJKN, IJKB, IJKT
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: TRUE_CORNER
!-----------------------------------------------

      IMJK = UNDEFINED_I
      IPJK = UNDEFINED_I
      IJMK = UNDEFINED_I
      IJPK = UNDEFINED_I
      IJKM = UNDEFINED_I
      IJKP = UNDEFINED_I
      IJKW = UNDEFINED_I
      IJKE = UNDEFINED_I
      IJKS = UNDEFINED_I
      IJKN = UNDEFINED_I
      IJKB = UNDEFINED_I
      IJKT = UNDEFINED_I
      TRUE_CORNER = .FALSE.


      IF(IM1(I).NE.UNDEFINED_I) THEN
        TRUE_CORNER = .FALSE.
        TRUE_CORNER = TRUE_CORNER.OR.I_OF(IJK).EQ.IMIN1
        IF((WALL_AT(IJK).OR.FLOW_AT(IJK)).AND.TRUE_CORNER) THEN
           IMJK = IJK
        ELSE
           IMJK = BOUND_FUNIJK(IM1(I),J,K)
        ENDIF
!
!  IJKW
!
        IF (WALL_AT(IMJK)) THEN
           IJKW = IJK
        ELSE
           IJKW = IMJK
        ENDIF
      ENDIF

      IF(IP1(I).NE.UNDEFINED_I) THEN
        TRUE_CORNER = .FALSE.
        TRUE_CORNER = TRUE_CORNER.OR.I_OF(IJK).EQ.IMAX1
        IF((WALL_AT(IJK).OR.FLOW_AT(IJK)).AND.TRUE_CORNER) THEN
           IPJK = IJK
        ELSE
           IPJK = BOUND_FUNIJK(IP1(I),J,K)
        ENDIF
!
!  IJKE
!
        IF (WALL_AT(IPJK)) THEN
           IJKE = IJK
        ELSE
           IJKE = IPJK
        ENDIF
      ENDIF

      IF(JM1(J).NE.UNDEFINED_I) THEN
        TRUE_CORNER = .FALSE.
        TRUE_CORNER = TRUE_CORNER.OR.J_OF(IJK).EQ.JMIN1
        IF((WALL_AT(IJK).OR.FLOW_AT(IJK)).AND.TRUE_CORNER) THEN
           IJMK = IJK
        ELSE
           IJMK = BOUND_FUNIJK(I,JM1(J),K)
        ENDIF
!
!  IJKS
!
        IF (WALL_AT(IJMK)) THEN
           IJKS = IJK
        ELSE
           IJKS = IJMK
        ENDIF
      ENDIF

      IF(JP1(J).NE.UNDEFINED_I) THEN
        TRUE_CORNER = .FALSE.
        TRUE_CORNER = TRUE_CORNER.OR.J_OF(IJK).EQ.JMAX1
        IF((WALL_AT(IJK).OR.FLOW_AT(IJK)).AND.TRUE_CORNER) THEN
           IJPK = IJK
        ELSE
           IJPK = BOUND_FUNIJK(I,JP1(J),K)
        ENDIF
!
!  IJKN
!
        IF (WALL_AT(IJPK)) THEN
           IJKN = IJK
        ELSE
           IJKN = IJPK
        ENDIF
      ENDIF

      IF(KM1(K).NE.UNDEFINED_I) THEN
        TRUE_CORNER = .FALSE.
        TRUE_CORNER = TRUE_CORNER.OR.K_OF(IJK).EQ.KMIN1
        IF((WALL_AT(IJK).OR.FLOW_AT(IJK)).AND.TRUE_CORNER) THEN
           IJKM = IJK
        ELSE
           IJKM = BOUND_FUNIJK(I,J,KM1(K))
        ENDIF
!
!  IJKB
!
        IF (WALL_AT(IJKM)) THEN
           IJKB = IJK
        ELSE
           IJKB = IJKM
        ENDIF
      ENDIF

      IF(KP1(K).NE.UNDEFINED_I) THEN
        TRUE_CORNER = .FALSE.
        TRUE_CORNER = TRUE_CORNER.OR.K_OF(IJK).EQ.KMAX1
        IF((WALL_AT(IJK).OR.FLOW_AT(IJK)).AND.TRUE_CORNER) THEN
           IJKP = IJK
        ELSE
           IJKP = BOUND_FUNIJK(I,J,KP1(K))
        ENDIF
!
!  IJKT
!
        IF (WALL_AT(IJKP)) THEN
           IJKT = IJK
        ELSE
           IJKT = IJKP
        ENDIF
      ENDIF
!
      RETURN
      END SUBROUTINE SET_INDEX1A
