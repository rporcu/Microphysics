!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Parse_RESID_string(IER)                                C
!  Author: M. Syamlal                                 Date:12-MAY-97   C
!                                                                      C
!  Purpose: Initialize residuals                                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine parse_resid_string() &
         bind(C, name="parse_resid_string")

      use param1, only: undefined_c, undefined_i
      use residual, only: resid_grp_string, resid_index, max_resid_index, resid_x, resid_prefix
      use residual, only: resid_string, energy_grp, group_resid, ke_grp, nprefix, scalar_grp, hydro_grp, theta_grp
      use write_error_module, only: write_error

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!
!                      local index
      integer          L, L1
!
!                      error message
      CHARACTER(LEN=80)     LINE
!
!                      LOGICAL
      LOGICAL          STRING_DEFINED

!
!-----------------------------------------------
!
!  If user did not define any residual strings use default values
!

      IF(GROUP_RESID) THEN

         RESID_GRP_STRING(HYDRO_GRP)   = 'HYDRO'
         RESID_GRP_STRING(THETA_GRP)   = 'THETA'
         RESID_GRP_STRING(ENERGY_GRP)  = 'ENERGY'
         RESID_GRP_STRING(SCALAR_GRP)  = 'SCALAR'
         RESID_GRP_STRING(KE_GRP)      = 'K-EPS.'

         RESID_STRING = UNDEFINED_C
         RESID_INDEX(8,1) = UNDEFINED_I

         RETURN
      ENDIF

      STRING_DEFINED = .FALSE.
      DO L = 1, MAX_RESID_INDEX
         IF (RESID_STRING(L) /= UNDEFINED_C) STRING_DEFINED = .TRUE.
      END DO
      IF (.NOT.STRING_DEFINED) THEN
         RESID_STRING(1) = 'P0'
         RESID_STRING(2) = 'U0'
         RESID_STRING(3) = 'V0'
         RESID_STRING(4) = 'W0'
      ENDIF


      DO L = 1, MAX_RESID_INDEX
         RESID_INDEX(L,1) = UNDEFINED_I
         STRING_DEFINED = .FALSE.
         DO L1 = 1, NPREFIX

            IF (RESID_STRING(L)(1:1) == RESID_PREFIX(L1)) THEN
               RESID_INDEX(L,1) = L1
               STRING_DEFINED = .TRUE.
               EXIT
            ENDIF
         END DO
         IF (STRING_DEFINED) THEN

            RESID_INDEX(L,2) = ICHAR(RESID_STRING(L)(2:2)) - 48
            !print *,"RESID_INDEX(L,2) = ",RESID_INDEX(L,2)
            IF (RESID_INDEX(L,2)<0) THEN
               WRITE (LINE, '(A, A1, A, A4, A)') 'Error: Phase index ', &
                  RESID_STRING(L)(2:2), ' in RESID_STRING ', RESID_STRING(L), &
                  ' out of bounds.'
               CALL WRITE_ERROR ('PARSE_RESID_STRING', LINE, 1)
            ENDIF

            IF (RESID_STRING(L)(1:1)=='G' .AND. RESID_INDEX(L,2)==0) THEN

               WRITE (LINE, '(A, A1, A, A4, A)') 'Error: Phase index ', &
                  RESID_STRING(L)(2:2), ' in RESID_STRING ', RESID_STRING(L), &
                  ' should be > 0.'
               CALL WRITE_ERROR ('PARSE_RESID_STRING', LINE, 1)
            ENDIF
            IF (RESID_STRING(L)(1:1) == 'X') RESID_INDEX(L,1) = RESID_X + &
               (ICHAR(RESID_STRING(L)(3:3)) - 48) * 10 &
              + ICHAR(RESID_STRING(L)(4:4)) - 48 - 1
         ENDIF
      END DO
      
      end subroutine parse_resid_string
