!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DISPLAY_RESID(NIT, IER)                                !
!  Author: M. Syamlal                                 Date: 8-JUL-96   !
!                                                                      !
!  Purpose: Display residuals                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DISPLAY_RESID(NIT)

      USE residual, only: GROUP_RESID

      IMPLICIT NONE

! iteration number
      INTEGER, INTENT(IN) :: NIT

! Print Location of Max_Resid
!      LOGICAL,PARAMETER:: Print_ijk=.FALSE.


      IF(GROUP_RESID) THEN
         CALL DISPLAY_GROUP_RESID
      ELSE
         CALL DISPLAY_FIELD_RESID
      ENDIF

      RETURN

      contains


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DISPLAY_FIELD_RESID(NIT, IER)                          !
!  Author: M. Syamlal                                 Date: 8-JUL-96   !
!                                                                      !
!  Purpose: Display residuals for each field variable.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DISPLAY_FIELD_RESID

      use param1, only: UNDEFINED_I
      use residual, only: RESID_STRING, RESID_INDEX, RESID

      use error_manager

      IMPLICIT NONE

      INTEGER :: LL, LC, LS, LE

      IF(NIT == 1) THEN
         WRITE(ERR_MSG(1)(1:5),'("  Nit")')
         LC=1
         DO LL = 1, 8
            IF (RESID_INDEX(LL,1) /= UNDEFINED_I) THEN
               LS= 6+10*(LC-1)
               LE= 5+10*(LC)
               WRITE(ERR_MSG(1)(LS:LE),'(5X,A4)') RESID_STRING(LL)
               LC=LC+1
            ENDIF
         END DO
         IF(RESID_INDEX(8,1) == UNDEFINED_I) THEN
            LS= 6+10*(LC-1)
            LE= 5+10*(LC)
            WRITE(ERR_MSG(1)(LS:LE),'(2X,A7)') 'Max res'
         ENDIF
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)
      ENDIF


      WRITE(ERR_MSG(1)(1:5),'(I5)') NIT
      LC=1
      DO LL = 1, 8
         IF(RESID_INDEX(LL,1) /= UNDEFINED_I) THEN
            LS= 6+10*(LC-1)
            LE= 5+10*(LC)
            WRITE(ERR_MSG(1)(LS:LE),'(2X,1PG8.1)') &
               RESID(RESID_INDEX(LL,1))
            LC=LC+1
         ENDIF
      ENDDO
      IF(RESID_INDEX(8,1) == UNDEFINED_I) THEN
         LS= 6+10*(LC-1)
         LE= 3+10*(LC)
         WRITE(ERR_MSG(1)(LS:LE),'(4X,A4)') RESID_STRING(8)
      ENDIF
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)

      RETURN

      END SUBROUTINE DISPLAY_FIELD_RESID



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DISPLAY_GROUP_RESID(NIT, IER)                          !
!  Author: M. Syamlal                                 Date: 8-JUL-96   !
!                                                                      !
!  Purpose: Display residuals grouped by equation type.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DISPLAY_GROUP_RESID

      use residual, only: RESID_STRING
      use residual, only: RESID_GRP, RESID_GRP_STRING
      use residual, only: HYDRO_GRP, THETA_GRP

      use error_manager

      IMPLICIT NONE


      INTEGER :: LC, LS, LE

      IF (NIT == 1) THEN
         WRITE(ERR_MSG(1)(1:5),'("  Nit")')
         LC=1

         LS= 6+10*(LC-1)
         LE= 5+10*(LC)
         WRITE(ERR_MSG(1)(LS:LE),1000) RESID_GRP_STRING(HYDRO_GRP)
         LC=LC+1

         LS= 6+10*(LC-1)
         LE= 5+10*(LC)
         WRITE(ERR_MSG(1)(LS:LE),1000) 'Max res '

         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)
      ENDIF

 1000 FORMAT(3X,A7)



     WRITE(ERR_MSG(1)(1:5),'(I5)') NIT
      LC=1

      LS= 6+10*(LC-1)
      LE= 5+10*(LC)
      WRITE(ERR_MSG(1)(LS:LE),1100) RESID_GRP(HYDRO_GRP)
      LC=LC+1


      LS= 6+10*(LC-1)
      LE= 3+10*(LC)
      WRITE(ERR_MSG(1)(LS:LE),'(4X,A4)') RESID_STRING(8)

 1100 FORMAT(2x,1PG8.1)

      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)

      RETURN
      END SUBROUTINE DISPLAY_GROUP_RESID

      END SUBROUTINE DISPLAY_RESID
