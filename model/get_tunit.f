module tunit_module

   CONTAINS

      SUBROUTINE GET_TUNIT(TLEFT, TUNIT)

      implicit none

      DOUBLE PRECISION, INTENT(INOUT) :: TLEFT
      CHARACTER(LEN=4) :: TUNIT

      IF (TLEFT < 3600.0d0) THEN
         TUNIT = 's'
      ELSE
         TLEFT = TLEFT/3600.0d0
         TUNIT = 'h'
         IF (TLEFT >= 24.) THEN
            TLEFT = TLEFT/24.0d0
            TUNIT = 'days'
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE GET_TUNIT

end module tunit_module
