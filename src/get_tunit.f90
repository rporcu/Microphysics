module tunit_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   CONTAINS

      SUBROUTINE GET_TUNIT(TLEFT, TUNIT)

      implicit none

      real(c_real), INTENT(inout) :: TLEFT
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
