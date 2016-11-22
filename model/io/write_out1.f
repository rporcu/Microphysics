!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_OUT1                                             C
!  Author: P. Nicoletti                               Date: 03-DEC-91  C
!                                                                      C
!  Purpose: write out the field variables to standard output           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_OUT1

      USE compar
      USE fldvar
      USE funits
      USE geometry
      USE param
      USE param1
      USE physprop
      USE run

      IMPLICIT NONE

      double precision, allocatable :: array1(:)    !//d


      if (myPE == PE_IO) then
         allocate (array1(ijkmax3))     !//d
      else
         allocate (array1(1))           !//d
      end if
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1000) CHAR(12), TIME
      ! call gather (P_g,array1,root)    !//
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'P_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1100) CHAR(12), TIME
      ! call gather (EP_g,array1,root)    !//
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'EP_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1200) CHAR(12), TIME
      ! call gather (RO_g,array1,root)    !//
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'RO_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1800) CHAR(12), TIME
      ! call gather (U_g,array1,root)    !//
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'U_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1900) CHAR(12), TIME
      ! call gather (V_g,array1,root)    !//
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'V_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 2000) CHAR(12), TIME
      ! call gather (W_g,array1,root)    !//
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'W_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, '(/1X,1A1)') CHAR(12)
      IF (CALL_USR) CALL USR_WRITE_OUT1

      deallocate(array1)  !//
!
!             form feed character = CHAR(12)
      if (myPE == PE_IO) WRITE (UNIT_OUT, '(/1X,1A1)') CHAR(12)
      IF (CALL_USR) CALL USR_WRITE_OUT1
      RETURN
 1000 FORMAT(1X,A1,/5X,'--- Gas pressure (P_g) at time ',G12.5,' ---',2/)
 1100 FORMAT(1X,A1,/5X,'--- Void fraction (EP_g) at time ',G12.5,' ---',2/)
 1200 FORMAT(1X,A1,/5X,'--- Gas density (RO_g) at time ',G12.5,' ---',2/)
 1800 FORMAT(1X,A1,/5X,'--- X-component of gas velocity (U_g) at time ',G12.5,&
         ' ---',2/)
 1900 FORMAT(1X,A1,/5X,'--- Y-component of gas velocity (V_g) at time ',G12.5,&
         ' ---',2/)
 2000 FORMAT(1X,A1,/5X,'--- Z-component of gas velocity (W_g) at time ',G12.5,&
         ' ---',2/)
      END SUBROUTINE WRITE_OUT1
