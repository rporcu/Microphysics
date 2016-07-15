!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_OUT1                                             C
!  Author: P. Nicoletti                               Date: 03-DEC-91  C
!                                                                      C
!  Purpose: write out the field variables to standard output           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_OUT1

      USE param
      USE param1
      USE physprop
      USE fldvar
      USE run
      USE funits
      USE compar             !//d
      USE mpi_utility        !//d

      IMPLICIT NONE

      INTEGER :: LC, N

      double precision, allocatable :: array1(:)    !//d


      if (myPE == PE_IO) then
         allocate (array1(ijkmax3))     !//d
      else
         allocate (array1(1))           !//d
      end if
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1000) CHAR(12), TIME
      call gather (P_g,array1,root)    !//
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'P_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1100) CHAR(12), TIME
      call gather (EP_g,array1,root)    !//
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'EP_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1200) CHAR(12), TIME
      call gather (RO_g,array1,root)    !//
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'RO_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1800) CHAR(12), TIME
      call gather (U_g,array1,root)    !//
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'U_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1900) CHAR(12), TIME
      call gather (V_g,array1,root)    !//
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'V_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 2000) CHAR(12), TIME
      call gather (W_g,array1,root)    !//
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
 1400 FORMAT(1X,A1,/5X,'--- Solids Phase-',I1,' density x volume',&
         ' fraction (ROP_s) at time ',G12.5,' ---',2/)
 1500 FORMAT(1X,A1,/5X,'--- Gas temperature (T_g) at time ',G12.5,' ---',2/)
 1600 FORMAT(1X,A1,/5X,'--- Solids Phase-',I1,' temperature (T_s)',' at time ',&
         G12.5,' ---',2/)
 1710 FORMAT(1X,A1,/5X,'--- Mass fraction of gas species (X_g) ',&
         I2,' at time ',&
         G12.5,' ---',2/)
 1720 FORMAT(1X,A1,/5X,'--- Mass fraction of solids-',I1,' species (X_s)',I2,&
         ' at time ',G12.5,' ---',2/)
 1800 FORMAT(1X,A1,/5X,'--- X-component of gas velocity (U_g) at time ',G12.5,&
         ' ---',2/)
 1900 FORMAT(1X,A1,/5X,'--- Y-component of gas velocity (V_g) at time ',G12.5,&
         ' ---',2/)
 2000 FORMAT(1X,A1,/5X,'--- Z-component of gas velocity (W_g) at time ',G12.5,&
         ' ---',2/)
 2100 FORMAT(1X,A1,/5X,'--- X-component of Solids Phase-',I1,&
         ' velocity (U_s) at time ',G12.5,' ---',2/)
 2200 FORMAT(1X,A1,/5X,'--- Y-component of Solids Phase-',I1,&
         ' velocity (V_s) at time ',G12.5,' ---',2/)
 2300 FORMAT(1X,A1,/5X,'--- Z-component of Solids Phase-',I1,&
         ' velocity (W_s) at time ',G12.5,' ---',2/)
 2500 FORMAT(1X,A1,/5X,'--- Scalar Field-',I2, ' (Scalar) at time ',G12.5,' ---',2/)
 2600 FORMAT(1X,A1,/5X,'--- Turbulence Field-', ' (K-Epsilon) at time ',G12.5,' ---',2/)
      END SUBROUTINE WRITE_OUT1
