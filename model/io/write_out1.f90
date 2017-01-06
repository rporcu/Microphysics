MODULE WRITE_OUT1_MODULE

   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_OUT1                                             C
!  Author: P. Nicoletti                               Date: 03-DEC-91  C
!                                                                      C
!  Purpose: write out the field variables to standard output           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_OUT1(time,slo,shi,ep_g,p_g,ro_g)

      USE compar, only: iend3, jend3, kend3
      USE compar, only: mype, pe_io
      USE funits, only: unit_out
      USE out_array_mod, only: out_array
      USE run, only: call_usr

      IMPLICIT NONE

      integer, intent(in   ) :: slo(3),shi(3)

      INTERFACE
         SUBROUTINE USR_WRITE_OUT1
         END SUBROUTINE USR_WRITE_OUT1
      END INTERFACE

      real(c_real), INTENT(IN) :: time

      real(c_real), INTENT(IN) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN) :: p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), INTENT(IN) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

      real(c_real), allocatable :: array1(:,:,:)

      if (myPE == PE_IO) then
         allocate ( array1(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3)) )
      else
         allocate (array1(1,1,1))
      end if
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1000) CHAR(12), TIME
      array1 = P_g
      ! call gather (P_g,array1,root)
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'P_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1100) CHAR(12), TIME
      array1 = EP_g
      ! call gather (EP_g,array1,root)
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'EP_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, 1200) CHAR(12), TIME
      array1 = RO_g
      ! call gather (RO_g,array1,root)
      if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'RO_g')
!
      !if (myPE == PE_IO) WRITE (UNIT_OUT, 1800) CHAR(12), TIME
      !array1 = U_g
      ! call gather (U_g,array1,root)
      !if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'U_g')
!
      !if (myPE == PE_IO) WRITE (UNIT_OUT, 1900) CHAR(12), TIME
      !array1 = V_g
      ! call gather (V_g,array1,root)
      !if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'V_g')
!
      !if (myPE == PE_IO) WRITE (UNIT_OUT, 2000) CHAR(12), TIME
      !array1 = W_g
      ! call gather (W_g,array1,root)
      !if (myPE == PE_IO) CALL OUT_ARRAY (array1, 'W_g')
!
      if (myPE == PE_IO) WRITE (UNIT_OUT, '(/1X,1A1)') CHAR(12)
      IF (CALL_USR) CALL USR_WRITE_OUT1

      deallocate(array1)
!
!             form feed character = CHAR(12)
      if (myPE == PE_IO) WRITE (UNIT_OUT, '(/1X,1A1)') CHAR(12)
      IF (CALL_USR) CALL USR_WRITE_OUT1
      RETURN
 1000 FORMAT(1X,A1,/5X,'--- Gas pressure (P_g) at time ',G12.5,' ---',2/)
 1100 FORMAT(1X,A1,/5X,'--- Void fraction (EP_g) at time ',G12.5,' ---',2/)
 1200 FORMAT(1X,A1,/5X,'--- Gas density (RO_g) at time ',G12.5,' ---',2/)
      END SUBROUTINE WRITE_OUT1
END MODULE WRITE_OUT1_MODULE
