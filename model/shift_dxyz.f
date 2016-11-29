!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SHIFT_DXYZ                                             C
!  Purpose:  shift the data in the dx,dy,dz arrays from 1:IMAX to      C
!            IMIN1:IMAX1,  1:JMAX to JMIN1:JMAX1 ,                     C
!            1:KMAX to KMIN1:KMAX1                                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 03-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX, IMAX1, IMAX2, JMAX, JMAX1, JMAX2, KMAX  C
!                        KMAX1 , KMAX2, IMIN1, JMIN1, KMIN1, NO_I,     C
!                        NO_J, NO_K                                    C
!  Variables modified:  DX, DY, DZ                                     C
!                                                                      C
!  Local variables: LC                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SHIFT_DXYZ
      USE param
      USE param1
      USE geometry
      IMPLICIT NONE
      RETURN
      END SUBROUTINE SHIFT_DXYZ
