MODULE DEPRECATED_OR_UNKNOWN_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: DEPRECATED_OR_UNKNOWN                               !
!     Author: J.Musser                                Date:  5-SEPT-14 !
!                                                                      !
!     Purpose: This routine is called when a keyword was not matched   !
!     to any of the keywords in the namelist files. This routine       !
!     reports if the keyword was deprecated or incorrect.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEPRECATED_OR_UNKNOWN(LINE_NO, INPUT)


      use error_manager, only: finl_err_msg, flush_err_msg, init_err_msg, ival

      IMPLICIT NONE

      integer, INTENT(IN) :: LINE_NO
      CHARACTER(len=*), INTENT(IN) :: INPUT

! 2015-2 Deprecated list:
!-----------------------------------------------------------------------
!     NAMELIST / DEP_2015_2 / SOME_KEYWORD

! 2015-2 Release Deprecated keywords.
!      STRING=''; STRING = '&DEP_2015_2 '//trim(adjustl(INPUT))//'/'
!      READ(STRING,NML=DEP_2015_1,IOSTAT=IOS)
!      IF(IOS == 0) CALL DEPRECATED(LINE_NO, INPUT, '2015-2')

! Everything else...  This should be the last call in this routine.
      CALL UNKNOWN_KEYWORD(LINE_NO, INPUT)


      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: UNKNOWN_KEYWORD                                     !
!     Author: J.Musser                                Date:  5-SEPT-14 !
!                                                                      !
!     Purpose: Write the error message for deprecated keywords.        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE UNKNOWN_KEYWORD(LINE_NO, INPUT)

      integer, INTENT(IN) :: LINE_NO
      CHARACTER(len=*), INTENT(IN) :: INPUT

      WRITE(*,2000) trim(iVAL(LINE_NO)), trim(INPUT)

      stop 19832

 2000 FORMAT(//1X,70('*')/' From: UNKNOWN_KEYWORD',/' Error 2000: ',   &
         'Unable to process line ',A,' of the mfix.dat file.',2/3x,    &
         A,2/1x,'Possible causes are',/3x,'* Incorrect or illegal ',   &
         'keyword format',/3x,'* Unknown or mistyped name',/3x,'* ',   &
         'The mentioned item is too small (array overflow).', 2/1x,    &
         'Please see the user documentation and update the mfix.dat ', &
         'file. ',/1X,70('*')//)


      END SUBROUTINE UNKNOWN_KEYWORD


      END SUBROUTINE DEPRECATED_OR_UNKNOWN
END MODULE DEPRECATED_OR_UNKNOWN_MODULE
