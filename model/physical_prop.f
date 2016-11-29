!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP                                           C
!  Purpose: Calculate the indicated physical properties that vary      C
!           with time if directed to do so by the corresponding flag   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Mods for MFIX 2.0 (old name CALC_PHYSPROP)                 C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow SI                                                   C
!  Author: S. Dartevelle                              Date: 01-Jul-02  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!    Perry, R.H., and C.H. Chilton, Chemical Engineer's Handbook, 5th  C
!      edition, McGraw-Hill Kogakusha, Tokyo, 1973.                    C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP(IER, LEVEL, ro_g, p_g, ep_g, rop_g, ro_g0)

      use compar
      use funits
      use geometry
      use param1
      use physprop

      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      implicit none

      double precision, intent(inout) ::  ro_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) ::   p_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) ::  ep_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) :: ro_g0

! Dummy arguments
!-----------------------------------------------------------------------
! Global error Flag.
      INTEGER, intent(inout) :: IER
      INTEGER, intent(in) :: LEVEL

! Local variables
!-----------------------------------------------------------------------
! Arrays for storing errors:
! 100 - Negative gas phase density
! 101 - Negative solids phase density
! 10x - Unclassified
! 900 - Invalid temperature in calc_CpoR
      INTEGER :: Err_l(0:numPEs-1)  ! local
      INTEGER :: Err_g(0:numPEs-1)  ! global

! Initialize error flags.
      Err_l = 0

! Calculate density only. This is invoked several times within iterate,
! making it the most frequently called.
      if(LEVEL == 0) then
         if(RO_G0 == UNDEFINED) CALL PHYSICAL_PROP_ROg(ro_g, p_g, ep_g, rop_g)

! Calculate everything except density. This is called at the start of
! each iteration.
      elseif(LEVEL == 1) then

! Calculate everything. This is invoked via calc_coeff_all as part of
! the initialization (before starting the time march) and at the start
! of each step step thereafter.
      elseif(LEVEL == 2) then
         if(RO_G0 == UNDEFINED) CALL PHYSICAL_PROP_ROg(ro_g, p_g, ep_g, rop_g)
      endif


! In case of negative density force exit from the physical property
! calculation routine and reduce the time step
      err_g = err_l
      ! CALL global_all_sum(Err_l, Err_g)
      IER = maxval(Err_g)
      if(IER == 0) return


! Error handeling. - Local.
!-----------------------------------------------------------------------
! An invalid temperature was found by calc_CpoR. This is a fatal run-
! time error and forces a call to MFIX_EXIT.
      IF(IER == 901 .OR. IER == 902) then
         if(myPE == PE_IO) then
            write(*,2000) IER
            write(UNIT_LOG,2000) IER
         endif
         CALL MFIX_EXIT(myPE)
      ENDIF

      return

 2000 FORMAT(/1X,70('*')/' From: PHYSICAL_PROP',/' Fatal Error 2000:', &
         ' calc_CpoR reporetd an invalid temperature: 0x0', I3/,       &
         'See Cp.log for details. Calling MFIX_EXIT.',/1X,70('*')/)

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_ROg                                       C
!  Purpose: Calculate the gas phase density.                           C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_ROg (ro_g, p_g, ep_g, rop_g)

! Global variables:
!-----------------------------------------------------------------------
      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3

      use fldvar, only: MW_AVG

! Maximum value for molecular weight (divided by one)
      use toleranc, only: OMW_MAX
! Run time flag for generating negative gas density log files
      use run, only: REPORT_NEG_DENSITY
! Equation of State - GAS
      use eos, only: EOSG

      use functions

      implicit none

      ! Gas phase density (compressible):  ro_g
      ! Gas phase pressure              :   p_g
      ! Gas phase volume fraction       :  ep_g
      ! Gas phase material density      : rop_g

      double precision, intent(inout) ::  ro_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) ::   p_g(istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) ::  ep_g(istart3:iend3,jstart3:jend3,kstart3:kend3)

! Local Variables:
!-----------------------------------------------------------------------
! Loop indices
      INTEGER :: I,J,K   ! Computational cell
! Flag to write log header
      LOGICAL :: wHeader
!......................................................................!

! Initialize:
      wHeader = .TRUE.

      DO K = kstart3, kend3
        DO J = jstart3, jend3
          DO I = istart3, iend3
            IF(.NOT.WALL_AT(i,j,k)) THEN

              RO_G(i,j,k) = EOSG(MW_AVG,P_G(I,J,K),293.15d0)
              ROP_G(i,j,k) = RO_G(i,j,k)*EP_G(I,J,K)

              IF(RO_G(I,J,K) < ZERO) THEN
                 Err_l(myPE) = 100
                 IF(REPORT_NEG_DENSITY)CALL ROgErr_LOG(i, j, k, wHeader)
              ENDIF
            ENDIF

          ENDDO
        ENDDO
      ENDDO


      RETURN
      END SUBROUTINE PHYSICAL_PROP_ROg



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ROgErr_LOG
!  Purpose: Record information about the location and conditions that  C
!           resulted in a negative gas phase density.                  C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ROgErr_LOG(i, j, k, tHeader)

! Simulation time
      use run, only: TIME
! Gas phase density (compressible).
      use fldvar, only: RO_g
! Gas phase pressure.
      use fldvar, only: P_g
      use functions, only: funijk

      INTEGER, intent(in) :: i, j, k
      LOGICAL, intent(inout) :: tHeader

      LOGICAL :: lExists
      CHARACTER(LEN=255) :: lFile
      INTEGER, parameter :: lUnit = 4868
      LOGICAL, save :: fHeader = .TRUE.


      lFile = '';
      if(numPEs > 1) then
         write(lFile,"('ROgErr_',I4.4,'.log')") myPE
      else
         write(lFile,"('ROgErr.log')")
      endif
      inquire(file=trim(lFile),exist=lExists)
      if(lExists) then
         open(lUnit,file=trim(adjustl(lFile)),                         &
            status='old', position='append')
      else
         open(lUnit,file=trim(adjustl(lFile)), status='new')
      endif

      if(fHeader) then
         write(lUnit,1000)
         fHeader = .FALSE.
      endif

      if(tHeader) then
         write(lUnit,"(/2x,'Simulation time: ',g12.5)") TIME
         tHeader = .FALSE.
      endif

      write(lUnit,1001)  i, j, k
      write(lUnit,"(6x,A,1X,g12.5)",ADVANCE='NO') 'RO_g:', RO_g(i,j,k)
      write(lUnit,"(2x,A,1X,g12.5)",ADVANCE='NO') 'P_g:', P_g(i,j,k)

      close(lUnit)

      RETURN

 1000 FORMAT(2X,'One or more cells have reported a negative gas',      &
         ' density (RO_g(i,j,k)). If',/2x,'this is a persistent issue,', &
         ' lower UR_FAC(1) in mfix.dat.')

 1001 FORMAT(/4X,'I: ',I4,'  J: ',I4,'  K: ',I4)

      END SUBROUTINE ROgErr_LOG



      END SUBROUTINE PHYSICAL_PROP
