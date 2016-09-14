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
      SUBROUTINE PHYSICAL_PROP(IER, LEVEL)

      use compar
      use funits
      use geometry
      use indices
      use mpi_utility
      use param1
      use fldvar
      use physprop

      implicit none

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
         if(RO_G0 == UNDEFINED) CALL PHYSICAL_PROP_ROg

! Calculate everything except density. This is called at the start of
! each iteration.
      elseif(LEVEL == 1) then

! Calculate everything. This is invoked via calc_coeff_all as part of
! the initialization (before starting the time march) and at the start
! of each step step thereafter.
      elseif(LEVEL == 2) then
         if(RO_G0 == UNDEFINED) CALL PHYSICAL_PROP_ROg
      endif


! In case of negative density force exit from the physical property
! calculation routine and reduce the time step
      CALL global_all_sum(Err_l, Err_g)
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
      SUBROUTINE PHYSICAL_PROP_ROg


! Global variables:
!-----------------------------------------------------------------------
! Gas phase density (compressible).
      use fldvar, only: RO_g
! Gas phase pressure.
      use fldvar, only: P_g
! Gas phase volume fraction.
      use fldvar, only: EP_g
! Gas phase material density.
      use fldvar, only: ROP_g
! Maximum value for molecular weight (divided by one)
      use toleranc, only: OMW_MAX
! Run time flag for generating negative gas density log files
      use run, only: REPORT_NEG_DENSITY
! Equation of State - GAS
      use eos, only: EOSG

      use functions

      implicit none

! Local Variables:
!-----------------------------------------------------------------------
! Loop indicies
      INTEGER :: IJK   ! Computational cell
! Flag to write log header
      LOGICAL :: wHeader
!......................................................................!

! Initialize:
      wHeader = .TRUE.

      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(WALL_AT(IJK)) cycle IJK_LP

         RO_G(IJK) = EOSG(MW_AVG,P_G(IJK),293.15d0)
         ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)

         IF(RO_G(IJK) < ZERO) THEN
            Err_l(myPE) = 100
            IF(REPORT_NEG_DENSITY)CALL ROgErr_LOG(IJK, wHeader)
         ENDIF
      ENDDO IJK_LP


      RETURN
      END SUBROUTINE PHYSICAL_PROP_ROg



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: NegROg_LOG                                              C
!  Purpose: Record information about the location and conditions that  C
!           resulted in a negative gas phase density.                  C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ROgErr_LOG(IJK, tHeader)

! Simulation time
      use run, only: TIME
! Gas phase density (compressible).
      use fldvar, only: RO_g
! Gas phase pressure.
      use fldvar, only: P_g
      use cutcell

      INTEGER, intent(in) :: IJK
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

      write(lUnit,1001) IJK, I_OF(IJK), J_OF(IJK), K_OF(IJK)
      write(lUnit,"(6x,A,1X,g12.5)",ADVANCE='NO') 'RO_g:', RO_g(IJK)
      write(lUnit,"(2x,A,1X,g12.5)",ADVANCE='NO') 'P_g:', P_g(IJK)
      if(CARTESIAN_GRID) then
         write(lUnit,"(6x,A,1X,L1)",ADVANCE='NO') 'Cut Cell:', CUT_CELL_AT(IJK)
         write(lUnit,"(6x,A,1X,L1)") 'Small Cell:', SMALL_CELL_AT(IJK)
         write(lUnit,"(6x,'Coordinates (E/N/T): ',1X,3(2x, g17.8))") &
            xg_e(I_OF(IJK)), yg_n(J_of(ijk)), zg_t(k_of(ijk))
      endif

      close(lUnit)

      RETURN

 1000 FORMAT(2X,'One or more cells have reported a negative gas',      &
         ' density (RO_g(IJK)). If',/2x,'this is a persistent issue,', &
         ' lower UR_FAC(1) in mfix.dat.')

 1001 FORMAT(/4X,'IJK: ',I8,7X,'I: ',I4,'  J: ',I4,'  K: ',I4)

      END SUBROUTINE ROgErr_LOG



      END SUBROUTINE PHYSICAL_PROP
