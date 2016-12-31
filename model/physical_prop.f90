MODULE PHYSICAL_PROP_MODULE
   CONTAINS
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
      SUBROUTINE PHYSICAL_PROP(LEVEL, ro_g, p_g, ep_g, rop_g, flag)&
        bind(C, name="physical_prop")

      use funits, only: unit_log
      use param1, only: is_undefined, zero

      USE compar, only: istart3, jstart3, kstart3, iend3, jend3, kend3
      USE compar, only: myPE, PE_IO, numPEs
      USE exit_mod, only: mfix_exit
      use fld_const, only: ro_g0
      use iso_c_binding, only: c_double, c_int

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
! Global error Flag.
      integer(c_int), intent(in   ) :: level

      real(c_double), intent(inout) :: ro_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(  out) :: rop_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(in   ) :: p_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(in   ) :: ep_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)

      integer(c_int), intent(in   ) :: flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

! Local variables
!-----------------------------------------------------------------------
! Arrays for storing errors:
! 100 - Negative gas phase density
! 101 - Negative solids phase density
! 10x - Unclassified
! 900 - Invalid temperature in calc_CpoR
      INTEGER :: Err_l(0:numPEs-1)  ! local
      INTEGER :: Err_g(0:numPEs-1)  ! global
      INTEGER :: IER

! Initialize error flags.
      Err_l = 0

! Calculate density only. This is invoked several times within iterate,
! making it the most frequently called.
      if(LEVEL == 0) then
         if(IS_UNDEFINED(RO_G0)) &
            CALL PHYSICAL_PROP_ROg(ro_g, p_g, ep_g, rop_g, flag)

! Calculate everything except density. This is called at the start of
! each iteration.
      elseif(LEVEL == 1) then

! Calculate everything. This is invoked via calc_coeff_all as part of
! the initialization (before starting the time march) and at the start
! of each step step thereafter.
      elseif(LEVEL == 2) then
         if(IS_UNDEFINED(RO_G0)) &
            CALL PHYSICAL_PROP_ROg(ro_g, p_g, ep_g, rop_g, flag)
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
      !IF(IER == 901 .OR. IER == 902) then
         if(myPE == PE_IO) then
            write(*,2000) IER
            write(UNIT_LOG,2000) IER
         endif
         CALL MFIX_EXIT(myPE)
      !ENDIF

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
      SUBROUTINE PHYSICAL_PROP_ROg (ro_g, p_g, ep_g, rop_g, flag)

! Global variables:
!-----------------------------------------------------------------------
      USE compar   , only: istart3, jstart3, kstart3, iend3, jend3, kend3
      use fld_const, only: mw_avg

! Run time flag for generating negative gas density log files
      use run, only: REPORT_NEG_DENSITY
! Equation of State - GAS
      use eos, only: EOSG


      implicit none

      ! Gas phase density (compressible):  ro_g
      ! Gas phase pressure              :   p_g
      ! Gas phase volume fraction       :  ep_g
      ! Gas phase material density      : rop_g

      double precision, intent(inout) ::  ro_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(inout) :: rop_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) ::   p_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      double precision, intent(in   ) ::  ep_g&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      integer, intent(in   ) ::  flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)

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
            IF(flag(i,j,k,1)<100) THEN

              RO_G(i,j,k) = EOSG(MW_AVG,P_G(I,J,K),293.15d0)
              ROP_G(i,j,k) = RO_G(i,j,k)*EP_G(I,J,K)

              IF(RO_G(I,J,K) < ZERO) THEN
                 Err_l(myPE) = 100
              ENDIF
            ENDIF

          ENDDO
        ENDDO
      ENDDO


      RETURN
      END SUBROUTINE PHYSICAL_PROP_ROg






      END SUBROUTINE PHYSICAL_PROP
END MODULE PHYSICAL_PROP_MODULE
