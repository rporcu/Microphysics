module gs_mass_flux_mod
contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!  Purpose:  In the following subroutine the mass flux across a periodic
!            domain with pressure drop is held constant at a
!            user-specified value.  This module is activated only if
!            the user specifies a value for the keyword flux_g in the
!            mfix.dat file.
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   integer(c_int) function goal_seek_mFlux(NIT, gsmf, delp_n, mdot_n, &
      flux_ge, flux_gn, flux_gt, flag)&
      bind(C, name="goal_seek_mFlux")

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc, only: delp_x, delp_y, delp_z, flux_g
      USE param1, only: one
      USE compar   ,only: istart3, iend3, jstart3, jend3, kstart3, kend3, myPE, PE_IO
      USE exit_mod, only: mfix_exit
      USE geometry, only: axy, ayz, axz, cyclic_x_mf, cyclic_y_mf, cyclic_z_mf
      USE utilities, ONLY: mfix_isnan
      USE vavg_mod, ONLY: vavg_flux_g
      use iso_c_binding, only: c_double, c_int
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg

      IMPLICIT NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      integer(c_int), intent(inout) :: nit
      integer(c_int), intent(inout) :: gsmf
      real(c_double), intent(inout) :: delp_n, mdot_n

      real(c_double), intent(inout) :: flux_ge&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: flux_gn&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      real(c_double), intent(inout) :: flux_gt&
         (istart3:iend3,jstart3:jend3,kstart3:kend3)
      integer(c_int), intent(in   ) :: flag&
         (istart3:iend3,jstart3:jend3,kstart3:kend3,4)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER, PARAMETER :: MAXGSMF = 500
      DOUBLE PRECISION, PARAMETER :: omega = 0.9
      DOUBLE PRECISION, PARAMETER :: TOL = 1E-03

      DOUBLE PRECISION :: mdot_nm1, delp_nm1
      DOUBLE PRECISION :: delp_xyz

! Store previous values (only used for GSMF>1)
      mdot_nm1 = mdot_n
      delp_nm1 = delp_n

! Calculate the average gas mass flux and error
      IF(CYCLIC_X_MF)THEN
         delp_n = delp_x
         mdot_n = VAVG_Flux_G(flux_ge, ayz, flag)
      ELSEIF(CYCLIC_Y_MF)THEN
         delp_n = delp_y
         mdot_n = VAVG_Flux_G(flux_gn, axz, flag)
      ELSEIF(CYCLIC_Z_MF)THEN
         delp_n = delp_z
         mdot_n = VAVG_Flux_G(flux_gt, axy, flag)
      ELSE
         RETURN
      ENDIF

      GSMF = gsmf + 1

      IF(GSMF > MAXGSMF) THEN
         IF (myPE.EQ.PE_IO) write(*,5400) MAXGSMF
         CALL mfix_exit(0)
      ENDIF

      IF (mfix_isnan(mdot_n) .OR. mfix_isnan(delp_n)) THEN
         IF (myPE.eq.PE_IO) write(*,*) mdot_n, delp_n, &
            ' NaN being caught in GoalSeekMassFlux '
         stop  32232
         RETURN
      ENDIF


! correct delP
      if(GSMF > 1)then
! Fail-Safe Newton's method
         delp_xyz = delp_n - omega * (delp_n - delp_nm1) * &
            ((mdot_n - FLUX_g)/(mdot_nm1 - FLUX_g)) / &
            ((mdot_n - FLUX_g)/(mdot_nm1 - FLUX_g) - ONE)
      else
         delp_xyz = delp_n*0.99
      endif

! Check for convergence
      IF(abs((mdot_n - FLUX_g)/FLUX_g) >= TOL) THEN
         goal_seek_mFlux = 0
         NIT = 1
      ELSE
         goal_seek_mFlux = 1
         delp_n = delp_xyz

         Write(ERR_MSG,5500) GSMF, delp_n, mdot_n
         call flush_err_msg(header=.false., footer=.false.)
      ENDIF

5500  Format('Mass Flux Iterations:', I0,'   DelP=', &
      G12.5, ' Gas Flux=', G12.5)


      IF(CYCLIC_X_MF)THEN
         delp_x = delp_xyz
      ELSEIF(CYCLIC_Y_MF)THEN
         delp_y = delp_xyz
      ELSEIF(CYCLIC_Z_MF)THEN
         delp_z = delp_xyz
      ENDIF

      RETURN

5400 FORMAT(/1X,70('*')//' From: GoalSeekMassFlux',/&
      ' Message: Number of outer iterations exceeded ', I4,/1X,70('*')/)

   end function goal_seek_mflux
end module gs_mass_flux_mod
