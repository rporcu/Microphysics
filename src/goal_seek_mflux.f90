module gs_mass_flux_mod

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!  Purpose:  In the following subroutine the mass flux across a periodic
!            domain with pressure drop is held constant at a
!            user-specified value.  This module is activated only if
!            the user specifies a value for the keyword flux_g in the
!            mfix.dat file.
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   integer(c_int) function &
       goal_seek_mflux(ulo, uhi, vlo, vhi, wlo, whi, NIT, gsmf, delp_n, mdot_n, &
                       fluxX, fluxY, fluxZ, dx, dy, dz)&
      bind(C, name="goal_seek_mflux")

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use bc, only: delp_x, delp_y, delp_z, flux_g
      use param1, only: one
      use compar   ,only: myPE, PE_IO
      use exit_mod, only: mfix_exit
      use bc, only: cyclic_x_mf, cyclic_y_mf, cyclic_z_mf
      use utilities, ONLY: mfix_isnan
      use vavg_mod, ONLY: vavg_flux_g
      use error_manager, only: finl_err_msg, err_msg, flush_err_msg

      IMPLICIT NONE

      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

      integer(c_int), intent(inout) :: nit
      integer(c_int), intent(inout) :: gsmf
      real(c_real), intent(inout) :: delp_n, mdot_n

      real(c_real), intent(in   ) :: fluxX&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(in   ) :: fluxY&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(in   ) :: fluxZ&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: dx, dy, dz
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      integer, parameter :: MAXGSMF = 500
      real(c_real), PARAMETER :: omega = 0.9
      real(c_real), PARAMETER :: TOL = 1E-03

      real(c_real) :: mdot_nm1, delp_nm1
      real(c_real) :: delp_xyz
      real(c_real) :: axy, ayz, axz

      axy = dx*dy
      axz = dx*dz
      ayz = dy*dz

! Store previous values (only used for GSMF>1)
      mdot_nm1 = mdot_n
      delp_nm1 = delp_n

! Calculate the average gas mass flux and error
      if (cyclic_x_mf) then
         delp_n = delp_x
         mdot_n = vavg_flux_g(ulo, uhi, fluxX, ayz)
      elseif (cyclic_y_mf) then
         delp_n = delp_y
         mdot_n = vavg_flux_g(vlo, vhi, fluxY, axz)
      elseif (cyclic_z_mf) then
         delp_n = delp_z
         mdot_n = vavg_flux_g(wlo, whi, fluxZ, axy)
      else
         return
      end if

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

      if (cyclic_x_mf) then
         delp_x = delp_xyz
      elseif (cyclic_y_mf) then
         delp_y = delp_xyz
      elseif (cyclic_z_mf) then
         delp_z = delp_xyz
      end if

5400 FORMAT(/1X,70('*')//' From: GoalSeekMassFlux',/&
      ' Message: Number of outer iterations exceeded ', I4,/1X,70('*')/)

   end function goal_seek_mflux
end module gs_mass_flux_mod
