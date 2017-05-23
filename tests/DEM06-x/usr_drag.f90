!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: usr_drag                                               !
!                                                                      !
!  Purpose: Provide a hook for user defined drag law implementation.   !
!                                                                      !
!  This routine is called from inside fluid (TFM) and particle (DES)   !
!  loops. The fluid cell index (IJK) and phase (TFM) or particle index !
!  (DES) is passed.                                                    !
!                                                                      !
!  ***************************   WARNING   **************************  !
!  *----------------------------------------------------------------*  !
!  * The dummy arguments changed in the 2015-1 MFIX Release.        *  !
!  *                                                                *  !
!  *   1) Phase index (M) is now particle index (NP) for DES. This  *  !
!  *      is reflected in the name change M --> M_NP.               *  !
!  *                                                                *  !
!  *   2) The fluid velocity was added as a dummy argument. This    *  !
!  *      provides access to the interpolated gas velocity for      *  !
!  *      coupled DES simulations.                                  *  !
!  *                                                                *  !
!  ******************************************************************  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine drag_usr(i,j,k, M_NP, lDgA, EPg, Mug, ROg, &
                          Vrel, DPM, ROs, lUg, lVg, lWg)

      use amrex_fort_module, only  c_real => amrex_real
      use iso_c_binding    , only: c_int
      use usr              , only: Re, Cd

      implicit none

      ! Index of fluid cell:
      integer(c_int), intent(in) :: i,j,k

      ! DES SOLIDS --> Index of particle (NP); M = particle_phase(NP,5)
      integer(c_int), intent(in) :: M_NP

      ! drag coefficient
      real(c_real), intent(OUT) :: lDgA

      ! gas volume fraction
      real(c_real), intent(in) :: EPg

      ! gas laminar viscosity
      real(c_real), intent(in) :: Mug

      ! gas density
      real(c_real), intent(in) :: ROg

      ! Magnitude of gas-solids relative velocity
      real(c_real), intent(in) :: Vrel

      ! particle diameter of solids phase M or
      ! average particle diameter if PCF
      real(c_real), intent(in) :: DPM

      ! particle density of solids phase M
      real(c_real), intent(in) :: ROs

      ! fluid velocity components:
      ! o TFM: Averaged from faces to cell center
      ! o DES: Interpolated to the particle's position
      real(c_real), intent(in) :: lUg, lVg, lWg

!......................................................................!

      ! Drag force
      lDgA = 0.75d0*(ROg*Vrel/DPM) * Cd(Re(Vrel))

      end subroutine drag_usr
