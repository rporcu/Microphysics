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

      use amrex_fort_module, only: c_real => amrex_real
      use iso_c_binding    , only: c_int

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

    contains

!......................................................................!
!                                                                      !
!  Function name: Re                                                   !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate the particle Reynods Number.                     !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 22, Equation (67).                                             !
!......................................................................!
      real(c_real) function Re(Vs)

      use fld_const, only: RO_G0, Mu_g0

      implicit none

      ! Slip velocity magnitude
      real(c_real), intent(in) :: Vs

      Re = (RO_g0 * DPM * Vs) / Mu_g0

      end function Re

!......................................................................!
!                                                                      !
!  Function name: Cd                                                   !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate a single particle drag coefficient.              !
!                                                                      !
!  Ref: Schiller and Naumman. (1933) 'A drag coefficient correlation', !
!  Z.Ver. Deutsch Ing., pages 318-320.                                 !
!......................................................................!
      real(c_real) function Cd(Re)

      implicit none

      real(c_real), intent(in) :: Re  ! particle height.

      Cd = 0.0d0
      if (Re > 0.0d0) Cd = (24.0d0/Re)*(1.0d0 + 0.15d0*(Re**0.687d0))

      end function Cd

      end subroutine drag_usr
