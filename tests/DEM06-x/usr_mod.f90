      module usr

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      implicit none

      ! a dummy variable listed in usrnlst.inc
      real(c_real) DUMMY_DP

      ! Initial position and velocity
      real(c_real) :: RK4_POS(3)
      real(c_real) :: RK4_VEL(3)

      real(c_real) :: RK4_TIME = 0.0d0

      ! Fixed gas velocity: (m/sec)
      real(c_real), parameter :: Vg(3) = (/0.4d0, 0.0d0, 0.0d0/)

      contains

!......................................................................!
!                                                                      !
!  Subroutine name: RK4                                                !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Integrate the RK4 solution from RK4_TIME to T1.            !
!                                                                      !
!......................................................................!
      subroutine UPDATE_RK4_SOL(T1)

      implicit none

      real(c_real), INTENT(in) :: T1

      ! Target step size.
      real(c_real), parameter :: RK4_DT_DEFAULT = 1.0d-8

      ! Time step size adjusted for uniform integration.
      real(c_real) :: RK4_DT

      ! Number of integration steps.
      integer(c_int) :: RK4_STEPS

      ! Loop counter.
      integer(c_int) :: LC

      ! Calculate the value for the RK4 solutions.
      RK4_STEPS = max(1,ceiling((T1-RK4_TIME)/RK4_DT_DEFAULT))
      RK4_DT = (T1-RK4_TIME)/dble(RK4_STEPS)

      ! Take the total number of RK4 steps.
      do LC=1, RK4_STEPS
         CALL RK4_V2b3(RK4_DT, RK4_POS, RK4_VEL)
         RK4_TIME = RK4_TIME + RK4_DT
      end do

      end subroutine UPDATE_RK4_SOL

!......................................................................!
!                                                                      !
!  Subroutine name: RK4_V2b3                                           !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Execute one RK4 time step (DT) to calculate the new        !
!  position (POS) and velocity (VEL).                                  !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 22, Equation (67).                                             !
!......................................................................!
      subroutine RK4_V2b3(DT, POS, VEL)

      implicit none

      real(c_real), intent(in) :: DT

      real(c_real), intent(inout) :: POS(3)
      real(c_real), intent(inout) :: VEL(3)

      real(c_real) :: POS_K1(3), POS_K2(3), POS_K3(3), POS_K4(3)
      real(c_real) :: VEL_K1(3), VEL_K2(3), VEL_K3(3), VEL_K4(3)

      POS_K1 = DT * VEL
      VEL_K1 = DT * (Fb() + Fd(VEL))

      POS_K2 = DT * (VEL + VEL_K1/2.0d0)
      VEL_K2 = DT * (Fb() + Fd(VEL + VEL_K1/2.0d0))

      POS_K3 = DT * (VEL + VEL_K2/2.0d0)
      VEL_K3 = DT * (Fb() + Fd(VEL + VEL_K2/2.0d0))

      POS_K4 = DT * (VEL + VEL_K3)
      VEL_K4 = DT * (Fb() + Fd(VEL + VEL_K3))

      POS = POS + (POS_K1 + 2.0d0*POS_K2 + 2.0d0*POS_K3 + POS_K4)/6.0d0
      VEL = VEL + (VEL_K1 + 2.0d0*VEL_K2 + 2.0d0*VEL_K3 + VEL_K4)/6.0d0

      end subroutine RK4_V2b3

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

      use constant, only: D_p0
      use fld_const, only: RO_G0, Mu_g0

      implicit none

      ! Slip velocity magnitude
      real(c_real), intent(in) :: Vs

      Re = (RO_g0 * D_p0(1) * Vs) / Mu_g0

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

!......................................................................!
!                                                                      !
!  Function name: Fb                                                   !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate the upward bouyancy.                             !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 26, Equation (74).                                             !
!......................................................................!
      real(c_real) function Fb( )

      use constant , only: gravity, RO_s0
      use fld_const, only: RO_g0

      implicit none

      real(c_real), dimension(3) :: Fb

      Fb(:) = gravity(:) * (RO_s0(1)-RO_g0)/RO_s0(1)

      end function Fb

!......................................................................!
!                                                                      !
!  Function name: Fd                                                   !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate the drag.                                        !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 26, Equation (74).                                             !
!......................................................................!
      real(c_real) function Fd(Vp)

      use constant, only: RO_s0, D_p0
      use fld_const, only: RO_g0

      implicit none

      real(c_real), intent(in) :: Vp(3) ! particle velocity

      real(c_real) :: Fd(3)
      real(c_real) :: Vpg(3)
      real(c_real) :: Vs

      Vpg(1) = (Vp(1) - Vg(1))**2
      Vpg(2) = (Vp(2) - Vg(2))**2
      Vpg(3) = (Vp(3) - Vg(3))**2

      Vs = dsqrt(dot_product(Vp-Vg,Vp-Vg))

      Fd(:) = (3.0d0*RO_g0)/(4.0d0*RO_s0(1)*D_p0(1)) * Vpg(:)*Cd(Re(Vs))

      end function Fd

      end module usr
