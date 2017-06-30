module des_drag_gp_module

   use amrex_fort_module, only : c_real => amrex_real
   use iso_c_binding , only: c_int

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: des_drag_gp                                             C
!  Purpose: Calculate the gas-particle drag coefficient using          C
!           the gas velocity interpolated to the particle position     C
!           and the particle velocity.                                 C
!           Invoked from des_drag_gs and calc_des_drag_gs              C
!                                                                      C
!  Comments: The BVK drag model and all drag models with the           C
!            polydisperse correction factor (i.e., suffix _PCF)        C
!            require an average particle diameter. This has been       C
!            loosely defined for discrete particles based on their     C
!            solids phase                                              C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine des_drag_gp(slo, shi, NP, particle_vel, fluid_vel, EPg,&
           ro_g, mu_g, f_gp, i,j,k, des_radius,  pvol, ros, particle_phase)

      use drag  , only: drag_syam_obrien, drag_gidaspow, drag_gidaspow_blend,&
         drag_wen_yu, drag_koch_hill, drag_bvk
      use drag, only: syam_obrien, gidaspow, gidaspow_blend, bvk,&
         drag_type_enum, drag_type
      use drag, only: wen_yu, koch_hill, user_drag
      use drag, only: wen_yu_pcf, gidaspow_pcf, gidaspow_blend_pcf, koch_hill_pcf
      use param, only: one

      IMPLICIT NONE

      interface

         subroutine drag_usr(I,J,K, M_NP, lDgA, EPg, Mug, ROg, VREL, DPM, &
            ROs, lUg, lVg, lWg)

            use amrex_fort_module, only : c_real => amrex_real
            use iso_c_binding , only: c_int
            ! Index of fluid cell:
            integer, intent(in) :: I,J,K, M_NP
            real(c_real), intent(OUT) :: lDgA
            real(c_real), intent(in) :: EPg, Mug, ROg, VREL, DPM, ROs, lUg, lVg, lWg
         end subroutine drag_usr
      end interface

      ! indices, associated with current particle
      integer, intent(in   ) :: I, J, K

      real(c_real), intent(in) :: des_radius
      real(c_real), intent(in) :: pvol, ros

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      integer(c_int), intent(in   ) :: slo(3), shi(3)

      ! particle number id.
      integer , intent(in) :: NP

      ! particle velocity
      real(c_real), intent(in) :: particle_vel(3)

      ! fluid velocity interpolated to particle position
      real(c_real), intent(in) :: fluid_vel(3)

      ! Gas phase volume fraction.
      real(c_real), intent(in) :: EPg

      real(c_real), intent(in) :: ro_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(in) :: mu_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(out) :: F_gp

      integer, intent(in) :: particle_phase

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! solids phase index, associated with current particle
      ! integer :: M
! Slip velocity and its magnitude
      real(c_real) :: VSLP(3), VREL
! gas laminar viscosity redefined here to set viscosity at pressure
! boundaries
      real(c_real) :: Mu
! drag coefficient
      real(c_real) :: DgA
! correction factors for implementing polydisperse drag model
! proposed by van der Hoef et al. (2005)
      real(c_real) :: F_cor, tSUM
! average particle diameter in polydisperse systems
      real(c_real) :: DPA
! diameter ratio in polydisperse systems
      real(c_real) :: Y_i
! total solids volume fraction
      real(c_real) :: PHIS
! aliases for void fraction, gas density, gas bulk density,
! solids volume fraction, particle diameter, particle density
      real(c_real) :: ROg, ROPg, DPM
!-----------------------------------------------

! solids phase index of current particle
      ! M = particle_phase
! Gas material and bulk densities
      ROg = RO_G(I,J,K)
      ROPg = RO_G(I,J,K) * EPg
! Laminar viscosity.
      Mu = MU_G(I,J,K)
! Slip velocity and its magnitude
      VSLP = fluid_vel - particle_vel
      VREL = SQRT(dot_product(VSLP, VSLP))
! assign variables for short dummy arguments

      DPM = 2.0d0*des_radius

! Total solids volume fraction.
      PHIS = ONE - EPg

! determine the drag coefficient
      SELECT CASE(DRAG_TYPE_ENUM)

      CASE (SYAM_OBRIEN)
         CALL DRAG_SYAM_OBRIEN(DgA,EPG,Mu,ROg,VREL,DPM)

      CASE (GIDASPOW)
         CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,DPM)

      CASE (GIDASPOW_BLend)
         CALL drag_gidaspow_blend(DgA,EPg,Mu,ROg,ROPg,VREL,DPM)

      CASE (WEN_YU)
         CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,DPM)

      CASE (KOCH_HILL)
         CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,DPM,DPM,phis)

      CASE (USER_DRAG)
         CALL DRAG_USR(I,J,K,NP,DgA,EPg,Mu,ROg,VREL,DPM,ROs, &
            fluid_vel(1), fluid_vel(2), fluid_vel(3))

      CASE (BVK)
         ! calculate the average particle diameter and particle ratio
         ! HACK HACK HACK HACK -- Dependence on rop_s was removed

         CALL DRAG_BVK(DgA,EPg,Mu,ROPg,VREL,DPM,DPM,phis)

      CASE DEFAULT
         WRITE (*, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
         stop 20013
      end select   ! end selection of drag_type

! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
      F_gp = DgA * pvol

      end subroutine des_drag_gp

end module des_drag_gp_module
