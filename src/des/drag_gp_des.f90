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
                             ro_g, mu_g, f_gp, i,j,k, des_radius,  pvol, particle_phase)

      use compar  , only: myPE
      use exit_mod, only: mfix_exit
      use drag  , only: drag_syam_obrien, drag_gidaspow, drag_gidaspow_blend,&
         drag_wen_yu, drag_koch_hill, drag_bvk
      use drag, only: syam_obrien, gidaspow, gidaspow_blend, bvk,&
         drag_type_enum, drag_type
      use drag, only: wen_yu, koch_hill, user_drag
      use drag, only: wen_yu_pcf, gidaspow_pcf, gidaspow_blend_pcf, koch_hill_pcf
      use param, only: one
      use constant, only: ro_s0
      use constant, only: D_p0

      IMPLICIT NONE

      interface

         subroutine drag_usr(I,J,K, M_NP, lDgA, EPg, Mug, ROg, VREL, DPM, &
            ROs, lUg, lVg, lWg)

            use amrex_fort_module, only : c_real => amrex_real
            use iso_c_binding , only: c_int

            ! Index of fluid cell:
            integer, intent(in) :: I,J,K

            ! TFM SOLIDS --> Index of phase (M)
            ! DES SOLIDS --> Index of particle (NP); M = particle_phase(NP,5)
            integer, intent(in) :: M_NP

            ! drag coefficient
            real(c_real), intent(OUT) :: lDgA

            ! gas volume fraction
            real(c_real), intent(in) :: EPg

            ! gas laminar viscosity
            real(c_real), intent(in) :: Mug

            ! gas density
            real(c_real), intent(in) :: ROg

            ! Magnitude of gas-solids relative velocity
            real(c_real), intent(in) :: VREL

            ! particle diameter of solids phase M or
            ! average particle diameter if PCF
            real(c_real), intent(in) :: DPM

            ! particle density of solids phase M
            real(c_real), intent(in) :: ROs

            ! fluid velocity components:
            ! o TFM: Averaged from faces to cell center
            ! o DES: Interpolated to the particle's position
            real(c_real), intent(in) :: lUg, lVg, lWg
         end subroutine drag_usr
      end interface

      ! indices, associated with current particle
      integer, intent(in   ) :: I, J, K

      real(c_real), intent(in) :: des_radius
      real(c_real), intent(in) :: pvol

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
      integer :: M
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
      M = particle_phase
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
         CALL DRAG_GIDASPOW_BLend(DgA,EPg,Mu,ROg,ROPg,VREL,DPM)

      CASE (WEN_YU)
         CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,DPM)

      CASE (KOCH_HILL)
         CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,DPM,DPM,phis)

      CASE (USER_DRAG)
         CALL DRAG_USR(I,J,K,NP,DgA,EPg,Mu,ROg,VREL,DPM,RO_S0(M), &
            fluid_vel(1), fluid_vel(2), fluid_vel(3))

      CASE DEFAULT

! calculate the average particle diameter and particle ratio
! HACK HACK HACK HACK -- Dependence on rop_s was removed
         tSUM = ONE/D_p0(1)

         DPA = ONE / tSUM
         Y_i = DPM * tSUM

         select CASE(DRAG_TYPE_ENUM)


         CASE (GIDASPOW_PCF)
            CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,DPA)
         CASE (GIDASPOW_BLend_PCF)
            CALL DRAG_GIDASPOW_BLend(DgA,EPg,Mu,ROg,ROPg,VREL,DPA)
         CASE (WEN_YU_PCF)
            CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,DPA)
         CASE (KOCH_HILL_PCF)
            CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,DPM,DPA,phis)
         CASE (BVK)
            CALL DRAG_BVK(DgA,EPg,Mu,ROPg,VREL,DPM,DPA,phis)

         CASE DEFAULT
            WRITE (*, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
            CALL mfix_exit(myPE)
         end SELECT   ! end selection of drag_type


! Modify drag coefficient to account for possible corrections and for
! differences between Model B and Model A; see erratum Beetstra (2007)
         F_cor = Y_i
         DgA = ONE/(Y_i*Y_i) * DgA * F_cor

      end select   ! end selection of drag_type

! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
      F_gp = DgA * pvol

      end subroutine des_drag_gp

end module des_drag_gp_module
